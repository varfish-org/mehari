use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use crate::common::contig::ContigManager;
use crate::db::keys;
use anyhow::{Error, anyhow};
use prost::Message;
use rocksdb::{DBWithThreadMode, MultiThreaded};

#[derive(Debug)]
pub struct CustomDbAnnotator {
    db: DBWithThreadMode<MultiThreaded>,
    contig_manager: Arc<ContigManager>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CustomDbResult {
    pub db_name: String,
    pub fields: HashMap<String, String>,
}

impl CustomDbAnnotator {
    pub fn new(db: DBWithThreadMode<MultiThreaded>, contig_manager: Arc<ContigManager>) -> Self {
        Self { db, contig_manager }
    }

    pub fn get_fields(&self) -> Result<Vec<String>, Error> {
        let cf_meta = self
            .db
            .cf_handle("meta")
            .ok_or_else(|| anyhow!("meta CF not found"))?;
        if let Some(fields_bytes) = self.db.get_cf(&cf_meta, b"fields")? {
            let fields_str = String::from_utf8(fields_bytes)?;
            Ok(fields_str.split(',').map(|s| s.to_string()).collect())
        } else {
            Ok(vec![])
        }
    }

    pub(crate) fn from_path(
        path: impl AsRef<Path>,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        tracing::info!("Opening custom database at {}", path.as_ref().display());
        let options = rocksdb::Options::default();
        let db_custom =
            rocksdb::DB::open_cf_for_read_only(&options, &path, ["meta", "generic"], false)?;
        Ok(Self::new(db_custom, contig_manager))
    }

    pub fn annotate_record_custom(
        &self,
        db_name: &str,
        key: &[u8],
    ) -> Result<Option<CustomDbResult>, Error> {
        if let Some(raw_value) = self
            .db
            .get_cf(self.db.cf_handle("generic").as_ref().unwrap(), key)?
        {
            let record = crate::pbs::seqvars::GenericLookupRecord::decode(&raw_value[..])?;
            Ok(Some(CustomDbResult {
                db_name: db_name.to_string(),
                fields: record.fields,
            }))
        } else {
            Ok(None)
        }
    }

    pub fn annotate(
        &self,
        db_name: &str,
        vcf_var: &keys::Var,
    ) -> anyhow::Result<Option<CustomDbResult>> {
        // Normalize chrom to primary name to build standard key
        let chrom_std = self
            .contig_manager
            .get_primary_name(&vcf_var.chrom)
            .cloned()
            .unwrap_or_else(|| vcf_var.chrom.clone());

        let normalized_var = keys::Var {
            chrom: chrom_std,
            pos: vcf_var.pos,
            reference: vcf_var.reference.clone(),
            alternative: vcf_var.alternative.clone(),
        };

        let key: Vec<u8> = normalized_var.into();
        self.annotate_record_custom(db_name, &key)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use temp_testdir::TempDir;

    #[test]
    fn test_custom_db_import_and_annotate_tsv() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let input_path = temp.join("custom_test.tsv.gz");
        let output_path = temp.join("custom_db");

        // Write a mock TSV file
        let tsv_content = "Chrom\tPosition\tRef\tAlt\tScore1\tScore2\
1\t10000\tC\tT\tA1\tB1\
chr1\t10005\tA\tG\tA2\tB2";
        crate::db::test_utils::write_indexed_file(&input_path, tsv_content)?;

        let common_args = crate::common::Args {
            verbose: clap_verbosity_flag::Verbosity::new(0, 0),
        };
        let create_args = crate::db::generic::Args {
            assembly: "GRCh38".to_string(),
            input: vec![input_path],
            output: output_path.clone(),
            db_name: "custom_db".to_string(),
            format: "tsv".to_string(),
            col_chrom: "Chrom".to_string(),
            col_pos: "Position".to_string(),
            col_ref: "Ref".to_string(),
            col_alt: "Alt".to_string(),
            col_values: Some(vec!["Score1".to_string(), "Score2".to_string()]),
            vcf_info_fields: None,
            batch_size: 1000,
        };

        crate::db::generic::run(&common_args, &create_args)?;

        let contig_manager = Arc::new(ContigManager::new("GRCh38"));
        let annotator = CustomDbAnnotator::from_path(&output_path, contig_manager)?;

        // Check fields read from meta CF
        let fields = annotator.get_fields()?;
        assert_eq!(fields, vec!["Score1", "Score2"]);

        // Query 1
        let var1 = keys::Var {
            chrom: "chr1".to_string(),
            pos: 10000,
            reference: "C".to_string(),
            alternative: "T".to_string(),
        };
        let res1 = annotator.annotate("custom_db", &var1)?.unwrap();
        assert_eq!(res1.db_name, "custom_db");
        assert_eq!(res1.fields.get("Score1").unwrap(), "A1");
        assert_eq!(res1.fields.get("Score2").unwrap(), "B1");

        // Query 2
        let var2 = keys::Var {
            chrom: "1".to_string(),
            pos: 10005,
            reference: "A".to_string(),
            alternative: "G".to_string(),
        };
        let res2 = annotator.annotate("custom_db", &var2)?.unwrap();
        assert_eq!(res2.fields.get("Score1").unwrap(), "A2");
        assert_eq!(res2.fields.get("Score2").unwrap(), "B2");

        Ok(())
    }
}
