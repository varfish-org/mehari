use std::path::Path;
use std::sync::Arc;

use crate::common::contig::ContigManager;
use crate::db::keys;
use anyhow::Error;
use prost::Message;
use rocksdb::{DBWithThreadMode, MultiThreaded};

#[derive(Debug)]
pub struct CaddAnnotator {
    db: DBWithThreadMode<MultiThreaded>,
    contig_manager: Arc<ContigManager>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CaddResult {
    pub raw_score: f32,
    pub phred: f32,
}

impl CaddAnnotator {
    pub fn new(db: DBWithThreadMode<MultiThreaded>, contig_manager: Arc<ContigManager>) -> Self {
        Self { db, contig_manager }
    }

    pub(crate) fn from_path(
        path: impl AsRef<Path>,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        tracing::info!("Opening CADD database at {}", path.as_ref().display());
        let options = rocksdb::Options::default();
        let db_cadd = rocksdb::DB::open_cf_for_read_only(&options, &path, ["meta", "cadd"], false)?;
        Ok(Self::new(db_cadd, contig_manager))
    }

    pub fn annotate_record_cadd(&self, key: &[u8]) -> Result<Option<CaddResult>, Error> {
        if let Some(raw_value) = self
            .db
            .get_cf(self.db.cf_handle("cadd").as_ref().unwrap(), key)?
        {
            let record = crate::pbs::seqvars::CaddRecord::decode(&raw_value[..])?;
            Ok(Some(CaddResult {
                raw_score: record.raw_score,
                phred: record.phred,
            }))
        } else {
            Ok(None)
        }
    }

    pub fn annotate(&self, vcf_var: &keys::Var) -> anyhow::Result<Option<CaddResult>> {
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
        self.annotate_record_cadd(&key)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use temp_testdir::TempDir;

    #[test]
    fn test_cadd_import_and_annotate() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let input_path = temp.join("cadd_test.tsv.gz");
        let output_path = temp.join("cadd_db");

        let tsv_content = "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\
1\t10000\tC\tT\t1.5\t15.0\
chr1\t10005\tA\tG\t2.3\t20.5";
        crate::db::test_utils::write_indexed_file(&input_path, tsv_content)?;

        let common_args = crate::common::Args {
            verbose: clap_verbosity_flag::Verbosity::new(0, 0),
        };
        let create_args = crate::db::cadd::Args {
            assembly: "GRCh38".to_string(),
            input: vec![input_path],
            output: output_path.clone(),
            batch_size: 1000,
        };

        crate::db::cadd::run(&common_args, &create_args)?;

        let contig_manager = Arc::new(ContigManager::new("GRCh38"));
        let annotator = CaddAnnotator::from_path(&output_path, contig_manager)?;

        let var1 = keys::Var {
            chrom: "chr1".to_string(),
            pos: 10000,
            reference: "C".to_string(),
            alternative: "T".to_string(),
        };
        let res1 = annotator.annotate(&var1)?.unwrap();
        assert_eq!(res1.raw_score, 1.5);
        assert_eq!(res1.phred, 15.0);

        let var2 = keys::Var {
            chrom: "1".to_string(),
            pos: 10005,
            reference: "A".to_string(),
            alternative: "G".to_string(),
        };
        let res2 = annotator.annotate(&var2)?.unwrap();
        assert_eq!(res2.raw_score, 2.3);
        assert_eq!(res2.phred, 20.5);

        let var3 = keys::Var {
            chrom: "1".to_string(),
            pos: 99999,
            reference: "C".to_string(),
            alternative: "G".to_string(),
        };
        assert!(annotator.annotate(&var3)?.is_none());

        Ok(())
    }
}
