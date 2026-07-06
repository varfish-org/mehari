use std::path::Path;
use std::sync::Arc;

use crate::common::contig::ContigManager;
use crate::db::keys;
use crate::pbs::seqvars::DbsnpRecord;
use anyhow::Error;
use prost::Message;
use rocksdb::{DBWithThreadMode, MultiThreaded};
use rustc_hash::FxHashMap;

#[derive(Debug)]
pub struct DbsnpAnnotator {
    db: DBWithThreadMode<MultiThreaded>,
    contig_manager: Arc<ContigManager>,
    contig_dict: FxHashMap<String, u32>,
}

impl DbsnpAnnotator {
    pub fn new(
        db: DBWithThreadMode<MultiThreaded>,
        contig_manager: Arc<ContigManager>,
        contig_dict: FxHashMap<String, u32>,
    ) -> Self {
        Self {
            db,
            contig_manager,
            contig_dict,
        }
    }

    pub(crate) fn from_path(
        path: impl AsRef<Path>,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        tracing::info!("Opening dbSNP database at {}", path.as_ref().display());
        let db_dbsnp = crate::db::open_db_for_read(path.as_ref(), "dbsnp")?;

        let contig_dict = {
            let cf_meta = db_dbsnp
                .cf_handle("meta")
                .ok_or_else(|| anyhow::anyhow!("meta CF not found"))?;

            match db_dbsnp.get_cf(&cf_meta, b"contig_dictionary")? {
                Some(bytes) => serde_json::from_slice(&bytes)?,
                None => FxHashMap::default(),
            }
        };

        Ok(Self::new(db_dbsnp, contig_manager, contig_dict))
    }

    pub fn annotate_record_dbsnp(&self, key: &[u8]) -> Result<Option<DbsnpRecord>, Error> {
        if let Some(raw_value) = self
            .db
            .get_cf(self.db.cf_handle("dbsnp").as_ref().unwrap(), key)?
        {
            Ok(Some(DbsnpRecord::decode(&raw_value[..])?))
        } else {
            Ok(None)
        }
    }

    pub fn annotate(&self, vcf_var: &keys::Var) -> anyhow::Result<Option<DbsnpRecord>> {
        // Normalize chrom to primary name to build standard key
        let chrom_std = self
            .contig_manager
            .get_primary_name(&vcf_var.chrom)
            .cloned()
            .unwrap_or_else(|| vcf_var.chrom.clone());

        if let Some(&chrom_id) = self.contig_dict.get(&chrom_std) {
            let key = vcf_var.encode_with_id(chrom_id);
            self.annotate_record_dbsnp(&key)
        } else {
            Ok(None)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::db::CommonPipelineArgs;
    use temp_testdir::TempDir;

    #[test]
    fn test_dbsnp_import_and_annotate() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let input_path = temp.join("dbsnp_test.vcf.gz");
        let output_path = temp.join("dbsnp_db");

        let vcf_content = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10000\trs123\tC\tT\t.\t.\t.
chr1\t10005\trs456;rs789\tA\tG,C\t.\t.\t.";

        crate::db::test_utils::write_indexed_file(&input_path, vcf_content)?;

        let common_args = crate::common::Args {
            verbose: clap_verbosity_flag::Verbosity::new(0, 0),
        };
        let create_args = crate::db::dbsnp::Args {
            common: CommonPipelineArgs {
                assembly: "GRCh38".to_string(),
                input: vec![input_path],
                output: output_path.clone(),
                batch_size: 1000,
                no_progress: false,
                threads: 1,
            },
        };

        crate::db::dbsnp::run(&common_args, &create_args)?;

        let contig_manager = Arc::new(ContigManager::new("GRCh38"));
        let annotator = DbsnpAnnotator::from_path(&output_path, contig_manager)?;

        // Test basic single allele lookup
        let var1 = keys::Var {
            chrom: "chr1".to_string(),
            pos: 10000,
            reference: "C".to_string(),
            alternative: "T".to_string(),
        };
        let res1 = annotator.annotate(&var1)?.unwrap();
        assert_eq!(res1.rs_id, "rs123");
        assert_eq!(res1.allele, "T");

        // Test multi-allelic separation lookup allele 1 (G)
        let var2 = keys::Var {
            chrom: "1".to_string(),
            pos: 10005,
            reference: "A".to_string(),
            alternative: "G".to_string(),
        };
        let res2 = annotator.annotate(&var2)?.unwrap();
        assert_eq!(res2.rs_id, "rs456;rs789");
        assert_eq!(res2.allele, "G");

        // Test multi-allelic separation lookup allele 2 (C)
        let var3 = keys::Var {
            chrom: "1".to_string(),
            pos: 10005,
            reference: "A".to_string(),
            alternative: "C".to_string(),
        };
        let res3 = annotator.annotate(&var3)?.unwrap();
        assert_eq!(res3.rs_id, "rs456;rs789");
        assert_eq!(res3.allele, "C");

        Ok(())
    }
}
