use std::path::Path;
use std::sync::Arc;

use crate::common::contig::ContigManager;
use crate::db::keys;
use anyhow::Error;
use prost::Message;
use rocksdb::{DBWithThreadMode, MultiThreaded};

#[derive(Debug)]
pub struct DbsnpAnnotator {
    db: DBWithThreadMode<MultiThreaded>,
    contig_manager: Arc<ContigManager>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DbsnpResult {
    pub allele: String,
    pub rs_id: String,
}

impl DbsnpAnnotator {
    pub fn new(db: DBWithThreadMode<MultiThreaded>, contig_manager: Arc<ContigManager>) -> Self {
        Self { db, contig_manager }
    }

    pub(crate) fn from_path(
        path: impl AsRef<Path>,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        tracing::info!("Opening dbSNP database at {}", path.as_ref().display());
        let options = rocksdb::Options::default();
        let db_dbsnp =
            rocksdb::DB::open_cf_for_read_only(&options, &path, ["meta", "dbsnp"], false)?;
        Ok(Self::new(db_dbsnp, contig_manager))
    }

    pub fn annotate_record_dbsnp(&self, key: &[u8]) -> Result<Option<DbsnpResult>, Error> {
        if let Some(raw_value) = self
            .db
            .get_cf(self.db.cf_handle("dbsnp").as_ref().unwrap(), key)?
        {
            let record = crate::pbs::seqvars::DbsnpRecord::decode(&raw_value[..])?;
            Ok(Some(DbsnpResult {
                allele: record.allele,
                rs_id: record.rs_id,
            }))
        } else {
            Ok(None)
        }
    }

    pub fn annotate(&self, vcf_var: &keys::Var) -> anyhow::Result<Option<DbsnpResult>> {
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

        self.annotate_record_dbsnp(&key)
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
                quiet: false,
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
