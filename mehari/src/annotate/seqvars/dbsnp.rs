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

    /// Open a dbSNP database built with `annonars dbsnp import`.
    pub(crate) fn from_path(
        path: impl AsRef<Path>,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        tracing::info!("Opening dbSNP database at {}", path.as_ref().display());
        let options = rocksdb::Options::default();
        let db_dbsnp =
            rocksdb::DB::open_cf_for_read_only(&options, &path, ["meta", "dbsnp_data"], false)?;

        Ok(Self::new(db_dbsnp, contig_manager))
    }

    pub fn annotate_record_dbsnp(&self, key: &[u8]) -> Result<Option<DbsnpResult>, Error> {
        if let Some(raw_value) = self
            .db
            .get_cf(self.db.cf_handle("dbsnp_data").as_ref().unwrap(), key)?
        {
            let record = annonars::dbsnp::pbs::Record::decode(&raw_value[..])?;
            Ok(Some(DbsnpResult {
                allele: record.alt_allele,
                rs_id: format!("rs{}", record.rs_id),
            }))
        } else {
            Ok(None)
        }
    }

    pub fn annotate(&self, vcf_var: &keys::Var) -> anyhow::Result<Option<DbsnpResult>> {
        // Only attempt lookups into RocksDB for canonical contigs.
        if let Some(chrom) = self.contig_manager.get_annonars_chrom(&vcf_var.chrom) {
            let key: Vec<u8> = annonars::common::keys::Var::from(
                chrom,
                vcf_var.pos,
                &vcf_var.reference,
                &vcf_var.alternative,
            )
            .into();
            return self.annotate_record_dbsnp(&key);
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use temp_testdir::TempDir;

    #[test]
    fn test_dbsnp_import_and_annotate() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let input_path = temp.join("dbsnp_test.vcf.gz");
        let output_path = temp.join("dbsnp_db");

        let vcf_content = "##fileformat=VCFv4.2
##reference=GRCh38.p14
##dbSNP_BUILD_ID=156
##INFO=<ID=RS,Number=1,Type=Integer,Description=\"dbSNP ID\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10000\trs123\tC\tT\t.\t.\tRS=123
1\t10005\trs456\tA\tG,C\t.\t.\tRS=456\n";

        let mut writer = noodles_bgzf::io::Writer::new(File::create(&input_path)?);
        writer.write_all(vcf_content.as_bytes())?;
        writer.finish()?;
        let index = noodles::vcf::fs::index(&input_path)?;
        noodles::tabix::fs::write(input_path.with_added_extension("tbi"), &index)?;

        let common_args = annonars::common::cli::Args {
            verbose: Default::default(),
        };
        let import_args = annonars::dbsnp::cli::import::Args {
            genome_release: annonars::common::cli::GenomeRelease::Grch38,
            path_in_vcf: input_path.to_string_lossy().into(),
            path_out_rocksdb: output_path.to_string_lossy().into(),
            tbi_window_size: 100_000_000,
            cf_name: "dbsnp_data".to_string(),
            cf_name_by_rsid: "dbsnp_by_rsid".to_string(),
            path_wal_dir: None,
        };
        annonars::dbsnp::cli::import::run(&common_args, &import_args)?;

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

        // Test lookup by the RefSeq accession of chr1
        let var1_ac = keys::Var {
            chrom: "NC_000001.11".to_string(),
            ..var1
        };
        let res1_ac = annotator.annotate(&var1_ac)?.unwrap();
        assert_eq!(res1_ac.rs_id, "rs123");

        // Test multi-allelic separation lookup allele 1 (G)
        let var2 = keys::Var {
            chrom: "1".to_string(),
            pos: 10005,
            reference: "A".to_string(),
            alternative: "G".to_string(),
        };
        let res2 = annotator.annotate(&var2)?.unwrap();
        assert_eq!(res2.rs_id, "rs456");
        assert_eq!(res2.allele, "G");

        // Test multi-allelic separation lookup allele 2 (C)
        let var3 = keys::Var {
            chrom: "1".to_string(),
            pos: 10005,
            reference: "A".to_string(),
            alternative: "C".to_string(),
        };
        let res3 = annotator.annotate(&var3)?.unwrap();
        assert_eq!(res3.rs_id, "rs456");
        assert_eq!(res3.allele, "C");

        Ok(())
    }
}
