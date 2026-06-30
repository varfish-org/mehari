use std::path::Path;
use std::sync::Arc;

use crate::common::contig::ContigManager;
use annonars::common::keys;
use anyhow::Error;
use prost::Message;
use rocksdb::{DBWithThreadMode, MultiThreaded};

#[derive(Debug)]
pub struct SpliceAiAnnotator {
    db: DBWithThreadMode<MultiThreaded>,
    contig_manager: Arc<ContigManager>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SpliceAiPrediction {
    pub allele: String,
    pub symbol: String,
    pub ds_ag: f32,
    pub ds_al: f32,
    pub ds_dg: f32,
    pub ds_dl: f32,
    pub dp_ag: i32,
    pub dp_al: i32,
    pub dp_dg: i32,
    pub dp_dl: i32,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SpliceAiResult {
    pub predictions: Vec<SpliceAiPrediction>,
}

impl SpliceAiAnnotator {
    pub fn new(db: DBWithThreadMode<MultiThreaded>, contig_manager: Arc<ContigManager>) -> Self {
        Self { db, contig_manager }
    }

    pub(crate) fn from_path(
        path: impl AsRef<Path>,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        tracing::info!("Opening SpliceAI database at {}", path.as_ref().display());
        let options = rocksdb::Options::default();
        let db_spliceai =
            rocksdb::DB::open_cf_for_read_only(&options, &path, ["meta", "spliceai"], false)?;
        Ok(Self::new(db_spliceai, contig_manager))
    }

    pub fn annotate_record_spliceai(&self, key: &[u8]) -> Result<Option<SpliceAiResult>, Error> {
        if let Some(raw_value) = self
            .db
            .get_cf(self.db.cf_handle("spliceai").as_ref().unwrap(), key)?
        {
            let record = crate::pbs::seqvars::SpliceAiRecord::decode(&raw_value[..])?;
            let predictions = record
                .predictions
                .into_iter()
                .map(|p| SpliceAiPrediction {
                    allele: p.allele,
                    symbol: p.symbol,
                    ds_ag: p.ds_ag,
                    ds_al: p.ds_al,
                    ds_dg: p.ds_dg,
                    ds_dl: p.ds_dl,
                    dp_ag: p.dp_ag,
                    dp_al: p.dp_al,
                    dp_dg: p.dp_dg,
                    dp_dl: p.dp_dl,
                })
                .collect();
            Ok(Some(SpliceAiResult { predictions }))
        } else {
            Ok(None)
        }
    }

    pub fn annotate(&self, vcf_var: &keys::Var) -> anyhow::Result<Option<SpliceAiResult>> {
        if self
            .contig_manager
            .is_canonical_alias(vcf_var.chrom.as_str())
        {
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
            return self.annotate_record_spliceai(&key);
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
    fn test_spliceai_import_and_annotate() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let input_path = temp.join("spliceai_test.vcf");
        let output_path = temp.join("spliceai_db");

        let mut file = File::create(&input_path)?;
        writeln!(file, "##fileformat=VCFv4.2")?;
        writeln!(
            file,
            "##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAI prediction\">"
        )?;
        writeln!(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
        writeln!(
            file,
            "1\t10000\t.\tC\tT\t.\t.\tSpliceAI=T|GENEA|0.10|0.20|0.30|0.40|10|20|30|40"
        )?;
        writeln!(
            file,
            "chr1\t10005\t.\tA\tG,C\t.\t.\tSpliceAI=G|GENEB|0.05|0.00|0.15|0.00|5|0|15|0,C|GENEC|0.50|0.00|0.00|0.00|2|0|0|0"
        )?;
        file.flush()?;

        let common_args = crate::common::Args {
            verbose: clap_verbosity_flag::Verbosity::new(0, 0),
        };
        let create_args = crate::db::create_spliceai::Args {
            assembly: "GRCh38".to_string(),
            input: vec![input_path],
            output: output_path.clone(),
            batch_size: 1000,
        };

        crate::db::create_spliceai::run(&common_args, &create_args)?;

        let contig_manager = Arc::new(ContigManager::new("GRCh38"));
        let annotator = SpliceAiAnnotator::from_path(&output_path, contig_manager)?;

        let var1 = keys::Var {
            chrom: "chr1".to_string(),
            pos: 10000,
            reference: "C".to_string(),
            alternative: "T".to_string(),
        };
        let res1 = annotator.annotate(&var1)?.unwrap();
        assert_eq!(res1.predictions.len(), 1);
        assert_eq!(res1.predictions[0].symbol, "GENEA");
        assert_eq!(res1.predictions[0].ds_ag, 0.10);
        assert_eq!(res1.predictions[0].dp_ag, 10);

        let var2 = keys::Var {
            chrom: "1".to_string(),
            pos: 10005,
            reference: "A".to_string(),
            alternative: "G".to_string(),
        };
        let res2 = annotator.annotate(&var2)?.unwrap();
        assert_eq!(res2.predictions.len(), 1);
        assert_eq!(res2.predictions[0].symbol, "GENEB");
        assert_eq!(res2.predictions[0].ds_dg, 0.15);

        let var3 = keys::Var {
            chrom: "1".to_string(),
            pos: 10005,
            reference: "A".to_string(),
            alternative: "C".to_string(),
        };
        let res3 = annotator.annotate(&var3)?.unwrap();
        assert_eq!(res3.predictions.len(), 1);
        assert_eq!(res3.predictions[0].symbol, "GENEC");
        assert_eq!(res3.predictions[0].ds_ag, 0.50);

        Ok(())
    }
}
