use std::path::Path;
use std::sync::Arc;

use crate::common::contig::ContigManager;
use crate::db::keys;
use annonars::common::contig::ContigDict;
use anyhow::Error;
use prost::Message;
use rocksdb::{DBWithThreadMode, MultiThreaded};

#[derive(Debug)]
pub struct SpliceAiAnnotator {
    db: DBWithThreadMode<MultiThreaded>,
    contig_manager: Arc<ContigManager>,
    contig_dict: ContigDict,
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
    pub fn new(
        db: DBWithThreadMode<MultiThreaded>,
        contig_manager: Arc<ContigManager>,
        contig_dict: ContigDict,
    ) -> Self {
        Self {
            db,
            contig_manager,
            contig_dict,
        }
    }

    /// Open a SpliceAI track database built with `annonars spliceai import`.
    pub(crate) fn from_path(
        path: impl AsRef<Path>,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        tracing::info!("Opening SpliceAI database at {}", path.as_ref().display());
        let options = rocksdb::Options::default();
        let db_spliceai =
            rocksdb::DB::open_cf_for_read_only(&options, &path, ["meta", "spliceai"], false)?;
        let contig_dict = annonars::seqvars::read_contig_dict(&db_spliceai)?;

        Ok(Self::new(db_spliceai, contig_manager, contig_dict))
    }

    /// Look up `key`; the stored predictions do not carry their allele, so `allele` is filled in.
    pub fn annotate_record_spliceai(
        &self,
        key: &[u8],
        allele: &str,
    ) -> Result<Option<SpliceAiResult>, Error> {
        if let Some(raw_value) = self
            .db
            .get_cf(self.db.cf_handle("spliceai").as_ref().unwrap(), key)?
        {
            let record = annonars::pbs::seqvars::base::SpliceAiRecord::decode(&raw_value[..])?;
            let predictions = record
                .predictions
                .into_iter()
                .map(|p| SpliceAiPrediction {
                    allele: allele.to_string(),
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
        // Normalize chrom to primary name to build standard key
        let chrom_std = self
            .contig_manager
            .get_primary_name(&vcf_var.chrom)
            .cloned()
            .unwrap_or_else(|| vcf_var.chrom.clone());

        if let Some(chrom_id) = self.contig_dict.id_of(&chrom_std) {
            let annonars_var: annonars::common::keys::Var = vcf_var.into();
            let key = annonars_var.encode_with_id(chrom_id);
            self.annotate_record_spliceai(&key, &vcf_var.alternative)
        } else {
            Ok(None)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use temp_testdir::TempDir;

    #[test]
    fn test_spliceai_import_and_annotate() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let fai_path = temp.join("ref.fa.fai");
        let input_path = temp.join("spliceai_test.vcf");
        let output_path = temp.join("spliceai_db");

        std::fs::write(&fai_path, "chr1\t248956422\t6\t60\t61\n")?;
        let vcf_content = "##fileformat=VCFv4.2\n\
                                ##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAI prediction\">\n\
                                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                                1\t10000\t.\tC\tT\t.\t.\tSpliceAI=T|GENEA|0.10|0.20|0.30|0.40|10|20|30|40\n\
                                chr1\t10005\t.\tA\tG,C\t.\t.\tSpliceAI=G|GENEB|0.05|0.00|0.15|0.00|5|0|15|0,C|GENEC|0.50|0.00|0.00|0.00|2|0|0|0\n";
        std::fs::write(&input_path, vcf_content)?;

        let common_args = annonars::common::cli::Args {
            verbose: Default::default(),
        };
        let import_args = annonars::spliceai::cli::import::Args {
            path_in_vcf: input_path.to_string_lossy().into(),
            path_reference_fai: fai_path.to_string_lossy().into(),
            assembly: "GRCh38".to_string(),
            path_out_rocksdb: output_path.to_string_lossy().into(),
            cf_name: "spliceai".to_string(),
            path_wal_dir: None,
        };
        annonars::spliceai::cli::import::run(&common_args, &import_args)?;

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
        assert_eq!(res3.predictions[0].allele, "C");
        assert_eq!(res3.predictions[0].symbol, "GENEC");
        assert_eq!(res3.predictions[0].ds_ag, 0.50);

        Ok(())
    }
}
