use crate::common::Args as CommonArgs;
use crate::db::PipelineConfig;
use crate::db::keys::Var;
use crate::pbs::seqvars::{SpliceAiPrediction, SpliceAiRecord};
use anyhow::Error;
use clap::Parser;
use prost::Message;
use std::collections::{HashMap, HashSet};

/// Arguments for the SpliceAI database construction command.
#[derive(Parser, Debug, Clone)]
#[command(about = "Construct SpliceAI score RocksDB database", long_about = None)]
pub struct Args {
    #[command(flatten)]
    pub common: crate::db::CommonPipelineArgs,
}

pub mod cli {
    pub use super::Args;
}

pub fn run(_common: &CommonArgs, args: &Args) -> Result<(), Error> {
    let config = PipelineConfig {
        assembly: &args.common.assembly,
        input: &args.common.input,
        output: &args.common.output,
        batch_size: args.common.batch_size,
        quiet: args.common.quiet,
        threads: args.common.threads,
        db_type: "spliceai",
        schema_version: "1.0",
        extra_meta: HashMap::new(),
    };

    crate::db::run_vcf_pipeline(
        config,
        None::<fn(&mut noodles::vcf::Header)>,
        |record, contig_manager| {
            let mut kvs = Vec::new();
            let chrom = record.reference_sequence_name();
            let pos = match record.variant_start() {
                Some(start) => start.get() as i32,
                None => return Ok((kvs, HashSet::new())),
            };

            let spliceai_val = match record.info().get("SpliceAI").flatten() {
                Some(val) => val,
                None => return Ok((kvs, HashSet::new())),
            };

            let spliceai_str = match crate::db::get_info_string(spliceai_val) {
                Some(s) => s,
                None => return Ok((kvs, HashSet::new())),
            };

            let chrom_std = contig_manager
                .get_primary_name(chrom)
                .cloned()
                .unwrap_or_else(|| chrom.to_string());
            let reference = record.reference_bases();
            let mut predictions_by_allele: HashMap<String, Vec<SpliceAiPrediction>> =
                HashMap::new();

            for pred_str in spliceai_str.split(',') {
                let fields: Vec<&str> = pred_str.split('|').collect();
                if fields.len() < 10 {
                    continue;
                }

                let allele = fields[0].to_string();
                let prediction = SpliceAiPrediction {
                    allele: allele.clone(),
                    symbol: fields[1].to_string(),
                    ds_ag: fields[2].parse()?,
                    ds_al: fields[3].parse()?,
                    ds_dg: fields[4].parse()?,
                    ds_dl: fields[5].parse()?,
                    dp_ag: fields[6].parse()?,
                    dp_al: fields[7].parse()?,
                    dp_dg: fields[8].parse()?,
                    dp_dl: fields[9].parse()?,
                };
                predictions_by_allele
                    .entry(allele)
                    .or_default()
                    .push(prediction);
            }

            for (allele, predictions) in predictions_by_allele {
                let var = Var {
                    chrom: chrom_std.clone(),
                    pos,
                    reference: reference.to_string(),
                    alternative: allele,
                };
                let key: Vec<u8> = var.clone().into();
                let record_pb = SpliceAiRecord { predictions };

                let mut value = Vec::new();
                record_pb.encode(&mut value)?;

                let var_label = format!(
                    "{}:{}{}>{}",
                    var.chrom, var.pos, var.reference, var.alternative
                );
                kvs.push((key, value, var_label));
            }
            Ok((kvs, HashSet::new()))
        },
    )
}
