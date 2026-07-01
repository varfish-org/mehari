use crate::common::Args as CommonArgs;
use crate::db::PipelineConfig;
use crate::db::keys::Var;
use crate::pbs::seqvars::DbsnpRecord;
use anyhow::Error;
use clap::Parser;
use itertools::Itertools;
use noodles::vcf::variant::record::Ids;
use prost::Message;
use std::collections::{HashMap, HashSet};

/// Arguments for the dbSNP database construction command.
#[derive(Parser, Debug, Clone)]
#[command(about = "Construct dbSNP RocksDB database", long_about = None)]
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
        no_progress: args.common.no_progress,
        threads: args.common.threads,
        db_type: "dbsnp",
        schema_version: "3.0",
        extra_meta: HashMap::new(),
    };

    let header_modifier = |header: &mut noodles::vcf::Header| {
        if let Some(info) = header.infos_mut().get_mut("RS") {
            *info.type_mut() = noodles::vcf::header::record::value::map::info::Type::String;
        }
    };

    crate::db::run_vcf_pipeline(config, Some(header_modifier), |record, contig_manager| {
        let mut kvs = Vec::new();
        let chrom = record.reference_sequence_name();
        let pos = match record.variant_start() {
            Some(start) => start.get() as i32,
            None => return Ok((kvs, HashSet::new())),
        };

        let rs_id = if record.ids().is_empty() {
            "".to_string()
        } else {
            record.ids().iter().join(";")
        };
        if rs_id.is_empty() {
            return Ok((kvs, HashSet::new()));
        }

        let chrom_std = contig_manager
            .get_primary_name(chrom)
            .cloned()
            .unwrap_or_else(|| chrom.to_string());
        let reference = record.reference_bases();

        for alt in record.alternate_bases().as_ref() {
            let alt_str = alt.to_string();
            let var = Var {
                chrom: chrom_std.clone(),
                pos,
                reference: reference.to_string(),
                alternative: alt_str.clone(),
            };
            let key: Vec<u8> = var.clone().into();
            let record_pb = DbsnpRecord {
                allele: alt_str,
                rs_id: rs_id.clone(),
            };

            let mut value = Vec::new();
            record_pb.encode(&mut value)?;

            let var_label = format!(
                "{}:{}{}>{}",
                var.chrom, var.pos, var.reference, var.alternative
            );
            kvs.push((key, value, var_label));
        }
        Ok((kvs, HashSet::new()))
    })
}
