use crate::common::Args as CommonArgs;
use crate::db::PipelineConfig;
use crate::db::keys::Var;
use crate::pbs::seqvars::CaddRecord;
use anyhow::Error;
use clap::Parser;
use prost::Message;
use std::collections::{HashMap, HashSet};
use std::fs::File;

/// Arguments for the CADD database construction command.
#[derive(Parser, Debug, Clone)]
#[command(about = "Construct CADD score RocksDB database", long_about = None)]
pub struct Args {
    #[command(flatten)]
    pub common: crate::db::CommonPipelineArgs,
}

pub mod cli {
    pub use super::Args;
}

#[derive(serde::Deserialize, Clone)]
struct CaddRecordRow {
    chrom: String,
    pos: i32,
    r#ref: String,
    alt: String,
    raw_score: f32,
    phred: f32,
}

pub fn run(_common: &CommonArgs, args: &Args) -> Result<(), Error> {
    let config = PipelineConfig {
        assembly: &args.common.assembly,
        input: &args.common.input,
        output: &args.common.output,
        batch_size: args.common.batch_size,
        quiet: args.common.quiet,
        threads: args.common.threads,
        db_type: "cadd",
        schema_version: "1.0",
        extra_meta: HashMap::new(),
    };

    let open_reader = |path: &std::path::Path| {
        let file = File::open(path)?;
        let (reader, _) = niffler::get_reader(Box::new(file))?;
        let rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_reader(reader);
        Ok((rdr, csv::StringRecord::new()))
    };

    crate::db::run_tsv_pipeline(config, open_reader, |record, _, contig_manager| {
        let row: CaddRecordRow = record.deserialize(None)?;
        let chrom_std = contig_manager
            .get_primary_name(&row.chrom)
            .cloned()
            .unwrap_or_else(|| row.chrom.clone());

        let var = Var {
            chrom: chrom_std,
            pos: row.pos,
            reference: row.r#ref.clone(),
            alternative: row.alt.clone(),
        };
        let key: Vec<u8> = var.clone().into();
        let record_pb = CaddRecord {
            raw_score: row.raw_score,
            phred: row.phred,
        };

        let mut value = Vec::new();
        record_pb.encode(&mut value)?;

        let var_label = format!(
            "{}:{}{}>{}",
            var.chrom, var.pos, var.reference, var.alternative
        );
        Ok((vec![(key, value, var_label)], HashSet::new()))
    })
}
