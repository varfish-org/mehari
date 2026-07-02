use crate::common::Args as CommonArgs;
use crate::db::keys::Var;
use crate::db::{ContigIdMap, PipelineConfig};
use crate::pbs::seqvars::CaddRecord;
use anyhow::Error;
use clap::Parser;
use prost::Message;
use std::collections::HashSet;
use std::fs::File;
use std::sync::Arc;

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
        no_progress: args.common.no_progress,
        threads: args.common.threads,
        db_type: "cadd",
        schema_version: "1.0",
        extra_meta: std::collections::HashMap::new(),
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

    let chrom_to_id: ContigIdMap = ContigIdMap::default();
    let chrom_to_id_closure = Arc::clone(&chrom_to_id);

    crate::db::run_tsv_pipeline(config, open_reader, move |record, _, contig_manager| {
        let row: CaddRecordRow = record.deserialize(None)?;
        let (chrom_std, chrom_id) =
            crate::db::get_or_intern_chrom(&row.chrom, contig_manager, &chrom_to_id_closure);

        let var = Var::new(
            chrom_std.clone(),
            row.pos,
            row.r#ref.clone(),
            row.alt.clone(),
        );
        let key = var.encode_with_id(chrom_id);

        let record_pb = CaddRecord {
            raw_score: row.raw_score,
            phred: row.phred,
        };

        let mut value = Vec::new();
        record_pb.encode(&mut value)?;

        let var_label = format!("{}:{}{}>{}", chrom_std, row.pos, row.r#ref, row.alt);
        Ok((vec![(key, value, var_label)], HashSet::new()))
    })?;

    tracing::info!("Writing contig index metadata mapping into the meta CF...");
    let options = rocksdb::Options::default();
    let db = rocksdb::DB::open_cf(&options, &args.common.output, vec!["meta", "cadd"])?;
    let cf_meta = db
        .cf_handle("meta")
        .ok_or_else(|| anyhow::anyhow!("meta CF not found"))?;

    let map_guard = chrom_to_id.read().unwrap();
    let serialized_dict = serde_json::to_vec(&*map_guard)?;
    db.put_cf(&cf_meta, b"contig_dictionary", serialized_dict)?;

    Ok(())
}
