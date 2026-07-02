use crate::common::Args as CommonArgs;
use crate::db::PipelineConfig;
use crate::db::keys::Var;
use crate::pbs::seqvars::DbsnpRecord;
use anyhow::Error;
use clap::Parser;
use itertools::Itertools;
use noodles::vcf::variant::record::Ids;
use prost::Message;
use rustc_hash::FxHashMap;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, RwLock};

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

    let chrom_to_id = Arc::new(RwLock::new(FxHashMap::default()));
    let chrom_to_id_closure = Arc::clone(&chrom_to_id);

    crate::db::run_vcf_pipeline(
        config,
        Some(header_modifier),
        move |record, contig_manager| {
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

            let (chrom_std, chrom_id) =
                crate::db::get_or_intern_chrom(chrom, contig_manager, &chrom_to_id_closure);

            let reference = record.reference_bases();

            for alt in record.alternate_bases().as_ref() {
                let alt_str = alt.to_string();
                let var = Var::new(
                    chrom_std.clone(),
                    pos,
                    reference.to_string(),
                    alt_str.clone(),
                );
                let key = var.encode_with_id(chrom_id);
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
        },
    )?;

    tracing::info!("Writing dbSNP contig index metadata mapping into the meta CF...");
    let options = rocksdb::Options::default();
    let db = rocksdb::DB::open_cf(&options, &args.common.output, vec!["meta", "dbsnp"])?;
    let cf_meta = db
        .cf_handle("meta")
        .ok_or_else(|| anyhow::anyhow!("meta CF not found"))?;

    let map_guard = chrom_to_id.read().unwrap();
    let serialized_dict = serde_json::to_vec(&*map_guard)?;
    db.put_cf(&cf_meta, b"contig_dictionary", serialized_dict)?;

    Ok(())
}
