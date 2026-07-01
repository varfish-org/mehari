use clap::Parser;
use rayon::prelude::*;
use std::path::PathBuf;
use std::time::Instant;

use crate::common::Args as CommonArgs;
use crate::common::contig::ContigManager;
use crate::db::{DbWriter, finalize_db, open_db, open_vcf_reader};
use crate::pbs::seqvars::DbsnpRecord;
use annonars::common::keys::Var;
use anyhow::{Error, anyhow};
use itertools::Itertools;
use noodles::vcf::variant::record::Ids;
use prost::Message;

/// Command line arguments for `db dbsnp create` subcommand.
#[derive(Parser, Debug, Clone)]
#[command(about = "Construct dbSNP RocksDB database", long_about = None)]
pub struct Args {
    /// Assembly to use (e.g. GRCh37 or GRCh38)
    #[arg(long, required = true)]
    pub assembly: String,

    /// Path to dbSNP VCF/VCF.gz input file(s)
    #[arg(long, required = true)]
    pub input: Vec<PathBuf>,

    /// Path to output RocksDB directory
    #[arg(long, required = true)]
    pub output: PathBuf,

    /// Number of rows to write in a single RocksDB batch
    #[arg(long, default_value = "100000")]
    pub batch_size: usize,
}

pub mod cli {
    pub use super::Args;
}

pub fn run(_common: &CommonArgs, args: &Args) -> Result<(), Error> {
    tracing::info!("Creating dbSNP RocksDB database at {:?}", args.output);
    let start_time = Instant::now();
    let contig_manager = ContigManager::new(&args.assembly);

    let db = open_db(&args.output, "dbsnp")?;
    let mut writer = DbWriter::new(&db, "dbsnp", args.batch_size)?;

    for input_file in &args.input {
        tracing::info!("Processing input file: {:?}", input_file);
        let (mut reader, header) = open_vcf_reader(input_file)?;

        // dbSnp's rs IDs are larger than the VCF spec max of 32bit signed integers,
        // so we modify the type to String here.
        let mut modified_header = header.clone();
        if let Some(info) = modified_header.infos_mut().get_mut("RS") {
            *info.type_mut() = noodles::vcf::header::record::value::map::info::Type::String;
        }

        let mut chunk = Vec::with_capacity(args.batch_size);

        for result in reader.record_bufs(&modified_header) {
            chunk.push(result?);

            if chunk.len() == args.batch_size {
                write_chunk(&mut writer, &chunk, &contig_manager)?;
                chunk.clear();
            }
        }

        if !chunk.is_empty() {
            write_chunk(&mut writer, &chunk, &contig_manager)?;
            chunk.clear();
        }
    }

    let written = writer.flush()?;

    // Commit general metadata block
    let cf_meta = db
        .cf_handle("meta")
        .ok_or_else(|| anyhow!("meta CF not found"))?;
    db.put_cf(&cf_meta, b"db_type", b"dbsnp")?;
    db.put_cf(&cf_meta, b"assembly", args.assembly.as_bytes())?;
    db.put_cf(&cf_meta, b"schema_version", b"3.0")?;

    tracing::info!(
        "Successfully completed lightweight dbSNP import of {} records in {:?}",
        written,
        start_time.elapsed()
    );
    finalize_db(&db, &["dbsnp", "meta"])?;

    Ok(())
}

fn write_chunk(
    writer: &mut DbWriter,
    chunk: &[noodles::vcf::variant::RecordBuf],
    contig_manager: &ContigManager,
) -> Result<(), Error> {
    let processed_records: Vec<Result<Vec<(Vec<u8>, Vec<u8>, String)>, Error>> = chunk
        .par_iter()
        .map(|record| {
            let mut kvs = Vec::new();
            let chrom = record.reference_sequence_name();
            let pos = match record.variant_start() {
                Some(start) => start.get() as i32,
                None => return Ok(kvs),
            };
            let reference = record.reference_bases();
            let alternative = record.alternate_bases();

            let rs_id = if record.ids().is_empty() {
                "".to_string()
            } else {
                record.ids().iter().join(";")
            };

            if rs_id.is_empty() {
                return Ok(kvs);
            }

            let chrom_std = contig_manager
                .get_primary_name(chrom)
                .cloned()
                .unwrap_or_else(|| chrom.to_string());

            // Explode multi-allelic lines into individual variant keys
            for alt in alternative.as_ref() {
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
            Ok(kvs)
        })
        .collect();

    for batch_result in processed_records {
        for (key, value, var_label) in batch_result? {
            writer.put(&key, &value, &var_label)?;
        }
    }

    Ok(())
}
