use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;
use std::time::Instant;

use crate::common::Args as CommonArgs;
use crate::common::contig::ContigManager;
use crate::db::{DbWriter, finalize_db, get_total_records_from_tabix, open_db, open_vcf_reader};
use crate::pbs::seqvars::{SpliceAiPrediction, SpliceAiRecord};
use annonars::common::keys::Var;
use anyhow::{Error, anyhow};
use prost::Message;

/// Command line arguments for `db spliceai create` subcommand.
#[derive(Parser, Debug, Clone)]
#[command(about = "Construct SpliceAI score RocksDB database", long_about = None)]
pub struct Args {
    /// Assembly to use (e.g. GRCh37 or GRCh38)
    #[arg(long, required = true)]
    pub assembly: String,

    /// Path to SpliceAI VCF/VCF.gz input file(s)
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
    tracing::info!(
        "Creating SpliceAI score RocksDB database at {:?}",
        args.output
    );
    let start_time = Instant::now();

    // Initialize ContigManager for chromosome name normalization
    let contig_manager = ContigManager::new(&args.assembly);

    let db = open_db(&args.output, "spliceai")?;
    let mut writer = DbWriter::new(&db, "spliceai", args.batch_size)?;

    for input_file in &args.input {
        tracing::info!("Processing input file: {:?}", input_file);

        // Determine total records using our Tabix helper
        let total_records = get_total_records_from_tabix(input_file).unwrap_or(0);

        let pb = if total_records > 0 {
            ProgressBar::new(total_records)
        } else {
            ProgressBar::new_spinner()
        };

        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")?
                .progress_chars("█▒░")
        );

        let (mut reader, header) = open_vcf_reader(input_file)?;

        let mut chunk = Vec::with_capacity(args.batch_size);

        for result in reader.record_bufs(&header) {
            chunk.push(result?);

            if chunk.len() == args.batch_size {
                let chunk_len = chunk.len() as u64;
                write_chunk(&mut writer, &chunk, &contig_manager)?;
                pb.inc(chunk_len);
                chunk.clear();
            }
        }

        if !chunk.is_empty() {
            let chunk_len = chunk.len() as u64;
            write_chunk(&mut writer, &chunk, &contig_manager)?;
            pb.inc(chunk_len);
            chunk.clear();
        }

        pb.finish_and_clear();
    }

    let written = writer.flush()?;

    let cf_meta = db
        .cf_handle("meta")
        .ok_or_else(|| anyhow!("meta CF not found"))?;
    db.put_cf(&cf_meta, b"db_type", b"spliceai")?;
    db.put_cf(&cf_meta, b"assembly", args.assembly.as_bytes())?;
    db.put_cf(&cf_meta, b"schema_version", b"1.0")?;

    tracing::info!(
        "Successfully completed SpliceAI import of {} records in {:?}",
        written,
        start_time.elapsed()
    );
    finalize_db(&db, &["spliceai", "meta"])?;

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

            let spliceai_val = match record.info().get("SpliceAI").flatten() {
                Some(val) => val,
                None => return Ok(kvs),
            };

            let spliceai_str = match spliceai_val {
                noodles::vcf::variant::record_buf::info::field::Value::String(s) => s.to_string(),
                noodles::vcf::variant::record_buf::info::field::Value::Array(
                    noodles::vcf::variant::record_buf::info::field::value::Array::String(arr),
                ) => arr.iter().flatten().cloned().collect::<Vec<_>>().join(","),
                _ => return Ok(kvs),
            };

            // Normalize chrom to primary name
            let chrom_std = contig_manager
                .get_primary_name(chrom)
                .cloned()
                .unwrap_or_else(|| chrom.to_string());

            // Parse individual predictions (multiple comma-separated entries are possible)
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
