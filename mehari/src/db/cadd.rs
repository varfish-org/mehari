use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::fs::File;
use std::path::PathBuf;
use std::time::Instant;

use crate::common::Args as CommonArgs;
use crate::common::contig::ContigManager;
use crate::db::{DbWriter, finalize_db, get_total_records_from_tabix, open_db};
use crate::pbs::seqvars::CaddRecord;
use annonars::common::keys::Var;
use anyhow::{Error, anyhow};
use prost::Message;

/// Command line arguments for `db cadd create` subcommand.
#[derive(Parser, Debug, Clone)]
#[command(about = "Construct CADD score RocksDB database", long_about = None)]
pub struct Args {
    /// Assembly to use (e.g. GRCh37 or GRCh38)
    #[arg(long, required = true)]
    pub assembly: String,

    /// Path to CADD TSV/TSV.gz input file(s)
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
    tracing::info!("Creating CADD score RocksDB database at {:?}", args.output);
    let start_time = Instant::now();

    // Initialize ContigManager for chromosome name normalization
    let contig_manager = ContigManager::new(&args.assembly);

    let db = open_db(&args.output, "cadd")?;
    let mut writer = DbWriter::new(&db, "cadd", args.batch_size)?;

    for input_file in &args.input {
        tracing::info!("Processing input file: {:?}", input_file);

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

        let file = File::open(input_file)?;
        let (reader, _format) = niffler::get_reader(Box::new(file))?;

        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_reader(reader);

        let mut chunk = Vec::with_capacity(args.batch_size);

        for result in rdr.deserialize::<CaddRecordRow>() {
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
    db.put_cf(&cf_meta, b"db_type", b"cadd")?;
    db.put_cf(&cf_meta, b"assembly", args.assembly.as_bytes())?;
    db.put_cf(&cf_meta, b"schema_version", b"1.0")?;

    tracing::info!(
        "Successfully completed CADD import of {} records in {:?}",
        written,
        start_time.elapsed()
    );
    finalize_db(&db, &["cadd", "meta"])?;

    Ok(())
}

fn write_chunk(
    writer: &mut DbWriter,
    chunk: &[CaddRecordRow],
    contig_manager: &ContigManager,
) -> Result<(), Error> {
    let processed_records: Vec<Result<(Vec<u8>, Vec<u8>, String), Error>> = chunk
        .par_iter()
        .map(|row| {
            // Normalize chrom to primary name
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
            Ok((key, value, var_label))
        })
        .collect();

    for item in processed_records {
        let (key, value, var_label) = item?;
        writer.put(&key, &value, &var_label)?;
    }

    Ok(())
}
