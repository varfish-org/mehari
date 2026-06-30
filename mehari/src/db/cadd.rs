use clap::Parser;
use std::fs::File;
use std::path::PathBuf;
use std::time::Instant;

use crate::common::Args as CommonArgs;
use crate::common::contig::ContigManager;
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

#[derive(serde::Deserialize)]
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

    // Open RocksDB database
    let mut options = rocksdb::Options::default();
    options.create_if_missing(true);
    options.create_missing_column_families(true);
    options.set_compression_type(rocksdb::DBCompressionType::Zstd);
    options.set_compression_options(-14, 19, 0, 0);

    let cfs = vec!["meta", "cadd"];
    let db = rocksdb::DB::open_cf(&options, &args.output, cfs)?;

    let cf_cadd = db
        .cf_handle("cadd")
        .ok_or_else(|| anyhow!("cadd CF not found"))?;
    let cf_meta = db
        .cf_handle("meta")
        .ok_or_else(|| anyhow!("meta CF not found"))?;

    let mut batch = rocksdb::WriteBatch::default();
    let mut count = 0;
    let mut written = 0;

    for input_file in &args.input {
        tracing::info!("Processing input file: {:?}", input_file);
        let file = File::open(input_file)?;
        let (reader, _format) = niffler::get_reader(Box::new(file))?;

        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_reader(reader);

        for result in rdr.deserialize::<CaddRecordRow>() {
            let row = result?;

            // Normalize chrom to primary name
            let chrom_std = contig_manager
                .get_primary_name(&row.chrom)
                .cloned()
                .unwrap_or_else(|| row.chrom.clone());

            let var = Var {
                chrom: chrom_std,
                pos: row.pos,
                reference: row.r#ref,
                alternative: row.alt,
            };
            let key: Vec<u8> = var.clone().into();

            // Check for clashing/duplicate keys
            if db.get_cf(&cf_cadd, &key)?.is_some() {
                tracing::warn!(
                    "Duplicate key found in database for variant: {}:{}{}>{}",
                    var.chrom,
                    var.pos,
                    var.reference,
                    var.alternative
                );
            }

            let record_pb = CaddRecord {
                raw_score: row.raw_score,
                phred: row.phred,
            };
            let mut value = Vec::new();
            record_pb.encode(&mut value)?;

            batch.put_cf(&cf_cadd, &key, &value);
            count += 1;

            if count % args.batch_size == 0 {
                db.write(batch)?;
                batch = rocksdb::WriteBatch::default();
                written += count;
                tracing::info!("Imported {} records...", written);
                count = 0;
            }
        }
    }

    if count > 0 {
        db.write(batch)?;
        written += count;
    }

    // Write metadata
    db.put_cf(&cf_meta, b"db_type", b"cadd")?;
    db.put_cf(&cf_meta, b"assembly", args.assembly.as_bytes())?;
    db.put_cf(&cf_meta, b"schema_version", b"1.0")?;

    tracing::info!(
        "Successfully completed CADD import of {} records in {:?}",
        written,
        start_time.elapsed()
    );

    Ok(())
}
