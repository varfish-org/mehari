use clap::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader};
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

fn open_tsv_reader(
    path: &std::path::Path,
) -> Result<(csv::Reader<Box<dyn std::io::Read>>, Vec<String>), Error> {
    let file = File::open(path)?;
    let decoder: Box<dyn std::io::Read> = if path.extension().and_then(|s| s.to_str()) == Some("gz")
    {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut buf_reader = BufReader::new(decoder);

    let mut header_line = String::new();
    let mut line = String::new();
    while buf_reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim();
        if trimmed.starts_with("##") {
            line.clear();
            continue;
        }
        if trimmed.starts_with('#') {
            header_line = trimmed.to_string();
            break;
        }
        header_line = trimmed.to_string();
        break;
    }

    if header_line.is_empty() {
        anyhow::bail!("No header found in TSV file: {:?}", path);
    }

    let clean_header = header_line.trim_start_matches('#').to_string();
    let headers: Vec<String> = clean_header
        .split('\t')
        .map(|s| s.trim().to_string())
        .collect();

    let rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(Box::new(buf_reader) as Box<dyn std::io::Read>);

    Ok((rdr, headers))
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
        let (mut rdr, headers) = open_tsv_reader(input_file)?;

        let clean_headers: Vec<String> = headers
            .iter()
            .map(|h| h.trim_start_matches('#').to_string())
            .collect();

        let col_chrom = clean_headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("chrom"))
            .ok_or_else(|| anyhow!("Missing Chrom column in {:?}", input_file))?;
        let col_pos = clean_headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("pos"))
            .ok_or_else(|| anyhow!("Missing Pos column in {:?}", input_file))?;
        let col_ref = clean_headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("ref"))
            .ok_or_else(|| anyhow!("Missing Ref column in {:?}", input_file))?;
        let col_alt = clean_headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("alt"))
            .ok_or_else(|| anyhow!("Missing Alt column in {:?}", input_file))?;
        let col_raw = clean_headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("rawscore"))
            .ok_or_else(|| anyhow!("Missing RawScore column in {:?}", input_file))?;
        let col_phred = clean_headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("phred"))
            .ok_or_else(|| anyhow!("Missing PHRED column in {:?}", input_file))?;

        for result in rdr.records() {
            let record = result?;
            if record.is_empty() || record[0].starts_with('#') {
                continue;
            }

            let chrom = &record[col_chrom];
            let pos: i32 = record[col_pos].parse()?;
            let reference = &record[col_ref];
            let alternative = &record[col_alt];
            let raw_score: f32 = record[col_raw].parse()?;
            let phred: f32 = record[col_phred].parse()?;

            // Normalize chrom to primary name
            let chrom_std = contig_manager
                .get_primary_name(chrom)
                .cloned()
                .unwrap_or_else(|| chrom.to_string());

            let var = Var {
                chrom: chrom_std,
                pos,
                reference: reference.to_string(),
                alternative: alternative.to_string(),
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

            let record_pb = CaddRecord { raw_score, phred };
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
