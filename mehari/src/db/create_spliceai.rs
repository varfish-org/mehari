use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::time::Instant;

use crate::common::Args as CommonArgs;
use crate::common::contig::ContigManager;
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

fn open_vcf_reader(
    path: &std::path::Path,
) -> Result<
    (
        noodles::vcf::io::Reader<Box<dyn std::io::BufRead>>,
        noodles::vcf::Header,
    ),
    Error,
> {
    let file = File::open(path)?;
    let decoder: Box<dyn std::io::Read> = if path.extension().and_then(|s| s.to_str()) == Some("gz")
    {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buf_reader = BufReader::new(decoder);
    let mut reader =
        noodles::vcf::io::Reader::new(Box::new(buf_reader) as Box<dyn std::io::BufRead>);
    let header = reader.read_header()?;
    Ok((reader, header))
}

pub fn run(_common: &CommonArgs, args: &Args) -> Result<(), Error> {
    tracing::info!(
        "Creating SpliceAI score RocksDB database at {:?}",
        args.output
    );
    let start_time = Instant::now();

    // Initialize ContigManager for chromosome name normalization
    let contig_manager = ContigManager::new(&args.assembly);

    // Open RocksDB database
    let mut options = rocksdb::Options::default();
    options.create_if_missing(true);
    options.create_missing_column_families(true);
    options.set_compression_type(rocksdb::DBCompressionType::Zstd);
    options.set_compression_options(-14, 19, 0, 0);

    let cfs = vec!["meta", "spliceai"];
    let db = rocksdb::DB::open_cf(&options, &args.output, cfs)?;

    let cf_spliceai = db
        .cf_handle("spliceai")
        .ok_or_else(|| anyhow!("spliceai CF not found"))?;
    let cf_meta = db
        .cf_handle("meta")
        .ok_or_else(|| anyhow!("meta CF not found"))?;

    let mut batch = rocksdb::WriteBatch::default();
    let mut count = 0;
    let mut written = 0;

    for input_file in &args.input {
        tracing::info!("Processing input file: {:?}", input_file);
        let (mut reader, header) = open_vcf_reader(input_file)?;

        for result in reader.record_bufs(&header) {
            let record = result?;
            let chrom = record.reference_sequence_name();
            let pos = match record.variant_start() {
                Some(start) => start.get() as i32,
                None => continue,
            };
            let reference = record.reference_bases();

            let spliceai_val = match record.info().get("SpliceAI").flatten() {
                Some(val) => val,
                None => continue,
            };

            let spliceai_str = match spliceai_val {
                noodles::vcf::variant::record_buf::info::field::Value::String(s) => s.to_string(),
                noodles::vcf::variant::record_buf::info::field::Value::Array(
                    noodles::vcf::variant::record_buf::info::field::value::Array::String(arr),
                ) => arr.iter().flatten().cloned().collect::<Vec<_>>().join(","),
                _ => continue,
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
                let symbol = fields[1].to_string();
                let ds_ag: f32 = fields[2].parse()?;
                let ds_al: f32 = fields[3].parse()?;
                let ds_dg: f32 = fields[4].parse()?;
                let ds_dl: f32 = fields[5].parse()?;
                let dp_ag: i32 = fields[6].parse()?;
                let dp_al: i32 = fields[7].parse()?;
                let dp_dg: i32 = fields[8].parse()?;
                let dp_dl: i32 = fields[9].parse()?;

                let prediction = SpliceAiPrediction {
                    allele: allele.clone(),
                    symbol,
                    ds_ag,
                    ds_al,
                    ds_dg,
                    ds_dl,
                    dp_ag,
                    dp_al,
                    dp_dg,
                    dp_dl,
                };

                predictions_by_allele
                    .entry(allele)
                    .or_default()
                    .push(prediction);
            }

            // For each allele, write its predictions to RocksDB
            for (allele, predictions) in predictions_by_allele {
                let var = Var {
                    chrom: chrom_std.clone(),
                    pos,
                    reference: reference.to_string(),
                    alternative: allele,
                };
                let key: Vec<u8> = var.clone().into();

                // Check for clashing/duplicate keys
                if db.get_cf(&cf_spliceai, &key)?.is_some() {
                    tracing::warn!(
                        "Duplicate key found in database for variant: {}:{}{}>{}",
                        var.chrom,
                        var.pos,
                        var.reference,
                        var.alternative
                    );
                }

                let record_pb = SpliceAiRecord { predictions };
                let mut value = Vec::new();
                record_pb.encode(&mut value)?;

                batch.put_cf(&cf_spliceai, &key, &value);
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
    }

    if count > 0 {
        db.write(batch)?;
        written += count;
    }

    // Write metadata
    db.put_cf(&cf_meta, b"db_type", b"spliceai")?;
    db.put_cf(&cf_meta, b"assembly", args.assembly.as_bytes())?;
    db.put_cf(&cf_meta, b"schema_version", b"1.0")?;

    tracing::info!(
        "Successfully completed SpliceAI import of {} records in {:?}",
        written,
        start_time.elapsed()
    );

    Ok(())
}
