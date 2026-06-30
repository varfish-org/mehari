use clap::Parser;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::time::Instant;

use crate::common::Args as CommonArgs;
use crate::common::contig::ContigManager;
use crate::pbs::seqvars::GenericLookupRecord;
use annonars::common::keys::Var;
use anyhow::{Error, anyhow};
use prost::Message;

/// Command line arguments for `db generic create` subcommand.
#[derive(Parser, Debug, Clone)]
#[command(about = "Construct generic lookup RocksDB database", long_about = None)]
pub struct Args {
    /// Assembly to use (e.g. GRCh37 or GRCh38)
    #[arg(long, required = true)]
    pub assembly: String,

    /// Path to input TSV or VCF file(s) (can be gzipped)
    #[arg(long, required = true)]
    pub input: Vec<PathBuf>,

    /// Path to output RocksDB directory
    #[arg(long, required = true)]
    pub output: PathBuf,

    /// Name of the custom database (used as prefix in INFO field names)
    #[arg(long, required = true)]
    pub db_name: String,

    /// Format of the input file ("tsv" or "vcf")
    #[arg(long, required = true)]
    pub format: String,

    /// TSV: Chromosome column header name
    #[arg(long, default_value = "chrom")]
    pub col_chrom: String,

    /// TSV: Position column header name
    #[arg(long, default_value = "pos")]
    pub col_pos: String,

    /// TSV: Reference allele column header name
    #[arg(long, default_value = "ref")]
    pub col_ref: String,

    /// TSV: Alternative allele column header name
    #[arg(long, default_value = "alt")]
    pub col_alt: String,

    /// TSV: List of value column header names to store. If omitted, all other columns are stored.
    #[arg(long)]
    pub col_values: Option<Vec<String>>,

    /// VCF: List of INFO field keys to store. If omitted, all INFO fields are stored.
    #[arg(long)]
    pub vcf_info_fields: Option<Vec<String>>,

    /// Number of rows to write in a single RocksDB batch
    #[arg(long, default_value = "100000")]
    pub batch_size: usize,
}

pub mod cli {
    pub use super::Args;
}

fn open_tsv_reader(
    path: &std::path::Path,
) -> Result<(csv::Reader<Box<dyn std::io::Read>>, csv::StringRecord), Error> {
    let file = File::open(path)?;
    let (reader, _format) = niffler::get_reader(Box::new(file))?;
    let mut buf_reader = BufReader::new(reader);

    // Skip any leading lines starting with '##'
    let mut header_line = String::new();
    let mut line = String::new();
    while buf_reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim();
        if trimmed.starts_with("##") {
            line.clear();
            continue;
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
    let headers_record = csv::StringRecord::from(headers);

    let rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .comment(Some(b'#'))
        .from_reader(Box::new(buf_reader) as Box<dyn std::io::Read>);

    Ok((rdr, headers_record))
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
    let (reader, _format) = niffler::get_reader(Box::new(file))?;
    let buf_reader = BufReader::new(reader);
    let mut reader =
        noodles::vcf::io::Reader::new(Box::new(buf_reader) as Box<dyn std::io::BufRead>);
    let header = reader.read_header()?;
    Ok((reader, header))
}

fn get_info_string(val: &noodles::vcf::variant::record_buf::info::field::Value) -> Option<String> {
    match val {
        noodles::vcf::variant::record_buf::info::field::Value::String(s) => Some(s.to_string()),
        noodles::vcf::variant::record_buf::info::field::Value::Array(
            noodles::vcf::variant::record_buf::info::field::value::Array::String(arr),
        ) => Some(arr.iter().flatten().cloned().collect::<Vec<_>>().join(",")),
        noodles::vcf::variant::record_buf::info::field::Value::Integer(i) => Some(i.to_string()),
        noodles::vcf::variant::record_buf::info::field::Value::Float(f) => Some(f.to_string()),
        noodles::vcf::variant::record_buf::info::field::Value::Flag => Some("true".to_string()),
        _ => None,
    }
}

pub fn run(_common: &CommonArgs, args: &Args) -> Result<(), Error> {
    tracing::info!(
        "Creating generic lookup RocksDB database at {:?}",
        args.output
    );
    let start_time = Instant::now();

    let contig_manager = ContigManager::new(&args.assembly);

    let mut options = rocksdb::Options::default();
    options.create_if_missing(true);
    options.create_missing_column_families(true);
    options.set_compression_type(rocksdb::DBCompressionType::Zstd);
    options.set_compression_options(-14, 19, 0, 0);

    let cfs = vec!["meta", "generic"];
    let db = rocksdb::DB::open_cf(&options, &args.output, cfs)?;

    let cf_generic = db
        .cf_handle("generic")
        .ok_or_else(|| anyhow!("generic CF not found"))?;
    let cf_meta = db
        .cf_handle("meta")
        .ok_or_else(|| anyhow!("meta CF not found"))?;

    let mut batch = rocksdb::WriteBatch::default();
    let mut count = 0;
    let mut written = 0;

    let is_vcf = args.format.to_lowercase() == "vcf";
    let mut seen_keys = HashSet::new();

    for input_file in &args.input {
        tracing::info!("Processing input file: {:?}", input_file);
        if is_vcf {
            let (mut reader, header) = open_vcf_reader(input_file)?;
            for result in reader.record_bufs(&header) {
                let record = result?;
                let chrom = record.reference_sequence_name();
                let pos = match record.variant_start() {
                    Some(start) => start.get() as i32,
                    None => continue,
                };
                let reference = record.reference_bases();
                let alternative = record.alternate_bases();

                let mut fields = HashMap::new();
                if let Some(keys) = &args.vcf_info_fields {
                    for k in keys {
                        if let Some(val) = record.info().get(k).flatten() {
                            if let Some(v_str) = get_info_string(val) {
                                fields.insert(k.clone(), v_str);
                            }
                        }
                    }
                } else {
                    for (k, val) in record.info().as_ref() {
                        if let Some(val_inner) = val {
                            if let Some(v_str) = get_info_string(val_inner) {
                                fields.insert(k.to_string(), v_str);
                            }
                        }
                    }
                }

                for k in fields.keys() {
                    seen_keys.insert(k.clone());
                }

                let chrom_std = contig_manager
                    .get_primary_name(chrom)
                    .cloned()
                    .unwrap_or_else(|| chrom.to_string());

                for alt in alternative.as_ref() {
                    let var = Var {
                        chrom: chrom_std.clone(),
                        pos,
                        reference: reference.to_string(),
                        alternative: alt.to_string(),
                    };
                    let key: Vec<u8> = var.clone().into();

                    if db.get_cf(&cf_generic, &key)?.is_some() {
                        tracing::warn!(
                            "Duplicate key found in database for variant: {}:{}{}>{}",
                            var.chrom,
                            var.pos,
                            var.reference,
                            var.alternative
                        );
                    }

                    let record_pb = GenericLookupRecord {
                        fields: fields.clone(),
                    };
                    let mut value = Vec::new();
                    record_pb.encode(&mut value)?;

                    batch.put_cf(&cf_generic, &key, &value);
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
        } else {
            // TSV
            let (mut rdr, headers_record) = open_tsv_reader(input_file)?;

            for result in rdr.records() {
                let record = result?;
                let record_map: HashMap<String, String> =
                    record.deserialize(Some(&headers_record))?;

                let chrom = record_map.get(&args.col_chrom).ok_or_else(|| {
                    anyhow!(
                        "Missing Chromosome column '{}' in {:?}",
                        args.col_chrom,
                        input_file
                    )
                })?;
                let pos: i32 = record_map
                    .get(&args.col_pos)
                    .ok_or_else(|| {
                        anyhow!(
                            "Missing Position column '{}' in {:?}",
                            args.col_pos,
                            input_file
                        )
                    })?
                    .parse()?;
                let reference = record_map.get(&args.col_ref).ok_or_else(|| {
                    anyhow!(
                        "Missing Reference column '{}' in {:?}",
                        args.col_ref,
                        input_file
                    )
                })?;
                let alternative = record_map.get(&args.col_alt).ok_or_else(|| {
                    anyhow!(
                        "Missing Alternative column '{}' in {:?}",
                        args.col_alt,
                        input_file
                    )
                })?;

                let chrom_std = contig_manager
                    .get_primary_name(chrom)
                    .cloned()
                    .unwrap_or_else(|| chrom.to_string());

                let mut fields = HashMap::new();
                if let Some(vals) = &args.col_values {
                    for v in vals {
                        if let Some(val) = record_map.get(v) {
                            fields.insert(v.clone(), val.clone());
                        } else {
                            anyhow::bail!("Value column '{}' not found in {:?}", v, input_file);
                        }
                    }
                } else {
                    for (k, v) in &record_map {
                        if k != &args.col_chrom
                            && k != &args.col_pos
                            && k != &args.col_ref
                            && k != &args.col_alt
                        {
                            fields.insert(k.clone(), v.clone());
                        }
                    }
                }

                for k in fields.keys() {
                    seen_keys.insert(k.clone());
                }

                let var = Var {
                    chrom: chrom_std,
                    pos,
                    reference: reference.to_string(),
                    alternative: alternative.to_string(),
                };
                let key: Vec<u8> = var.clone().into();

                if db.get_cf(&cf_generic, &key)?.is_some() {
                    tracing::warn!(
                        "Duplicate key found in database for variant: {}:{}{}>{}",
                        var.chrom,
                        var.pos,
                        var.reference,
                        var.alternative
                    );
                }

                let record_pb = GenericLookupRecord { fields };
                let mut value = Vec::new();
                record_pb.encode(&mut value)?;

                batch.put_cf(&cf_generic, &key, &value);
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

    db.put_cf(&cf_meta, b"db_type", b"generic")?;
    db.put_cf(&cf_meta, b"db_name", args.db_name.as_bytes())?;
    db.put_cf(&cf_meta, b"assembly", args.assembly.as_bytes())?;
    db.put_cf(&cf_meta, b"schema_version", b"1.0")?;

    let mut field_names: Vec<String> = seen_keys.into_iter().collect();
    field_names.sort();
    if !field_names.is_empty() {
        db.put_cf(&cf_meta, b"fields", field_names.join(",").as_bytes())?;
    }

    tracing::info!(
        "Successfully completed generic import of {} records in {:?}",
        written,
        start_time.elapsed()
    );

    Ok(())
}
