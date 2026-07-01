use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use crate::common::Args as CommonArgs;
use crate::common::contig::ContigManager;
use crate::db::{DbWriter, finalize_db, get_total_records_from_tabix, open_db, open_vcf_reader};
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

    let db = open_db(&args.output, "generic")?;
    let mut writer = DbWriter::new(&db, "generic", args.batch_size)?;

    let is_vcf = args.format.to_lowercase() == "vcf";
    let mut global_seen_keys = HashSet::new();

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

        if is_vcf {
            let (mut reader, header) = open_vcf_reader(input_file)?;
            let mut chunk = Vec::with_capacity(args.batch_size);

            for result in reader.record_bufs(&header) {
                chunk.push(result?);

                if chunk.len() == args.batch_size {
                    let chunk_len = chunk.len() as u64;
                    let fields_seen = write_vcf_chunk(&mut writer, &chunk, &contig_manager, args)?;
                    global_seen_keys.extend(fields_seen);
                    pb.inc(chunk_len);
                    chunk.clear();
                }
            }
            if !chunk.is_empty() {
                let chunk_len = chunk.len() as u64;
                let fields_seen = write_vcf_chunk(&mut writer, &chunk, &contig_manager, args)?;
                global_seen_keys.extend(fields_seen);
                pb.inc(chunk_len);
            }
        } else {
            let (mut rdr, headers_record) = open_tsv_reader(input_file)?;
            let headers_arc = Arc::new(headers_record);
            let mut chunk = Vec::with_capacity(args.batch_size);

            for result in rdr.records() {
                chunk.push(result?);

                if chunk.len() == args.batch_size {
                    let chunk_len = chunk.len() as u64;
                    let fields_seen =
                        write_tsv_chunk(&mut writer, &chunk, &headers_arc, &contig_manager, args)?;
                    global_seen_keys.extend(fields_seen);
                    pb.inc(chunk_len);
                    chunk.clear();
                }
            }
            if !chunk.is_empty() {
                let chunk_len = chunk.len() as u64;
                let fields_seen =
                    write_tsv_chunk(&mut writer, &chunk, &headers_arc, &contig_manager, args)?;
                global_seen_keys.extend(fields_seen);
                pb.inc(chunk_len);
            }
        }

        pb.finish_and_clear();
    }

    let written = writer.flush()?;

    let cf_meta = db
        .cf_handle("meta")
        .ok_or_else(|| anyhow!("meta CF not found"))?;
    db.put_cf(&cf_meta, b"db_type", b"generic")?;
    db.put_cf(&cf_meta, b"db_name", args.db_name.as_bytes())?;
    db.put_cf(&cf_meta, b"assembly", args.assembly.as_bytes())?;
    db.put_cf(&cf_meta, b"schema_version", b"1.0")?;

    let mut field_names: Vec<String> = global_seen_keys.into_iter().collect();
    field_names.sort();
    if !field_names.is_empty() {
        db.put_cf(&cf_meta, b"fields", field_names.join(",").as_bytes())?;
    }

    tracing::info!(
        "Successfully completed generic import of {} records in {:?}",
        written,
        start_time.elapsed()
    );
    finalize_db(&db, &["generic", "meta"])?;

    Ok(())
}

fn write_vcf_chunk(
    writer: &mut DbWriter,
    chunk: &[noodles::vcf::variant::RecordBuf],
    contig_manager: &ContigManager,
    args: &Args,
) -> Result<HashSet<String>, Error> {
    let processed: Vec<Result<(Vec<(Vec<u8>, Vec<u8>, String)>, HashSet<String>), Error>> = chunk
        .par_iter()
        .map(|record| {
            let mut kvs = Vec::new();
            let mut local_keys = HashSet::new();

            let chrom = record.reference_sequence_name();
            let pos = match record.variant_start() {
                Some(start) => start.get() as i32,
                None => return Ok((kvs, local_keys)),
            };
            let reference = record.reference_bases();
            let alternative = record.alternate_bases();

            let mut fields = HashMap::new();
            if let Some(keys) = &args.vcf_info_fields {
                for k in keys {
                    if let Some(val) = record.info().get(k).flatten()
                        && let Some(v_str) = get_info_string(val)
                    {
                        fields.insert(k.clone(), v_str);
                    }
                }
            } else {
                for (k, val) in record.info().as_ref() {
                    if let Some(val_inner) = val
                        && let Some(v_str) = get_info_string(val_inner)
                    {
                        fields.insert(k.to_string(), v_str);
                    }
                }
            }

            for k in fields.keys() {
                local_keys.insert(k.clone());
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

                let record_pb = GenericLookupRecord {
                    fields: fields.clone(),
                };
                let mut value = Vec::new();
                record_pb.encode(&mut value)?;

                let var_label = format!(
                    "{}:{}{}>{}",
                    var.chrom, var.pos, var.reference, var.alternative
                );
                kvs.push((key, value, var_label));
            }

            Ok((kvs, local_keys))
        })
        .collect();

    let mut chunk_seen = HashSet::new();
    for res in processed {
        let (kvs, local_keys) = res?;
        chunk_seen.extend(local_keys);
        for (key, value, var_label) in kvs {
            writer.put(&key, &value, &var_label)?;
        }
    }
    Ok(chunk_seen)
}

fn write_tsv_chunk(
    writer: &mut DbWriter,
    chunk: &[csv::StringRecord],
    headers_record: &csv::StringRecord,
    contig_manager: &ContigManager,
    args: &Args,
) -> Result<HashSet<String>, Error> {
    let processed: Vec<Result<((Vec<u8>, Vec<u8>, String), HashSet<String>), Error>> = chunk
        .par_iter()
        .map(|record| {
            let record_map: HashMap<String, String> = record.deserialize(Some(headers_record))?;

            let chrom = record_map
                .get(&args.col_chrom)
                .ok_or_else(|| anyhow!("Missing Chromosome column"))?;
            let pos: i32 = record_map
                .get(&args.col_pos)
                .ok_or_else(|| anyhow!("Missing Position column"))?
                .parse()?;
            let reference = record_map
                .get(&args.col_ref)
                .ok_or_else(|| anyhow!("Missing Reference column"))?;
            let alternative = record_map
                .get(&args.col_alt)
                .ok_or_else(|| anyhow!("Missing Alternative column"))?;

            let chrom_std = contig_manager
                .get_primary_name(chrom)
                .cloned()
                .unwrap_or_else(|| chrom.to_string());

            let mut fields = HashMap::new();
            if let Some(vals) = &args.col_values {
                for v in vals {
                    let val = record_map
                        .get(v)
                        .ok_or_else(|| anyhow!("Value column not found"))?;
                    fields.insert(v.clone(), val.clone());
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

            let local_keys: HashSet<String> = fields.keys().cloned().collect();

            let var = Var {
                chrom: chrom_std,
                pos,
                reference: reference.to_string(),
                alternative: alternative.to_string(),
            };
            let key: Vec<u8> = var.clone().into();

            let record_pb = GenericLookupRecord { fields };
            let mut value = Vec::new();
            record_pb.encode(&mut value)?;

            let var_label = format!(
                "{}:{}{}>{}",
                var.chrom, var.pos, var.reference, var.alternative
            );
            Ok(((key, value, var_label), local_keys))
        })
        .collect();

    let mut chunk_seen = HashSet::new();
    for res in processed {
        let ((key, value, var_label), local_keys) = res?;
        chunk_seen.extend(local_keys);
        writer.put(&key, &value, &var_label)?;
    }
    Ok(chunk_seen)
}
