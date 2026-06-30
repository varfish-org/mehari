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
    #[arg(long)]
    pub col_chrom: Option<String>,

    /// TSV: Position column header name
    #[arg(long)]
    pub col_pos: Option<String>,

    /// TSV: Reference allele column header name
    #[arg(long)]
    pub col_ref: Option<String>,

    /// TSV: Alternative allele column header name
    #[arg(long)]
    pub col_alt: Option<String>,

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

struct TsvColResolver {
    idx_chrom: usize,
    idx_pos: usize,
    idx_ref: usize,
    idx_alt: usize,
    value_cols: Vec<(String, usize)>,
}

impl TsvColResolver {
    fn new(
        headers: &[String],
        col_chrom: &Option<String>,
        col_pos: &Option<String>,
        col_ref: &Option<String>,
        col_alt: &Option<String>,
        col_values: &Option<Vec<String>>,
    ) -> Result<Self, Error> {
        let clean_headers: Vec<String> = headers
            .iter()
            .map(|h| h.trim_start_matches('#').to_string())
            .collect();

        let find_idx = |spec: &Option<String>, default_name: &str| -> Result<usize, Error> {
            if let Some(s) = spec {
                if let Some(idx) = clean_headers.iter().position(|h| h.eq_ignore_ascii_case(s)) {
                    return Ok(idx);
                }
                anyhow::bail!("Column not found: {}", s);
            }
            if let Some(idx) = clean_headers
                .iter()
                .position(|h| h.eq_ignore_ascii_case(default_name))
            {
                Ok(idx)
            } else {
                anyhow::bail!("Could not resolve column for {}", default_name);
            }
        };

        let idx_chrom = find_idx(col_chrom, "chrom")?;
        let idx_pos = find_idx(col_pos, "pos")?;
        let idx_ref = find_idx(col_ref, "ref")?;
        let idx_alt = find_idx(col_alt, "alt")?;

        let mut value_cols = Vec::new();
        if let Some(vals) = col_values {
            for v in vals {
                if let Some(idx) = clean_headers.iter().position(|h| h.eq_ignore_ascii_case(v)) {
                    value_cols.push((headers[idx].clone(), idx));
                } else {
                    anyhow::bail!("Value column not found: {}", v);
                }
            }
        } else {
            for (idx, h) in headers.iter().enumerate() {
                if idx != idx_chrom && idx != idx_pos && idx != idx_ref && idx != idx_alt {
                    value_cols.push((h.clone(), idx));
                }
            }
        }

        Ok(Self {
            idx_chrom,
            idx_pos,
            idx_ref,
            idx_alt,
            value_cols,
        })
    }
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
            let (mut rdr, headers) = open_tsv_reader(input_file)?;
            let resolver = TsvColResolver::new(
                &headers,
                &args.col_chrom,
                &args.col_pos,
                &args.col_ref,
                &args.col_alt,
                &args.col_values,
            )?;

            for (name, _) in &resolver.value_cols {
                seen_keys.insert(name.clone());
            }

            for result in rdr.records() {
                let record = result?;
                if record.is_empty() || record[0].starts_with('#') {
                    continue;
                }

                if record.len() <= resolver.idx_chrom
                    || record.len() <= resolver.idx_pos
                    || record.len() <= resolver.idx_ref
                    || record.len() <= resolver.idx_alt
                {
                    continue;
                }

                let chrom = &record[resolver.idx_chrom];
                let pos: i32 = record[resolver.idx_pos].parse()?;
                let reference = &record[resolver.idx_ref];
                let alternative = &record[resolver.idx_alt];

                let chrom_std = contig_manager
                    .get_primary_name(chrom)
                    .cloned()
                    .unwrap_or_else(|| chrom.to_string());

                let mut fields = HashMap::new();
                for (name, idx) in &resolver.value_cols {
                    if *idx < record.len() {
                        fields.insert(name.clone(), record[*idx].to_string());
                    }
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
