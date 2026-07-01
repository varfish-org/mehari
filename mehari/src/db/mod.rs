//! Database construction and introspection tools.

use crate::pbs::txs::TxSeqDatabase;

pub mod cadd;
pub mod dbsnp;
pub mod generic;
pub mod spliceai;
pub mod transcripts;

use crate::common::contig::ContigManager;
use anyhow::{Error, anyhow};
use indicatif::{ProgressBar, ProgressStyle};
use noodles::csi::binning_index::ReferenceSequence;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Trait for transcript databases.
pub trait TranscriptDatabase {
    fn assembly(&self) -> String;
}

impl TranscriptDatabase for TxSeqDatabase {
    fn assembly(&self) -> String {
        let source_version = self
            .source_version
            .first()
            .expect("At least one source_version entry expected");

        // Prefer the new string field, fall back to deprecated enum field
        let assembly = if !source_version.assembly.trim().is_empty() {
            source_version.assembly.as_str()
        } else {
            // Fall back to deprecated enum field
            #[allow(deprecated)]
            match source_version.assembly_enum() {
                crate::pbs::txs::Assembly::Grch37 => "grch37",
                crate::pbs::txs::Assembly::Grch38 => "grch38",
                _ => "",
            }
        };

        match assembly {
            "grch37" | "grch37p10" => "grch37".into(),
            "grch38" => "grch38".into(),
            x => x.into(),
        }
    }
}

pub struct PipelineConfig<'a> {
    pub assembly: &'a str,
    pub input: &'a [PathBuf],
    pub output: &'a Path,
    pub batch_size: usize,
    pub db_type: &'a str,
    pub schema_version: &'a str,
    pub extra_meta: HashMap<String, String>,
}

pub fn open_db(path: &Path, data_cf: &str) -> Result<rocksdb::DB, Error> {
    let mut options = rocksdb::Options::default();
    options.create_if_missing(true);
    options.create_missing_column_families(true);

    let mut block_opts = rocksdb::BlockBasedOptions::default();
    block_opts.set_block_size(64 * 1024);
    options.set_block_based_table_factory(&block_opts);
    options.set_compression_type(rocksdb::DBCompressionType::Zstd);
    options.set_compression_options(-14, 19, 0, 16 * 1024);

    options.set_max_background_jobs(4);
    options.set_write_buffer_size(128 * 1024 * 1024);
    options.set_max_write_buffer_number(4);

    let cfs = vec!["meta", data_cf];
    Ok(rocksdb::DB::open_cf(&options, path, cfs)?)
}

pub fn open_vcf_reader(
    path: &Path,
) -> Result<
    (
        noodles::vcf::io::Reader<Box<dyn std::io::BufRead>>,
        noodles::vcf::Header,
    ),
    Error,
> {
    let file = File::open(path)?;
    let (reader, _) = niffler::get_reader(Box::new(file))?;
    let mut reader = noodles::vcf::io::Reader::new(
        Box::new(BufReader::new(reader)) as Box<dyn std::io::BufRead>
    );
    let header = reader.read_header()?;
    Ok((reader, header))
}

pub fn open_tsv_reader(
    path: &Path,
) -> Result<(csv::Reader<Box<dyn std::io::Read>>, csv::StringRecord), Error> {
    let file = File::open(path)?;
    let (reader, _) = niffler::get_reader(Box::new(file))?;
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
    let rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .comment(Some(b'#'))
        .from_reader(Box::new(buf_reader) as Box<dyn std::io::Read>);

    Ok((rdr, csv::StringRecord::from(headers)))
}

pub fn get_info_string(
    val: &noodles::vcf::variant::record_buf::info::field::Value,
) -> Option<String> {
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

pub fn run_vcf_pipeline<M, F>(
    config: PipelineConfig,
    mut header_modifier: Option<F>,
    mapper: M,
) -> Result<(), Error>
where
    F: FnMut(&mut noodles::vcf::Header),
    M: Fn(
            &noodles::vcf::variant::RecordBuf,
            &ContigManager,
        ) -> Result<(Vec<(Vec<u8>, Vec<u8>, String)>, HashSet<String>), Error>
        + Sync
        + Send,
{
    tracing::info!(
        "Creating {} RocksDB database at {:?}",
        config.db_type,
        config.output
    );
    let start_time = Instant::now();
    let contig_manager = ContigManager::new(config.assembly);
    let db = open_db(config.output, config.db_type)?;
    let mut writer = DbWriter::new(&db, config.db_type, config.batch_size)?;
    let mut global_seen_keys = HashSet::new();

    for input_file in config.input {
        tracing::info!("Processing input file: {:?}", input_file);
        let total_records = get_total_records_from_tabix(input_file).unwrap_or(0);
        let pb = if total_records > 0 {
            ProgressBar::new(total_records)
        } else {
            ProgressBar::new_spinner()
        };
        pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")?.progress_chars("█▒░"));

        let (mut reader, header) = open_vcf_reader(input_file)?;
        let mut modified_header = header;
        if let Some(ref mut modifier) = header_modifier {
            modifier(&mut modified_header);
        }

        let mut chunk = Vec::with_capacity(config.batch_size);
        for result in reader.record_bufs(&modified_header) {
            chunk.push(result?);
            if chunk.len() == config.batch_size {
                process_vcf_chunk(
                    &mut writer,
                    &chunk,
                    &contig_manager,
                    &mapper,
                    &mut global_seen_keys,
                    &pb,
                )?;
                chunk.clear();
            }
        }
        if !chunk.is_empty() {
            process_vcf_chunk(
                &mut writer,
                &chunk,
                &contig_manager,
                &mapper,
                &mut global_seen_keys,
                &pb,
            )?;
        }
        pb.finish_and_clear();
    }

    // Flush the writer here, releasing the borrow on `db`
    let written = writer.flush()?;

    // Pass `db` by reference instead of moving it
    finalize_pipeline(&db, written, config, global_seen_keys, start_time)
}

fn process_vcf_chunk<M>(
    writer: &mut DbWriter,
    chunk: &[noodles::vcf::variant::RecordBuf],
    contig_manager: &ContigManager,
    mapper: &M,
    global_seen_keys: &mut HashSet<String>,
    pb: &ProgressBar,
) -> Result<(), Error>
where
    M: Fn(
            &noodles::vcf::variant::RecordBuf,
            &ContigManager,
        ) -> Result<(Vec<(Vec<u8>, Vec<u8>, String)>, HashSet<String>), Error>
        + Sync
        + Send,
{
    let chunk_len = chunk.len() as u64;
    let processed: Vec<Result<(Vec<(Vec<u8>, Vec<u8>, String)>, HashSet<String>), Error>> = chunk
        .par_iter()
        .map(|rec| mapper(rec, contig_manager))
        .collect();
    for res in processed {
        let (kvs, local_keys) = res?;
        global_seen_keys.extend(local_keys);
        for (key, value, label) in kvs {
            writer.put(&key, &value, &label)?;
        }
    }
    pb.inc(chunk_len);
    Ok(())
}

pub fn run_tsv_pipeline<M, R>(
    config: PipelineConfig,
    open_reader: R,
    mapper: M,
) -> Result<(), Error>
where
    R: Fn(&Path) -> Result<(csv::Reader<Box<dyn std::io::Read>>, csv::StringRecord), Error>,
    M: Fn(
            &csv::StringRecord,
            &csv::StringRecord,
            &ContigManager,
        ) -> Result<(Vec<(Vec<u8>, Vec<u8>, String)>, HashSet<String>), Error>
        + Sync
        + Send,
{
    tracing::info!(
        "Creating {} RocksDB database at {:?}",
        config.db_type,
        config.output
    );
    let start_time = Instant::now();
    let contig_manager = ContigManager::new(config.assembly);
    let db = open_db(config.output, config.db_type)?;
    let mut writer = DbWriter::new(&db, config.db_type, config.batch_size)?;
    let mut global_seen_keys = HashSet::new();

    for input_file in config.input {
        tracing::info!("Processing input file: {:?}", input_file);
        let total_records = get_total_records_from_tabix(input_file).unwrap_or(0);
        let pb = if total_records > 0 {
            ProgressBar::new(total_records)
        } else {
            ProgressBar::new_spinner()
        };
        pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")?.progress_chars("█▒░"));

        let (mut rdr, headers_record) = open_reader(input_file)?;
        let headers_arc = std::sync::Arc::new(headers_record);
        let mut chunk = Vec::with_capacity(config.batch_size);

        for result in rdr.records() {
            chunk.push(result?);
            if chunk.len() == config.batch_size {
                process_tsv_chunk(
                    &mut writer,
                    &chunk,
                    &headers_arc,
                    &contig_manager,
                    &mapper,
                    &mut global_seen_keys,
                    &pb,
                )?;
                chunk.clear();
            }
        }
        if !chunk.is_empty() {
            process_tsv_chunk(
                &mut writer,
                &chunk,
                &headers_arc,
                &contig_manager,
                &mapper,
                &mut global_seen_keys,
                &pb,
            )?;
        }
        pb.finish_and_clear();
    }

    // Flush the writer here, releasing the borrow on `db`
    let written = writer.flush()?;

    // Pass `db` by reference instead of moving it
    finalize_pipeline(&db, written, config, global_seen_keys, start_time)
}

fn process_tsv_chunk<M>(
    writer: &mut DbWriter,
    chunk: &[csv::StringRecord],
    headers: &csv::StringRecord,
    contig_manager: &ContigManager,
    mapper: &M,
    global_seen_keys: &mut HashSet<String>,
    pb: &ProgressBar,
) -> Result<(), Error>
where
    M: Fn(
            &csv::StringRecord,
            &csv::StringRecord,
            &ContigManager,
        ) -> Result<(Vec<(Vec<u8>, Vec<u8>, String)>, HashSet<String>), Error>
        + Sync
        + Send,
{
    let chunk_len = chunk.len() as u64;
    let processed: Vec<Result<(Vec<(Vec<u8>, Vec<u8>, String)>, HashSet<String>), Error>> = chunk
        .par_iter()
        .map(|rec| mapper(rec, headers, contig_manager))
        .collect();
    for res in processed {
        let (kvs, local_keys) = res?;
        global_seen_keys.extend(local_keys);
        for (key, value, label) in kvs {
            writer.put(&key, &value, &label)?;
        }
    }
    pb.inc(chunk_len);
    Ok(())
}

fn finalize_pipeline(
    db: &rocksdb::DB,
    written: usize,
    config: PipelineConfig,
    global_seen_keys: HashSet<String>,
    start_time: Instant,
) -> Result<(), Error> {
    let cf_meta = db
        .cf_handle("meta")
        .ok_or_else(|| anyhow!("meta CF not found"))?;
    db.put_cf(&cf_meta, b"db_type", config.db_type.as_bytes())?;
    db.put_cf(&cf_meta, b"assembly", config.assembly.as_bytes())?;
    db.put_cf(
        &cf_meta,
        b"schema_version",
        config.schema_version.as_bytes(),
    )?;

    for (k, v) in config.extra_meta {
        db.put_cf(&cf_meta, k.as_bytes(), v.as_bytes())?;
    }
    let mut field_names: Vec<String> = global_seen_keys.into_iter().collect();
    field_names.sort();
    if !field_names.is_empty() {
        db.put_cf(&cf_meta, b"fields", field_names.join(",").as_bytes())?;
    }

    tracing::info!(
        "Successfully completed {} import of {} records in {:?}",
        config.db_type,
        written,
        start_time.elapsed()
    );
    finalize_db(db, &[config.db_type, "meta"])?;
    Ok(())
}

pub struct DbWriter<'a> {
    db: &'a rocksdb::DB,
    cf: rocksdb::ColumnFamilyRef<'a>,
    batch: rocksdb::WriteBatch,
    batch_size: usize,
    count: usize,
    written: usize,
}

impl<'a> DbWriter<'a> {
    pub fn new(db: &'a rocksdb::DB, cf_name: &str, batch_size: usize) -> Result<Self, Error> {
        let cf = db
            .cf_handle(cf_name)
            .ok_or_else(|| anyhow!("CF '{}' not found", cf_name))?;
        Ok(Self {
            db,
            cf,
            batch: rocksdb::WriteBatch::default(),
            batch_size,
            count: 0,
            written: 0,
        })
    }

    pub fn put(&mut self, key: &[u8], value: &[u8], var_label: &str) -> Result<(), Error> {
        if self.db.get_cf(&self.cf, key)?.is_some() {
            tracing::warn!("Duplicate key found in database for variant: {}", var_label);
        }
        self.batch.put_cf(&self.cf, key, value);
        self.count += 1;

        if self.count >= self.batch_size {
            let active_batch = std::mem::take(&mut self.batch);
            self.db.write(active_batch)?;
            self.written += self.count;
            tracing::info!("Imported {} records...", self.written);
            self.count = 0;
        }
        Ok(())
    }

    pub fn flush(mut self) -> Result<usize, Error> {
        if self.count > 0 {
            self.db.write(self.batch)?;
            self.written += self.count;
        }
        Ok(self.written)
    }
}

pub fn finalize_db(db: &rocksdb::DB, column_families: &[&str]) -> Result<(), Error> {
    tracing::info!("Running final database compaction...");
    for cf_name in column_families {
        let cf = db
            .cf_handle(cf_name)
            .ok_or_else(|| anyhow!("CF '{}' not found", cf_name))?;
        db.compact_range_cf(&cf, None::<&[u8]>, None::<&[u8]>);
    }
    tracing::info!("Compaction complete!");
    Ok(())
}

pub fn get_total_records_from_tabix(path: &Path) -> Option<u64> {
    let tbi_path = path.to_path_buf();
    let mut os_str = tbi_path.into_os_string();
    os_str.push(".tbi");
    let tbi_path = std::path::PathBuf::from(os_str);

    if !tbi_path.exists() {
        return None;
    }

    let mut reader = noodles::tabix::io::Reader::new(File::open(tbi_path).ok()?);
    let index = reader.read_index().ok()?;

    let mut total_records = 0;
    for reference in index.reference_sequences() {
        if let Some(metadata) = reference.metadata() {
            total_records += metadata.mapped_record_count();
        }
    }

    Some(total_records)
}
