//! Database construction and introspection tools.

use crate::pbs::txs::TxSeqDatabase;

pub mod cadd;
pub mod dbsnp;
pub mod generic;
pub(crate) mod keys;
pub mod spliceai;
pub mod transcripts;

use crate::common::contig::ContigManager;
use anyhow::{Context, Error, anyhow};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressDrawTarget, ProgressStyle};
use noodles::csi::BinningIndex;
use noodles::csi::binning_index::ReferenceSequence;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
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
        noodles::vcf::io::Reader<Box<dyn BufRead>>,
        noodles::vcf::Header,
    ),
    Error,
> {
    let file = File::open(path)?;
    let (reader, _) = niffler::get_reader(Box::new(file))?;
    let mut reader =
        noodles::vcf::io::Reader::new(Box::new(BufReader::new(reader)) as Box<dyn BufRead>);
    let header = reader.read_header()?;
    Ok((reader, header))
}

pub fn open_tsv_reader(
    path: &Path,
) -> Result<(csv::Reader<Box<dyn Read>>, csv::StringRecord), Error> {
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
        .from_reader(Box::new(buf_reader) as Box<dyn Read>);

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

/// Core parallel map-reduce orchestration engine.
/// Handles Tabix space partitioning, Rayon execution, and atomic thread-local batch commits.
pub fn run_parallel_pipeline<R, M, Rec>(
    config: PipelineConfig,
    region_reader: R,
    mapper: M,
) -> Result<(), Error>
where
    Rec: Send + Sync,
    R: Fn(&Path, &str, usize, usize) -> Result<Vec<Rec>, Error> + Sync + Send,
    M: Fn(
            &Rec,
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

    let mut global_seen_keys = HashSet::new();
    let mut total_records_written = 0;
    let window_size = config.batch_size.max(500_000);

    for input_file in config.input {
        if !input_file.exists() {
            tracing::warn!("Input file does not exist, skipping: {:?}", input_file);
            continue;
        }

        let total_records = get_total_records_from_tabix(input_file)?;

        let tbi_path = input_file.with_added_extension("tbi");
        let index = noodles::tabix::fs::read(&tbi_path)
            .with_context(|| format!("Failed to read tabix index: {:?}", tbi_path.display()))?;
        let tabix_header = index
            .header()
            .ok_or_else(|| anyhow!("Missing tabix header"))?;

        let mut windows = Vec::new();
        for ref_name in tabix_header.reference_sequence_names() {
            let ref_name_str = ref_name.to_string();
            let mut begin = 0;
            let max_len = 300_000_000; // Safe default upper limit for chromosome/scaffold lengths
            while begin < max_len {
                let end = (begin + window_size).min(max_len);
                windows.push((ref_name_str.clone(), begin, end));
                begin += window_size;
            }
        }

        let pb = ProgressBar::new(windows.len() as u64);
        pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")?.progress_chars("█▒░"));
        if cfg!(test) {
            pb.set_draw_target(ProgressDrawTarget::hidden());
        }

        match total_records {
            Some(count) => {
                tracing::info!(
                    "Importing ~{} records across {} genomic windows for file: {:?}",
                    count,
                    windows.len(),
                    input_file
                );
            }
            None => {
                tracing::info!(
                    "Importing records across {} genomic windows for file: {:?}",
                    windows.len(),
                    input_file
                );
            }
        }

        let cf_data = db
            .cf_handle(config.db_type)
            .ok_or_else(|| anyhow!("CF '{}' not found", config.db_type))?;

        let window_results: Vec<Result<(usize, HashSet<String>), Error>> = windows
            .par_iter()
            .progress_with(pb)
            .map(|(chrom, begin, end)| {
                let mut local_seen_keys = HashSet::new();
                let mut local_batch = rocksdb::WriteBatch::default();
                let mut local_count = 0;

                let records = region_reader(input_file, chrom, *begin, *end)?;

                for rec in records {
                    let (kvs, local_keys) = mapper(&rec, &contig_manager)?;
                    local_seen_keys.extend(local_keys);
                    for (key, value, _label) in kvs {
                        local_batch.put_cf(&cf_data, &key, &value);
                    }
                    local_count += 1;
                }

                db.write(local_batch)?;
                Ok((local_count, local_seen_keys))
            })
            .collect();

        for res in window_results {
            let (count, local_keys) = res?;
            total_records_written += count;
            global_seen_keys.extend(local_keys);
        }
    }

    finalize_pipeline(
        &db,
        total_records_written,
        config,
        global_seen_keys,
        start_time,
    )
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
    if config.input.is_empty() {
        return Ok(());
    }

    // Extract structure schemas upfront
    let (_, header) = open_vcf_reader(&config.input[0])?;
    let mut modified_header = header;
    if let Some(ref mut modifier) = header_modifier {
        modifier(&mut modified_header);
    }
    let shared_header = std::sync::Arc::new(modified_header);

    run_parallel_pipeline(
        config,
        move |path, chrom, begin, end| {
            let mut reader =
                noodles::vcf::io::indexed_reader::Builder::default().build_from_path(path)?;
            let region = format!("{}:{}-{}", chrom, begin + 1, end).parse()?;

            let mut records = Vec::new();
            match reader.query(&shared_header, &region) {
                Ok(query) => {
                    for result in query.records() {
                        let record_buf = noodles::vcf::variant::RecordBuf::try_from_variant_record(
                            &shared_header,
                            &result?,
                        )?;
                        records.push(record_buf);
                    }
                }
                Err(e)
                    if e.to_string()
                        .contains("region reference sequence does not exist") => {}
                Err(e) => return Err(Error::from(e)),
            }
            Ok(records)
        },
        mapper,
    )
}

pub fn run_tsv_pipeline<M, R>(
    config: PipelineConfig,
    open_reader: R,
    mapper: M,
) -> Result<(), Error>
where
    R: Fn(&Path) -> Result<(csv::Reader<Box<dyn Read>>, csv::StringRecord), Error>,
    M: Fn(
            &csv::StringRecord,
            &csv::StringRecord,
            &ContigManager,
        ) -> Result<(Vec<(Vec<u8>, Vec<u8>, String)>, HashSet<String>), Error>
        + Sync
        + Send,
{
    if config.input.is_empty() {
        return Ok(());
    }

    let (_, headers_record) = open_reader(&config.input[0])?;
    let shared_headers = std::sync::Arc::new(headers_record);

    let mapper_wrapper = move |record: &csv::StringRecord, contig_manager: &ContigManager| {
        mapper(record, &shared_headers, contig_manager)
    };

    run_parallel_pipeline(
        config,
        move |path, chrom, begin, end| {
            let mut reader =
                noodles::tabix::io::indexed_reader::Builder::default().build_from_path(path)?;
            let region = format!("{}:{}-{}", chrom, begin + 1, end).parse()?;
            let mut records = Vec::new();

            match reader.query(&region) {
                Ok(query) => {
                    for line_result in query {
                        let line = line_result?;
                        let fields = line
                            .as_ref()
                            .split('\t')
                            .map(String::from)
                            .collect::<Vec<_>>();
                        records.push(csv::StringRecord::from(fields));
                    }
                }
                Err(e)
                    if e.to_string()
                        .contains("region reference sequence does not exist") => {}
                Err(e) => return Err(Error::from(e)),
            }
            Ok(records)
        },
        mapper_wrapper,
    )
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

pub fn get_total_records_from_tabix(path: &Path) -> anyhow::Result<Option<u64>> {
    let tbi_path = path.with_added_extension("tbi");

    if !tbi_path.exists() {
        tracing::warn!(path = %tbi_path.display(), "Index does not exist, skipping.");
        return Ok(None);
    }

    let index = noodles::tabix::fs::read(&tbi_path).with_context(|| {
        format!(
            "Failed to read tabix index for record counts at {:?}",
            tbi_path.display()
        )
    })?;

    let mut total_records = 0;
    let mut metadata_found = false;

    for reference in index.reference_sequences() {
        if let Some(metadata) = reference.metadata() {
            total_records += metadata.mapped_record_count();
            metadata_found = true;
        }
    }

    if metadata_found {
        Ok(Some(total_records))
    } else {
        Ok(None)
    }
}

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use noodles::core::Position;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use std::io::Write;

    /// Reusable helper to write a mock text-based database file (VCF or TSV)
    /// into valid BGZF compression alongside a compiled Tabix (.tbi) index.
    pub fn write_indexed_file(path: &Path, content: &str) -> anyhow::Result<()> {
        let mut lines = Vec::new();
        let mut chroms = std::collections::BTreeSet::new();
        let mut has_fileformat = false;
        let mut fileformat_line = String::new();

        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with("##fileformat=") {
                has_fileformat = true;
                fileformat_line = trimmed.to_string();
                continue;
            }
            if !trimmed.starts_with('#') {
                let fields: Vec<&str> = trimmed.split('\t').collect();
                if !fields.is_empty() {
                    chroms.insert(fields[0].to_string());
                }
            }
            lines.push(trimmed.to_string());
        }

        let out_file = File::create(path)?;
        let mut bgzf_writer = noodles_bgzf::io::Writer::new(out_file);

        if has_fileformat {
            bgzf_writer.write_all(fileformat_line.as_bytes())?;
            bgzf_writer.write_all(b"\n")?;
        }

        for chrom in &chroms {
            let contig_line = format!("##contig=<ID={}>", chrom);
            bgzf_writer.write_all(contig_line.as_bytes())?;
            bgzf_writer.write_all(b"\n")?;
        }

        let mut indexer = noodles::tabix::index::Indexer::default();
        indexer.set_header(noodles::csi::binning_index::index::header::Builder::vcf().build());

        let mut start_pos = bgzf_writer.virtual_position();

        for line in lines {
            bgzf_writer.write_all(line.as_bytes())?;
            bgzf_writer.write_all(b"\n")?;
            let end_pos = bgzf_writer.virtual_position();

            if line.starts_with('#') {
                start_pos = end_pos;
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 2 {
                let chrom = fields[0];
                if let Ok(pos_val) = fields[1].parse::<usize>()
                    && let Some(position) = Position::new(pos_val)
                {
                    indexer.add_record(
                        chrom,
                        position,
                        position,
                        Chunk::new(start_pos, end_pos),
                    )?;
                }
            }
            start_pos = end_pos;
        }
        bgzf_writer.finish()?;

        let index = indexer.build();
        let tbi_path = path.with_added_extension("tbi");
        let mut idx_writer = noodles::tabix::io::Writer::new(File::create(&tbi_path)?);
        idx_writer.write_index(&index)?;

        Ok(())
    }
}
