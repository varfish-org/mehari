//! Database construction and introspection tools.

use crate::pbs::txs::TxSeqDatabase;

pub mod cadd;
pub mod dbsnp;
pub(crate) mod keys;
pub mod spliceai;
pub mod transcripts;
pub mod unified;

use crate::common::contig::ContigManager;
use anyhow::{Context, Error, anyhow};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressDrawTarget, ProgressStyle};
use noodles::csi::BinningIndex;
use noodles::csi::binning_index::ReferenceSequence;
use rayon::prelude::*;
use rocksdb::{BottommostLevelCompaction, CompactOptions, DBWithThreadMode, MultiThreaded};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::time::Instant;
use thousands::Separable;

/// Reusable command-line options shared across all database builders.
#[derive(clap::Args, Debug, Clone)]
pub struct CommonPipelineArgs {
    /// Genome assembly version (e.g., "grch37" or "grch38").
    #[arg(long, required = true)]
    pub assembly: String,

    /// Path(s) to the input database file(s). All input files must share the same header/schema; mixed-header inputs are unsupported.
    #[arg(long, required = true)]
    pub input: Vec<PathBuf>,

    /// Path to the output RocksDB database directory.
    #[arg(long, required = true)]
    pub output: PathBuf,

    /// Number of records to chunk and commit per database write batch.
    #[arg(long, default_value = "100000")]
    pub batch_size: usize,

    /// Suppress progress bar rendering and status output.
    #[arg(long, short)]
    pub no_progress: bool,

    /// Number of threads to use for parallel execution (0 utilizes default Rayon pool).
    #[arg(long, short, default_value = "0")]
    pub threads: usize,
}

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
    pub no_progress: bool,
    pub threads: usize,
}

pub fn open_db_for_build(path: &Path, data_cf: &str) -> Result<rocksdb::DB, Error> {
    let options = rocksdb::Options::default();
    let options = tune_options_for_build(options, None);

    let cfs = vec!["meta", data_cf];
    Ok(rocksdb::DB::open_cf(&options, path, cfs)?)
}

pub fn open_db_for_read(
    path: &Path,
    data_cf: &str,
) -> Result<DBWithThreadMode<MultiThreaded>, Error> {
    let options = tune_options_for_read();
    let cfs = vec!["meta", data_cf];
    Ok(DBWithThreadMode::<MultiThreaded>::open_cf_for_read_only(
        &options, path, cfs, false,
    )?)
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

pub fn run_parallel_pipeline<R, M, Rec>(
    config: PipelineConfig,
    header_contig_lengths: HashMap<String, usize>,
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
    let db = open_db_for_build(config.output, config.db_type)?;

    let mut pool_builder = rayon::ThreadPoolBuilder::new();
    if config.threads > 0 {
        pool_builder = pool_builder.num_threads(config.threads);
    }
    let pool = pool_builder
        .build()
        .context("Failed to build local Rayon thread pool")?;

    let (total_records_written, global_seen_keys) = pool.install(|| -> Result<_, Error> {
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

                let max_len = header_contig_lengths
                    .get(&ref_name_str)
                    .copied()
                    .or_else(|| {
                        contig_manager
                            .get_contig_info(&ref_name_str)
                            .map(|info| info.length)
                    })
                    .unwrap_or(300_000_000);

                while begin < max_len {
                    let end = (begin + window_size).min(max_len);
                    windows.push((ref_name_str.clone(), begin, end));
                    begin += window_size;
                }
            }

            let pb = ProgressBar::new(windows.len() as u64);
            pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")?.progress_chars("█▒░"));

            // Respect either explicit quiet configuration flag or general test profile context
            if config.no_progress || cfg!(test) {
                pb.set_draw_target(ProgressDrawTarget::hidden());
            }

            match total_records {
                Some(count) => {
                    tracing::info!("Importing ~{} records across {} genomic windows for file: {:?}", count, windows.len(), input_file);
                }
                None => {
                    tracing::info!("Importing records across {} genomic windows for file: {:?}", windows.len(), input_file);
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

                    let mut batch_pair_count = 0;
                    for rec in records {
                        let (kvs, local_keys) = mapper(&rec, &contig_manager)?;
                        local_seen_keys.extend(local_keys);
                        for (key, value, _label) in kvs {
                            local_batch.put_cf(&cf_data, &key, &value);
                            batch_pair_count += 1;
                            if batch_pair_count >= config.batch_size {
                                db.write(local_batch)?;
                                local_batch = rocksdb::WriteBatch::default();
                                batch_pair_count = 0;
                            }
                        }
                        local_count += 1;
                    }
                    if batch_pair_count > 0 {
                        db.write(local_batch)?;
                    }
                    Ok((local_count, local_seen_keys))
                })
                .collect();

            for res in window_results {
                let (count, local_keys) = res?;
                total_records_written += count;
                global_seen_keys.extend(local_keys);
            }
        }
        Ok((total_records_written, global_seen_keys))
    })?;

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

    let (_, header) = open_vcf_reader(&config.input[0])?;
    let mut modified_header = header;
    if let Some(ref mut modifier) = header_modifier {
        modifier(&mut modified_header);
    }

    let mut header_contig_lengths = HashMap::new();
    for (name, contig_map) in modified_header.contigs() {
        if let Some(length) = contig_map.length() {
            header_contig_lengths.insert(name.to_string(), length);
        }
    }

    let shared_header = std::sync::Arc::new(modified_header);

    run_parallel_pipeline(
        config,
        header_contig_lengths,
        move |path, chrom, begin, end| {
            thread_local! {
                static VCF_CACHE: std::cell::RefCell<
                    Option<(PathBuf, noodles::vcf::io::IndexedReader<noodles_bgzf::io::Reader<File>>)>
                > = const { std::cell::RefCell::new(None) };
            }

            let region = format!("{}:{}-{}", chrom, begin + 1, end).parse()?;
            let mut records = Vec::new();

            VCF_CACHE.with(|cache| {
                let mut cache_borrow = cache.borrow_mut();
                let reader = match &mut *cache_borrow {
                    Some((cached_path, reader)) if cached_path == path => reader,
                    slot => {
                        let r = noodles::vcf::io::indexed_reader::Builder::default()
                            .build_from_path(path)?;
                        *slot = Some((path.to_path_buf(), r));
                        &mut slot.as_mut().unwrap().1
                    }
                };

                match reader.query(&shared_header, &region) {
                    Ok(query) => {
                        for result in query.records() {
                            let record_buf =
                                noodles::vcf::variant::RecordBuf::try_from_variant_record(
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
                Ok(())
            })?;

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
        HashMap::with_capacity(0),
        move |path, chrom, begin, end| {
            thread_local! {
                static TSV_CACHE: std::cell::RefCell<Option<(PathBuf, noodles::csi::io::IndexedReader<noodles_bgzf::io::Reader<File>, noodles::tabix::Index>)>> = const { std::cell::RefCell::new(None) };
            }

            let region = format!("{}:{}-{}", chrom, begin + 1, end).parse()?;
            let mut records = Vec::new();

            TSV_CACHE.with(|cache| {
                let mut cache_borrow = cache.borrow_mut();
                let reader = match &mut *cache_borrow {
                    Some((cached_path, reader)) if cached_path == path => reader,
                    slot => {
                        let r = noodles::tabix::io::indexed_reader::Builder::default()
                            .build_from_path(path)?;
                        *slot = Some((path.to_path_buf(), r));
                        &mut slot.as_mut().unwrap().1
                    }
                };

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
                Ok(())
            })?;

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
    finalize_db(db, &[config.db_type, "meta"], config.db_type)?;
    Ok(())
}

/// Profile optimized for rapid bulk data ingestion.
pub fn tune_options_for_build(
    options: rocksdb::Options,
    wal_dir: Option<&str>,
) -> rocksdb::Options {
    let mut options = options;

    options.create_if_missing(true);
    options.create_missing_column_families(true);
    options.prepare_for_bulk_load();

    options.set_max_background_jobs(16);
    options.set_max_subcompactions(8);
    options.increase_parallelism(8);

    options.set_write_buffer_size(1 << 30);
    options.set_target_file_size_base(1 << 30);

    if let Some(wal_dir) = wal_dir {
        options.set_wal_dir(wal_dir);
    }

    // Set lightweight compression for upper levels to protect disk bandwidth,
    // and heavy ZSTD compression for the bottom level.
    options.set_compression_per_level(&[
        rocksdb::DBCompressionType::None, // L0
        rocksdb::DBCompressionType::Lz4,  // L1
        rocksdb::DBCompressionType::Lz4,  // L2
        rocksdb::DBCompressionType::Lz4,  // L3
        rocksdb::DBCompressionType::Lz4,  // L4
        rocksdb::DBCompressionType::Lz4,  // L5
        rocksdb::DBCompressionType::Zstd, // Bottommost
    ]);

    options
}

/// Profile optimized purely for high-throughput, read-only lookups.
pub fn tune_options_for_read() -> rocksdb::Options {
    let mut options = rocksdb::Options::default();

    // Maximize threading capacity for concurrent file reading
    options.increase_parallelism(8);
    options.optimize_for_point_lookup(1 << 26); // 64MB block cache block hint

    // Configure high-efficiency partitioned indexes for your immutable layout
    let mut block_opts = rocksdb::BlockBasedOptions::default();
    block_opts.set_index_type(rocksdb::BlockBasedIndexType::TwoLevelIndexSearch);
    block_opts.set_bloom_filter(10.0, false); // 10 bits per key
    block_opts.set_partition_filters(true);
    block_opts.set_metadata_block_size(4096);
    block_opts.set_cache_index_and_filter_blocks(true);
    block_opts.set_pin_top_level_index_and_filter(true);
    block_opts.set_pin_l0_filter_and_index_blocks_in_cache(true);

    options.set_block_based_table_factory(&block_opts);
    options
}

/// Finalizes a freshly ingested WORM database by flushing memtables,
/// forcing a full multi-threaded bottommost compaction, and cleaning up WALs.
pub fn finalize_db(
    db: &DBWithThreadMode<MultiThreaded>,
    cf_names: &[&str],
    db_type: &str,
) -> anyhow::Result<()> {
    tracing::info!(
        "Starting database finalization protocol for [{}]...",
        db_type
    );
    let total_start = Instant::now();

    // 1. Proactively flush all memtables to disk.
    // This guarantees all ingested data is on disk as an SST before compaction begins.
    for cf_name in cf_names {
        if let Some(cf) = db.cf_handle(cf_name) {
            tracing::info!("Flushing memtable for Column Family: '{}'...", cf_name);
            db.flush_cf(&cf)?;
        }
    }

    // configure aggressive, exclusive manual compaction down to the bottommost level
    let mut compact_opt = CompactOptions::default();
    compact_opt.set_exclusive_manual_compaction(true);
    compact_opt.set_bottommost_level_compaction(BottommostLevelCompaction::Force);

    tracing::info!("Triggering global parallel bottommost compaction...");
    for cf_name in cf_names {
        if let Some(cf) = db.cf_handle(cf_name) {
            db.compact_range_cf_opt(&cf, None::<&[u8]>, None::<&[u8]>, &compact_opt);
        }
    }

    // poll compaction metrics until the database engine falls completely silent
    let compaction_start = Instant::now();
    let mut last_logged = compaction_start;

    loop {
        let pending = db
            .property_int_value(rocksdb::properties::COMPACTION_PENDING)?
            .unwrap_or(0);
        let running = db
            .property_int_value(rocksdb::properties::NUM_RUNNING_COMPACTIONS)?
            .unwrap_or(0);

        if pending == 0 && running == 0 {
            break;
        }

        if last_logged.elapsed() >= std::time::Duration::from_secs(5) {
            tracing::info!(
                "[{}] Still compacting... (Running workers: {}, Tasks pending: {}, Elapsed: {:?})",
                db_type,
                running,
                pending,
                compaction_start.elapsed()
            );
            last_logged = Instant::now();
        }

        std::thread::sleep(std::time::Duration::from_millis(250));
    }
    tracing::info!(
        "Compaction completed successfully in {:?}",
        compaction_start.elapsed()
    );

    // purge empty lingering WAL files to keep the directory clean
    if let Ok(entries) = std::fs::read_dir(db.path()) {
        for entry in entries.flatten() {
            if entry.path().extension() == Some(std::ffi::OsStr::new("log"))
                && let Ok(meta) = entry.metadata()
                && meta.len() == 0
            {
                let _ = std::fs::remove_file(entry.path());
            }
        }
    }

    // emit some stats
    if let Ok(Some(sst_size)) = db.property_int_value(rocksdb::properties::LIVE_SST_FILES_SIZE) {
        tracing::info!(
            "[{}] Database fully optimized. Total Live SST footprint: {} bytes. Total finalization runtime: {:?}",
            db_type,
            sst_size.separate_with_underscores(),
            total_start.elapsed()
        );
    }

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

/// Type alias for the thread-safe contig name to sequential ID lookup map.
pub type ContigIdMap = std::sync::Arc<std::sync::RwLock<rustc_hash::FxHashMap<String, u32>>>;

/// Standardizes a contig name and returns its matching interned `u32` ID.
/// If the contig hasn't been encountered yet, it registers a new sequential ID.
pub fn get_or_intern_contig(
    chrom: &str,
    contig_manager: &ContigManager,
    chrom_to_id: &ContigIdMap,
) -> (String, u32) {
    let chrom_std = contig_manager
        .get_primary_name(chrom)
        .cloned()
        .unwrap_or_else(|| chrom.to_string());

    let chrom_id = {
        let read_map = chrom_to_id.read().unwrap();
        read_map.get(&chrom_std).cloned()
    };

    let chrom_id = match chrom_id {
        Some(id) => id,
        None => {
            let mut write_map = chrom_to_id.write().unwrap();
            let next_id = write_map.len() as u32;

            if next_id >= (1 << 24) {
                panic!(
                    "Contig limit exceeded! 24-bit storage supports a maximum of 16777215 unique contigs."
                );
            }

            *write_map.entry(chrom_std.clone()).or_insert(next_id)
        }
    };

    (chrom_std, chrom_id)
}

pub fn write_contig_dictionary(
    output_path: &Path,
    db_type: &str,
    chrom_to_id: &ContigIdMap,
) -> Result<(), Error> {
    let options = rocksdb::Options::default();
    let options = tune_options_for_build(options, None);
    let db = rocksdb::DB::open_cf(&options, output_path, vec!["meta", db_type])?;
    let cf_meta = db
        .cf_handle("meta")
        .ok_or_else(|| anyhow::anyhow!("meta CF not found"))?;

    let map_guard = chrom_to_id.read().unwrap();
    let serialized_dict = serde_json::to_vec(&*map_guard)?;
    db.put_cf(&cf_meta, b"contig_dictionary", serialized_dict)?;
    finalize_db(&db, &["meta"], db_type)?;
    Ok(())
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
