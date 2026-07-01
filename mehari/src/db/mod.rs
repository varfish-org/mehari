//! Database construction and introspection tools.

use crate::pbs::txs::TxSeqDatabase;

pub mod cadd;
pub mod dbsnp;
pub mod generic;
pub mod spliceai;
pub mod transcripts;

use anyhow::{Error, anyhow};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Trait for transcript databases.
pub trait TranscriptDatabase {
    /// Get the assembly of the transcript database.
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
    let db = rocksdb::DB::open_cf(&options, path, cfs)?;
    Ok(db)
}

/// Unified compression-aware VCF reader using niffler
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
    let (reader, _format) = niffler::get_reader(Box::new(file))?;
    let buf_reader = BufReader::new(reader);
    let mut reader =
        noodles::vcf::io::Reader::new(Box::new(buf_reader) as Box<dyn std::io::BufRead>);
    let header = reader.read_header()?;
    Ok((reader, header))
}

/// Wrapper to manage batching, duplicate safety checks, and telemetry
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

        if self.count.is_multiple_of(self.batch_size) {
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

/// Run a full range compaction across specified column families
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
