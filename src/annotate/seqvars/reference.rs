use crate::common::contig::ContigManager;
use anyhow::anyhow;
use memmap2::Mmap;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub trait ReferenceReader {
    fn get(
        &self,
        ac: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> anyhow::Result<Option<Vec<u8>>>;
}

/// In-memory reference sequence access.
pub struct InMemoryFastaAccess {
    /// Maps canonical RefSeq accession to the full sequence.
    sequences: HashMap<String, Vec<u8>>,
}
impl InMemoryFastaAccess {
    pub fn from_path(
        path: impl AsRef<Path>,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        tracing::info!(
            "Reading reference FASTA into memory from {}",
            path.as_ref().display()
        );
        let reference_path = path.as_ref().to_path_buf();
        let reference_reader = bio::io::fasta::Reader::from_file(&reference_path)
            .expect("Failed to create FASTA reader");

        let mut sequences = HashMap::new();
        for record_result in reference_reader.records() {
            let record = record_result?;
            if let Some(accession) = contig_manager.get_accession(record.id()) {
                sequences.insert(accession.clone(), record.seq().to_ascii_uppercase());
            } else {
                tracing::warn!(
                    "Contig '{}' from FASTA file not found in assembly info; it will be ignored.",
                    record.id()
                );
            }
        }

        tracing::info!("...done reading reference FASTA into memory.");
        Ok(Self { sequences })
    }
}

impl ReferenceReader for InMemoryFastaAccess {
    fn get(
        &self,
        ac: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> anyhow::Result<Option<Vec<u8>>> {
        Ok(self.sequences.get(ac).map(|seq| {
            let (start, end) = match (start, end) {
                (Some(start), Some(end)) => (start, end),
                (Some(start), None) => (start, seq.len() as u64),
                (None, Some(end)) => (0, end),
                (None, None) => (0, seq.len() as u64),
            };

            seq[start as usize..end as usize].to_vec()
        }))
    }
}

/// Memory-mapped, indexed FASTA access.
pub struct UnbufferedIndexedFastaAccess {
    #[allow(dead_code)]
    path: PathBuf,
    mmap: Mmap,
    /// Maps the canonical RefSeq accession to the FAI index record.
    accession_to_index: HashMap<String, IndexRecord>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
struct IndexRecord {
    name: String,
    length: u64,
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

impl UnbufferedIndexedFastaAccess {
    pub fn from_path(
        path: impl AsRef<Path>,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let index_path = format!("{}.fai", path.to_str().ok_or(anyhow!("Invalid path"))?);
        tracing::info!("Reading reference index from {}", &index_path);
        let index_records: Vec<IndexRecord> = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(index_path)?
            .deserialize()
            .map(|record| {
                let record: IndexRecord = record.expect("Failed to read index record");
                record
            })
            .collect();

        let mut accession_to_index = HashMap::new();
        for record in index_records {
            if let Some(accession) = contig_manager.get_accession(&record.name) {
                accession_to_index.insert(accession.clone(), record);
            } else {
                tracing::warn!(
                    "Contig '{}' from FASTA index not found in assembly info; it will be ignored.",
                    record.name
                );
            }
        }

        let file = File::open(&path)?;
        // SAFETY: The file is opened in read-only mode;
        // however, if the underlying file is modified, this can still lead to undefined behavior.
        let mmap = unsafe { Mmap::map(&file)? };
        mmap.advise(memmap2::Advice::Sequential)?;

        Ok(Self {
            path,
            mmap,
            accession_to_index,
        })
    }
}

impl ReferenceReader for UnbufferedIndexedFastaAccess {
    fn get(
        &self,
        ac: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> anyhow::Result<Option<Vec<u8>>> {
        if let Some(index_record) = self.accession_to_index.get(ac) {
            let (start, end) = match (start, end) {
                (Some(start), Some(end)) => (start, end),
                (Some(start), None) => (start, index_record.length),
                (None, Some(end)) => (0, end),
                (None, None) => (0, index_record.length),
            };

            if start > index_record.length || end > index_record.length {
                return Err(anyhow!(
                    "Requested sequence part {ac}:{start}-{end} is out of bounds for sequence of length {length}",
                    ac = ac,
                    start = start,
                    end = end,
                    length = index_record.length
                ));
            }

            let num_bases = end - start;
            let start_line = start / index_record.line_bases;
            let end_line = end / index_record.line_bases;
            let num_lines = (end_line - start_line) + 1;
            let newline_bytes = index_record.line_bytes - index_record.line_bases;
            let num_newline_bytes = (num_lines - 1) * newline_bytes;

            let line_offset = start % index_record.line_bases;
            let line_start = start / index_record.line_bases * index_record.line_bytes;
            let offset = index_record.offset + line_start + line_offset;

            let start = offset as usize;
            let end = (offset + num_bases + num_newline_bytes) as usize;
            assert!(end <= self.mmap.len());
            let mut seq = Vec::with_capacity(num_bases as usize);
            for c in self.mmap[start..end].iter().filter_map(|&c| {
                if c != b'\n' && c != b'\r' {
                    Some(c.to_ascii_uppercase())
                } else {
                    None
                }
            }) {
                seq.push(c);
            }
            assert_eq!(seq.len(), num_bases as usize);

            Ok(Some(seq))
        } else {
            Ok(None)
        }
    }
}
