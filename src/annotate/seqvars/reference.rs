use anyhow::anyhow;
use memmap2::Mmap;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};

pub trait ReferenceReader {
    fn get(
        &self,
        ac: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> anyhow::Result<Option<Vec<u8>>>;
}

pub struct InMemoryFastaAccess {
    sequences: HashMap<String, Vec<u8>>,
}

impl InMemoryFastaAccess {
    pub fn from_path(path: impl AsRef<Path>) -> anyhow::Result<Self> {
        tracing::info!("Reading reference from {}", path.as_ref().display());
        let reference_path = path.as_ref().to_path_buf();
        let reference_reader = bio::io::fasta::Reader::from_file(&reference_path)
            .expect("Failed to create FASTA reader");
        let sequences = reference_reader
            .records()
            .map(|r| {
                let record = r.expect("Failed to read FASTA record");
                (record.id().to_string(), record.seq().to_ascii_uppercase())
            })
            .collect();
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

pub struct UnbufferedIndexedFastaAccess {
    #[allow(dead_code)]
    path: PathBuf,
    mmap: Mmap,
    index: HashMap<String, IndexRecord>,
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
    pub fn from_path(path: impl AsRef<Path>) -> anyhow::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let index_path = format!("{}.fai", path.to_str().ok_or(anyhow!("Invalid path"))?);
        tracing::info!("Reading reference index from {}", &index_path);
        let index = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(index_path)?
            .deserialize()
            .map(|record| {
                let record: IndexRecord = record.expect("Failed to read index record");
                (record.name.clone(), record)
            })
            .collect();
        let file = File::open(&path)?;
        // SAFETY: The file is opened in read-only mode;
        // however, if the underlying file is modified, this can still lead to undefined behavior.
        let mmap = unsafe { Mmap::map(&file)? };
        mmap.advise(memmap2::Advice::Sequential)?;

        Ok(Self { path, mmap, index })
    }
}

impl ReferenceReader for UnbufferedIndexedFastaAccess {
    fn get(
        &self,
        ac: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> anyhow::Result<Option<Vec<u8>>> {
        if let Some(index_record) = self.index.get(ac) {
            let (start, end) = match (start, end) {
                (Some(start), Some(end)) => (start, end),
                (Some(start), None) => (start, index_record.length),
                (None, Some(end)) => (0, end),
                (None, None) => (0, index_record.length),
            };

            assert!(start <= index_record.length);

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
