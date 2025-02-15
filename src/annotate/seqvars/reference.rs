use anyhow::anyhow;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
#[cfg(target_os = "windows")]
use std::io::{Read, Seek, SeekFrom};
#[cfg(not(target_os = "windows"))]
use std::os::unix::fs::FileExt;
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
    file: File,
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

        Ok(Self { path, file, index })
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

            let mut seq_with_newlines = vec![0; (num_bases + num_newline_bytes) as usize];

            let line_offset = start % index_record.line_bases;
            let line_start = start / index_record.line_bases * index_record.line_bytes;
            let offset = index_record.offset + line_start + line_offset;

            #[cfg(target_os = "windows")]
            {
                let mut reader = File::open(&self.path)?;
                reader.seek(SeekFrom::Start(offset))?;
                reader.read_exact(&mut seq_with_newlines)?;
            }

            #[cfg(not(target_os = "windows"))]
            {
                let reader = &self.file;
                reader.read_exact_at(&mut seq_with_newlines, offset)?;
            }

            fn remove_newlines_and_uppercase(data: &mut Vec<u8>) {
                let mut new_idx = 0;
                let len = data.len();
                for old_idx in 0..len {
                    let c = data[old_idx] as char;

                    if !(c == '\n' || c == '\r') {
                        data[new_idx] = c.to_ascii_uppercase() as u8;
                        new_idx += 1;
                    }
                }

                data.truncate(new_idx);
            }

            remove_newlines_and_uppercase(&mut seq_with_newlines);
            let seq = seq_with_newlines;
            assert_eq!(seq.len(), num_bases as usize);

            Ok(Some(seq))
        } else {
            Ok(None)
        }
    }
}
