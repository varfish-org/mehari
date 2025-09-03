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
            // Fix: For end position, we want the line containing the last base
            let end_line = if end == 0 { 0 } else { (end - 1) / index_record.line_bases };
            let num_lines = (end_line - start_line) + 1;

            let line_offset = start % index_record.line_bases;
            let line_start = start / index_record.line_bases * index_record.line_bytes;
            let offset = index_record.offset + line_start + line_offset;

            // Calculate how many bytes we need to read to get num_bases sequence bases
            // We need to account for newlines that will be filtered out
            let start_position_in_line = start % index_record.line_bases;
            let end_position_in_line = (end - 1) % index_record.line_bases;
            
            // If the range spans multiple lines, calculate total bytes including newlines
            let total_bytes_needed = if start_line == end_line {
                // Same line: just the bases we need
                num_bases
            } else {
                // Multiple lines: bases + newlines between them
                let bases_first_line = index_record.line_bases - start_position_in_line;
                let bases_last_line = end_position_in_line + 1;
                let full_lines_between = end_line - start_line - 1;
                let bases_middle_lines = full_lines_between * index_record.line_bases;
                let newlines_count = end_line - start_line; // One newline per line boundary
                
                bases_first_line + bases_middle_lines + bases_last_line + newlines_count
            };

            let start_byte = offset as usize;
            let end_byte = (offset + total_bytes_needed) as usize;
            
            // We might need to read a bit more to ensure we get all bases
            // This is a conservative approach to handle edge cases
            let end_byte_safe = std::cmp::min(end_byte + 1, self.mmap.len());
            assert!(end_byte_safe <= self.mmap.len());
            
            // Debug logging for assertion failure investigation
            let raw_bytes_len = end_byte - start_byte;
            tracing::debug!(
                "FASTA extraction debug - ac: {}, start: {:?}, end: {:?}",
                ac, start, end
            );
            tracing::debug!(
                "Index record - length: {}, line_bases: {}, line_bytes: {}, offset: {}",
                index_record.length, index_record.line_bases, index_record.line_bytes, index_record.offset
            );
            tracing::debug!(
                "Calculated values - num_bases: {}, start_line: {}, end_line: {}, num_lines: {}",
                num_bases, start_line, end_line, num_lines
            );
            tracing::debug!(
                "Byte calculations - total_bytes_needed: {}, raw_bytes_len: {}",
                total_bytes_needed, raw_bytes_len
            );
            tracing::debug!(
                "Position details - start_position_in_line: {}, end_position_in_line: {}",
                start_position_in_line, end_position_in_line
            );
            if start_line != end_line {
                let bases_first_line = index_record.line_bases - start_position_in_line;
                let bases_last_line = end_position_in_line + 1;
                let full_lines_between = end_line - start_line - 1;
                let bases_middle_lines = full_lines_between * index_record.line_bases;
                let newlines_count = end_line - start_line;
                tracing::debug!(
                    "Multi-line calculation - bases_first_line: {}, bases_middle_lines: {}, bases_last_line: {}, newlines_count: {}",
                    bases_first_line, bases_middle_lines, bases_last_line, newlines_count
                );
            }
            
            let mut seq = Vec::with_capacity(num_bases as usize);
            let mut newline_count = 0;
            let mut carriage_return_count = 0;
            let mut byte_index = start_byte;
            
            // Read bytes until we have enough sequence bases
            while seq.len() < num_bases as usize && byte_index < self.mmap.len() {
                let c = self.mmap[byte_index];
                if c == b'\n' {
                    newline_count += 1;
                } else if c == b'\r' {
                    carriage_return_count += 1;
                } else {
                    seq.push(c.to_ascii_uppercase());
                }
                byte_index += 1;
            }
            
            tracing::debug!(
                "Filtering results - seq.len(): {}, newline_count: {}, carriage_return_count: {}, bytes_read: {}",
                seq.len(), newline_count, carriage_return_count, byte_index - start_byte
            );
            
            if seq.len() != num_bases as usize {
                tracing::error!(
                    "Assertion would fail! Expected {} bases, got {} bases. Raw slice sample: {:?}",
                    num_bases,
                    seq.len(),
                    &self.mmap[start_byte..std::cmp::min(start_byte + 100, byte_index)]
                );
            }
            
            assert_eq!(seq.len(), num_bases as usize);

            Ok(Some(seq))
        } else {
            // Debug: Show available chromosome names when lookup fails
            tracing::debug!(
                "Chromosome '{}' not found. Available chromosomes: {:?}",
                ac,
                self.index.keys().collect::<Vec<_>>()
            );
            Ok(None)
        }
    }
}
