use crate::common::contig::ContigManager;
use anyhow::anyhow;
use memmap2::Mmap;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
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
    /// Maps any valid alias (accession, chr-prefixed, etc.) to the raw FASTA ID
    alias_to_id: FxHashMap<String, String>,
    /// Maps the raw FASTA ID to the full sequence
    sequences: FxHashMap<String, Vec<u8>>,
    /// Set of aliases that correspond to the circular mitochondrial genome
    circular_contigs: FxHashSet<String>,
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

        let mut sequences = FxHashMap::default();
        let mut alias_to_id = FxHashMap::default();
        let mut circular_contigs = FxHashSet::default();

        for record_result in reference_reader.records() {
            let record = record_result?;
            let raw_id = record.id().to_string();

            // Collect all possible aliases for this FASTA sequence
            let mut aliases = vec![raw_id.clone()];
            if let Some(info) = contig_manager.get_contig_info(&raw_id) {
                aliases.push(info.name_with_chr.clone());
                aliases.push(info.name_without_chr.clone());
            }
            if let Some(acc) = contig_manager.get_accession(&raw_id) {
                aliases.push(acc.to_string());
            }

            // Map all aliases to the raw_id and check for circularity
            for alias in aliases {
                alias_to_id.insert(alias.clone(), raw_id.clone());
                if contig_manager.is_mitochondrial_alias(&alias) {
                    circular_contigs.insert(alias);
                }
            }

            sequences.insert(raw_id, record.seq().to_ascii_uppercase());
        }

        tracing::info!("...done reading reference FASTA into memory.");
        Ok(Self {
            alias_to_id,
            sequences,
            circular_contigs,
        })
    }
}

impl ReferenceReader for InMemoryFastaAccess {
    fn get(
        &self,
        ac: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> anyhow::Result<Option<Vec<u8>>> {
        let raw_id = match self.alias_to_id.get(ac) {
            Some(id) => id,
            None => return Ok(None),
        };

        if let Some(seq) = self.sequences.get(raw_id) {
            let seq_len = seq.len() as u64;

            if let Some(s) = start
                && s >= seq_len
            {
                return Err(anyhow!(
                    "Requested start ({}) is out of bounds for sequence {} of length {}",
                    s,
                    ac,
                    seq_len
                ));
            }

            let (start, end) = match (start, end) {
                (Some(start), Some(end)) => (start, end),
                (Some(start), None) => (start, seq_len),
                (None, Some(end)) => (0, end),
                (None, None) => (0, seq_len),
            };

            // Safely wrap around if the contig is circular
            if self.circular_contigs.contains(ac) && end > seq_len {
                let mut result_seq = seq[start as usize..].to_vec();
                let wrap_around_end = (end - seq_len).min(seq_len);
                if wrap_around_end > 0 {
                    let wrapped_part = &seq[..wrap_around_end as usize];
                    result_seq.extend_from_slice(wrapped_part);
                }
                Ok(Some(result_seq))
            } else {
                let clamped_end = end.min(seq_len);
                if clamped_end != end {
                    tracing::trace!(
                        "Requested sequence part {ac}:{start}-{end} is out of bounds for sequence of length {length}; clamping to {clamped_end}",
                        ac = ac,
                        start = start,
                        end = end,
                        length = seq_len,
                        clamped_end = clamped_end
                    );
                }

                if clamped_end < start {
                    Err(anyhow!(
                        "Requested sequence part {ac}:{start}-{end} is invalid (start > end)",
                        ac = ac,
                        start = start,
                        end = end
                    ))
                } else if clamped_end == start {
                    Ok(Some(Vec::new()))
                } else {
                    Ok(Some(seq[start as usize..clamped_end as usize].to_vec()))
                }
            }
        } else {
            Ok(None)
        }
    }
}

/// Memory-mapped, indexed FASTA access.
pub struct UnbufferedIndexedFastaAccess {
    #[allow(dead_code)]
    path: PathBuf,
    mmap: Mmap,
    /// Maps any valid alias (accession, chr-prefixed, etc.) to the FAI index record.
    alias_to_index: FxHashMap<String, IndexRecord>,
    /// Set of aliases that correspond to the circular mitochondrial genome
    circular_contigs: FxHashSet<String>,
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
            .map(|record| record.expect("Failed to read index record"))
            .collect();

        let mut alias_to_index = FxHashMap::default();
        let mut circular_contigs = FxHashSet::default();

        for record in index_records {
            let mut aliases = vec![record.name.clone()];
            if let Some(info) = contig_manager.get_contig_info(&record.name) {
                aliases.push(info.name_with_chr.clone());
                aliases.push(info.name_without_chr.clone());
            }
            if let Some(acc) = contig_manager.get_accession(&record.name) {
                aliases.push(acc.to_string());
            }

            for alias in aliases {
                alias_to_index.insert(alias.clone(), record.clone());
                if contig_manager.is_mitochondrial_alias(&alias) {
                    circular_contigs.insert(alias);
                }
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
            alias_to_index,
            circular_contigs,
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
        if let Some(index_record) = self.alias_to_index.get(ac) {
            if let Some(s) = start
                && s >= index_record.length
            {
                return Err(anyhow!(
                    "Requested start ({}) is out of bounds for sequence {} of length {}",
                    s,
                    ac,
                    index_record.length
                ));
            }

            let (start, end) = match (start, end) {
                (Some(start), Some(end)) => (start, end),
                (Some(start), None) => (start, index_record.length),
                (None, Some(end)) => (0, end),
                (None, None) => (0, index_record.length),
            };

            let get_linear_slice = |sub_start: u64,
                                    sub_end: u64|
             -> Result<Vec<u8>, anyhow::Error> {
                let sub_end = sub_end.min(index_record.length);
                if sub_end <= sub_start {
                    if sub_end == sub_start {
                        return Ok(Vec::new());
                    }
                    return Err(anyhow!(
                        "Requested sequence part {ac}:{start}-{end} is invalid (start > end)",
                        ac = ac,
                        start = sub_start,
                        end = sub_end
                    ));
                }
                let num_bases = sub_end - sub_start;

                let start_line = sub_start / index_record.line_bases;
                let end_line = (sub_end - 1) / index_record.line_bases;

                let num_lines_spanned = end_line - start_line + 1;
                let newline_bytes_per_line = index_record.line_bytes - index_record.line_bases;
                let num_newline_bytes = (num_lines_spanned - 1) * newline_bytes_per_line;

                let start_offset_in_line = sub_start % index_record.line_bases;
                let bytes_from_line_starts = start_line * index_record.line_bytes;
                let offset = index_record.offset + bytes_from_line_starts + start_offset_in_line;

                let read_start = offset as usize;
                let read_end = (offset + num_bases + num_newline_bytes) as usize;
                assert!(read_end <= self.mmap.len());

                // Read the byte slice and filter out any newline characters.
                Ok(self.mmap[read_start..read_end]
                    .iter()
                    .filter_map(|&c| {
                        if c != b'\n' && c != b'\r' {
                            Some(c.to_ascii_uppercase())
                        } else {
                            None
                        }
                    })
                    .collect())
            };

            if self.circular_contigs.contains(ac) && end > index_record.length {
                let mut sequence = get_linear_slice(start, index_record.length)?;
                let wrap_around_end = (end - index_record.length).min(index_record.length);
                if wrap_around_end > 0 {
                    let wrapped_part = get_linear_slice(0, wrap_around_end)?;
                    sequence.extend(wrapped_part);
                }

                return Ok(Some(sequence));
            }

            // Clamp end to sequence length (end may be larger than sequence length because of padding).
            let clamped_end = end.min(index_record.length);
            if clamped_end != end {
                tracing::trace!(
                    "Requested sequence part {ac}:{start}-{end} is out of bounds for sequence of length {length}; clamping to {clamped_end}",
                    ac = ac,
                    start = start,
                    end = end,
                    length = index_record.length,
                    clamped_end = clamped_end
                );
            }

            let seq = get_linear_slice(start, clamped_end)?;
            Ok(Some(seq))
        } else {
            Ok(None)
        }
    }
}
