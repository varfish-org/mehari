//! Implementation of `hgvs` Provider interface based on protobuf.

use crate::annotate::cli::{TranscriptPickMode, TranscriptPickType};
use crate::db::create::Reason;
use crate::db::TranscriptDatabase;
use crate::{
    annotate::seqvars::csq::ALT_ALN_METHOD,
    pbs::txs::{GeneToTxId, Strand, Transcript, TranscriptTag, TxSeqDatabase},
};
use annonars::common::cli::CANONICAL;
use anyhow::anyhow;
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use biocommons_bioutils::assemblies::{Assembly, ASSEMBLY_INFOS};
use hgvs::{
    data::error::Error,
    data::{
        cdot::json::NCBI_ALN_METHOD,
        interface::{
            Provider as ProviderInterface, TxExonsRecord, TxForRegionRecord, TxIdentityInfo,
            TxInfoRecord, TxMappingOptionsRecord,
        },
    },
    sequences::{seq_md5, TranslationTable},
};
use itertools::Itertools;
use nom::AsChar;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::SeekFrom;

#[cfg(target_os = "windows")]
use std::io::{Read, Seek};
#[cfg(not(target_os = "windows"))]
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};

/// Mitochondrial accessions.
const MITOCHONDRIAL_ACCESSIONS: &[&str] = &[
    "NC_012920.1", // rCRS
    "NC_001807.4", // CRS
];

type IntervalTree = ArrayBackedIntervalTree<i32, u32>;

pub struct TxIntervalTrees {
    /// Mapping from contig accession to index in `trees`.
    pub contig_to_idx: HashMap<String, usize>,
    /// Interval tree to index in `TxSeqDatabase::tx_db::transcripts`, for each contig.
    pub trees: Vec<IntervalTree>,
}

impl TxIntervalTrees {
    pub fn new(db: &TxSeqDatabase) -> Self {
        let (contig_to_idx, trees) = Self::build_indices(db);
        Self {
            contig_to_idx,
            trees,
        }
    }

    fn build_indices(db: &TxSeqDatabase) -> (HashMap<String, usize>, Vec<IntervalTree>) {
        let assembly = db.assembly();
        let mut contig_to_idx = HashMap::new();
        let mut trees: Vec<IntervalTree> = Vec::new();

        let mut txs = 0;

        // Pre-create interval trees for canonical contigs.
        ASSEMBLY_INFOS[assembly].sequences.iter().for_each(|seq| {
            if CANONICAL.contains(&seq.name.as_str()) {
                let contig_idx = *contig_to_idx
                    .entry(seq.refseq_ac.clone())
                    .or_insert(trees.len());
                if contig_idx >= trees.len() {
                    trees.push(IntervalTree::new());
                }
            }
        });

        for (tx_id, tx) in db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts
            .iter()
            .enumerate()
        {
            for genome_alignment in &tx.genome_alignments {
                let contig = &genome_alignment.contig;
                if let Some(contig_idx) = contig_to_idx.get(contig) {
                    let mut start = i32::MAX;
                    let mut stop = i32::MIN;
                    for exon in &genome_alignment.exons {
                        start = std::cmp::min(start, exon.alt_start_i);
                        stop = std::cmp::max(stop, exon.alt_end_i);
                    }
                    trees[*contig_idx].insert(start..stop, tx_id as u32);
                }
            }

            txs += 1;
        }

        tracing::debug!("Loaded {} transcript", txs);
        trees.iter_mut().for_each(|t| t.index());

        (contig_to_idx, trees)
    }

    pub fn get_tx_for_region(
        &self,
        tx_seq_db: &TxSeqDatabase,
        alt_ac: &str,
        _alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, Error> {
        let contig_idx = *self
            .contig_to_idx
            .get(alt_ac)
            .ok_or(Error::NoTranscriptFound(alt_ac.to_string()))?;
        let query = start_i..end_i;
        let tx_idxs = self.trees[contig_idx].find(query);

        Ok(tx_idxs
            .iter()
            .map(|entry| {
                let tx = &tx_seq_db.tx_db.as_ref().expect("no tx_db?").transcripts
                    [*entry.data() as usize];
                assert_eq!(
                    tx.genome_alignments.len(),
                    1,
                    "Can only have one alignment in Mehari"
                );
                let alt_strand = tx.genome_alignments.first().unwrap().strand;
                TxForRegionRecord {
                    tx_ac: tx.id.clone(),
                    alt_ac: alt_ac.to_string(),
                    alt_strand: match Strand::try_from(alt_strand).expect("invalid strand") {
                        Strand::Plus => 1,
                        Strand::Minus => -1,
                        _ => unreachable!("invalid strand {}", alt_strand),
                    },
                    alt_aln_method: ALT_ALN_METHOD.to_string(),
                    start_i,
                    end_i,
                }
            })
            .collect())
    }
}

/// Configuration for constructing the `Provider`.
#[derive(Debug, Clone, Default, derive_builder::Builder)]
#[builder(pattern = "immutable")]
pub struct Config {
    /// Which kind of transcript to pick / restrict to. Default is not to pick at all.
    ///
    /// Depending on `--pick-transcript-mode`, if multiple transcripts match the selection,
    /// either the first one is kept or all are kept.
    #[builder(default)]
    pub pick_transcript: Vec<TranscriptPickType>,

    /// Determines how to handle multiple transcripts. Default is to keep all.
    ///
    /// When transcript picking is enabled via `--pick-transcript`,
    /// either keep the first one found or keep all that match.
    #[builder(default)]
    pub pick_transcript_mode: TranscriptPickMode,
}

/// Provider based on the protobuf `TxSeqDatabase`.
pub struct Provider {
    /// Database of transcripts and sequences as deserialized from protobuf.
    pub tx_seq_db: TxSeqDatabase,
    /// Interval trees for the tanscripts.
    pub tx_trees: TxIntervalTrees,
    /// Mapping from gene identifier to index in `TxSeqDatabase::tx_db::gene_to_tx`.
    gene_map: HashMap<String, u32>,
    /// Mapping from transcript accession to index in `TxSeqDatabase::tx_db::transcripts`.
    tx_map: HashMap<String, u32>,
    /// Mapping from sequence accession to index in `TxSeqDatabase::seq_db::seqs`.
    seq_map: HashMap<String, u32>,
    /// When transcript picking is enabled, contains the `GeneToTxIdx` entries
    /// for each gene; the order matches the one of `tx_seq_db.gene_to_tx`.
    picked_gene_to_tx_id: Option<Vec<GeneToTxId>>,

    // reference_sequences: HashMap<String, Vec<u8>>,
    reference_reader: Option<ReferenceReader>,

    /// The data version.
    data_version: String,
    /// The schema version.
    schema_version: String,
}

fn transcript_length(tx: &Transcript) -> i32 {
    let mut max_tx_length = 0;
    for genome_alignment in tx.genome_alignments.iter() {
        // We just count length in reference so we don't have to look
        // into the CIGAR string.
        let mut tx_length = 0;
        for exon_alignment in genome_alignment.exons.iter() {
            tx_length += exon_alignment.alt_cds_end_i() - exon_alignment.alt_cds_start_i();
        }
        if tx_length > max_tx_length {
            max_tx_length = tx_length;
        }
    }
    max_tx_length
}

impl Provider {
    /// Create a new `MehariProvider` from a `TxSeqDatabase`.
    ///
    /// # Arguments
    ///
    /// * `tx_seq_db` - The `TxSeqDatabase` to use.
    pub fn new(
        mut tx_seq_db: TxSeqDatabase,
        // reference: IndexedReader<noodles::fasta::io::BufReader<File>>,
        reference: Option<impl AsRef<Path>>,
        config: Config,
    ) -> Self {
        let tx_trees = TxIntervalTrees::new(&tx_seq_db);
        let gene_map = HashMap::from_iter(
            tx_seq_db
                .tx_db
                .as_ref()
                .expect("no tx_db?")
                .gene_to_tx
                .iter()
                // .filter(|gene| !gene.filtered.unwrap_or(false))
                .enumerate()
                .map(|(idx, entry)| (entry.gene_id.clone(), idx as u32)),
        );
        let tx_map = HashMap::from_iter(
            tx_seq_db
                .tx_db
                .as_ref()
                .expect("no tx_db?")
                .transcripts
                .iter()
                // .filter(|tx| !tx.filtered.unwrap_or(false))
                .enumerate()
                .map(|(idx, tx)| (tx.id.clone(), idx as u32)),
        );
        let seq_map = HashMap::from_iter(
            tx_seq_db
                .seq_db
                .as_ref()
                .expect("no seq_db?")
                .aliases
                .iter()
                .zip(
                    tx_seq_db
                        .seq_db
                        .as_ref()
                        .expect("no seq_db?")
                        .aliases_idx
                        .iter(),
                )
                .map(|(alias, idx)| (alias.clone(), *idx)),
        );

        // When transcript picking is enabled, restrict to ManeSelect and ManePlusClinical if
        // we have any such transcript.  Otherwise, fall back to the longest transcript.
        let picked_gene_to_tx_id = if !config.pick_transcript.is_empty() {
            if let Some(tx_db) = tx_seq_db.tx_db.as_mut() {
                // The new gene-to-txid mapping we will build.
                let mut new_gene_to_tx = Vec::new();

                // Process each gene.
                for entry in tx_db.gene_to_tx.iter() {
                    // First, determine whether we have any MANE transcripts.
                    let mut longest_tx_id = None;
                    let mut tx_tags = entry
                        .tx_ids
                        .iter()
                        .enumerate()
                        .filter_map(|(i, tx_id)| {
                            tx_map.get(tx_id).map(|tx_idx| {
                                let tx = &tx_db.transcripts[*tx_idx as usize];
                                let tags = tx
                                    .tags
                                    .iter()
                                    .filter_map(|tag| {
                                        match TranscriptTag::try_from(*tag).unwrap() {
                                            TranscriptTag::Unknown | TranscriptTag::Other => None,
                                            TranscriptTag::Basic => Some(TranscriptPickType::Basic),
                                            TranscriptTag::EnsemblCanonical => {
                                                Some(TranscriptPickType::EnsemblCanonical)
                                            }
                                            TranscriptTag::ManeSelect => {
                                                Some(TranscriptPickType::ManeSelect)
                                            }
                                            TranscriptTag::ManePlusClinical => {
                                                Some(TranscriptPickType::ManePlusClinical)
                                            }
                                            TranscriptTag::RefSeqSelect => {
                                                Some(TranscriptPickType::RefSeqSelect)
                                            }
                                            TranscriptTag::Selenoprotein => None,
                                            TranscriptTag::GencodePrimary => {
                                                Some(TranscriptPickType::GencodePrimary)
                                            }
                                        }
                                    })
                                    .collect_vec();
                                let length = transcript_length(tx);
                                if let Some((prev_i, prev_length)) = longest_tx_id.as_mut() {
                                    if length > *prev_length {
                                        *prev_i = i;
                                        *prev_length = length;
                                    }
                                } else {
                                    longest_tx_id = Some((i, length));
                                }
                                (tx_id, tags, length)
                            })
                        })
                        .collect_vec();
                    if let Some((i, _)) = longest_tx_id {
                        tx_tags[i].1.push(TranscriptPickType::Length);
                    }

                    let tx_ids = match config.pick_transcript_mode {
                        TranscriptPickMode::First => {
                            // only keep the first transcript that fulfills the transcript picking strategy,
                            // if any
                            let tx_id = config
                                .pick_transcript
                                .iter()
                                .filter_map(|pick| {
                                    tx_tags
                                        .iter()
                                        .find(|(_, tags, _)| tags.contains(pick))
                                        .map(|(tx_id, _, _)| tx_id)
                                })
                                .next();
                            if let Some(tx_id) = tx_id {
                                vec![tx_id.to_string()]
                            } else {
                                vec![]
                            }
                        }
                        TranscriptPickMode::All => {
                            // keep all transcripts that fulfill the transcript picking strategy
                            tx_tags
                                .iter()
                                .filter_map(|(tx_id, tags, _)| {
                                    tags.iter()
                                        .any(|tag| config.pick_transcript.contains(tag))
                                        .then_some(tx_id.to_string())
                                })
                                .collect()
                        }
                    };

                    let new_entry = if !tx_ids.is_empty() {
                        GeneToTxId {
                            gene_id: entry.gene_id.clone(),
                            tx_ids,
                            filtered: Some(false),
                            filter_reason: None,
                        }
                    } else {
                        tracing::trace!("no transcript found for gene {} with the chosen transcript picking strategy: {:?}", &entry.gene_id, &config.pick_transcript);
                        GeneToTxId {
                            gene_id: entry.gene_id.clone(),
                            tx_ids: vec![],
                            filtered: Some(true),
                            filter_reason: Some(Reason::NoTranscriptLeft as u32),
                        }
                    };

                    tracing::trace!(
                        "picked transcripts {:?} for gene {}",
                        new_entry.tx_ids,
                        new_entry.gene_id
                    );
                    new_gene_to_tx.push(new_entry);
                }

                Some(new_gene_to_tx)
            } else {
                None
            }
        } else {
            None
        };

        // TODO obtain or construct data_version and schema_version somehow
        // for now, these are just set to a combination of things that make this provider instance
        // unique, such that `data_version` and `schema_version` can be used as keys for
        // caching.
        let data_version = format!(
            "{}{}{:?}",
            tx_seq_db.version.as_ref().unwrap_or(&"".to_string()),
            tx_seq_db
                .source_version
                .iter()
                .map(|v| format!("{:#?}", v))
                .join(","),
            config
        );
        let schema_version = data_version.clone();

        // let reference_sequences = if let Some(reference_path) = reference {
        //     let reference_path = reference_path.as_ref().to_path_buf();
        //     let reference_reader = bio::io::fasta::Reader::from_file(&reference_path)
        //         .expect("Failed to create FASTA reader");
        //     reference_reader
        //         .records()
        //         .map(|r| {
        //             let record = r.expect("Failed to read FASTA record");
        //             (record.id().to_string(), record.seq().to_ascii_uppercase())
        //         })
        //         .collect()
        // } else {
        //     HashMap::new()
        // };

        let reference_reader = reference.map(ReferenceReader::new).transpose().unwrap();

        Self {
            tx_seq_db,
            tx_trees,
            gene_map,
            tx_map,
            seq_map,
            picked_gene_to_tx_id,
            reference_reader,
            data_version,
            schema_version,
        }
    }

    /// Return the assembly of the provider.
    ///
    /// # Returns
    ///
    /// The assembly of the provider.
    pub fn assembly(&self) -> Assembly {
        self.tx_seq_db.assembly()
    }

    /// Return whether transcript picking is enabled.
    pub fn transcript_picking(&self) -> bool {
        self.picked_gene_to_tx_id.is_some()
    }

    /// Return the picked transcript IDs for a gene.
    ///
    /// # Args
    ///
    /// * `gene_name` - The gene HGNC ID.
    ///
    /// # Returns
    ///
    /// The picked transcript IDs, or None if the gene is not found.
    pub fn get_picked_transcripts(&self, hgnc_id: &str) -> Option<Vec<String>> {
        self.gene_map.get(hgnc_id).and_then(|gene_idx| {
            let gene_to_tx = if let Some(picked_gene_to_tx_id) = self.picked_gene_to_tx_id.as_ref()
            {
                picked_gene_to_tx_id
            } else {
                &self.tx_seq_db.tx_db.as_ref().expect("no tx_db?").gene_to_tx
            };

            // tracing::trace!(
            //     "get_picked_transcripts({}) = {:?}",
            //     hgnc_id,
            //     &gene_to_tx[*gene_idx as usize].tx_ids
            // );
            let gene = &gene_to_tx[*gene_idx as usize];
            if let Some(true) = gene.filtered {
                None
            } else {
                Some(gene.tx_ids.clone())
            }
        })
    }

    /// Return `Transcript` for the given transcript accession.
    ///
    /// # Args
    ///
    /// * `tx_id` - The transcript accession.
    ///
    /// # Returns
    ///
    /// The `Transcript` for the given accession, or None if the accession was not found.
    pub fn get_tx(&self, tx_id: &str) -> Option<&Transcript> {
        self.tx_map.get(tx_id).and_then(|idx| {
            let tx = &self
                .tx_seq_db
                .tx_db
                .as_ref()
                .expect("no tx_db?")
                .transcripts[*idx as usize];
            if let Some(true) = tx.filtered {
                None
            } else {
                Some(tx)
            }
        })
    }

    pub fn reference_available(&self) -> bool {
        self.reference_reader.is_some()
    }

    pub fn build_chrom_to_acc(&self, assembly: Option<Assembly>) -> HashMap<String, String> {
        let acc_to_chrom: indexmap::IndexMap<String, String> =
            self.get_assembly_map(assembly.unwrap_or_else(|| self.assembly()));
        let mut chrom_to_acc = HashMap::new();
        for (acc, chrom) in &acc_to_chrom {
            let chrom = if chrom.starts_with("chr") {
                chrom.strip_prefix("chr").unwrap()
            } else {
                chrom
            };
            chrom_to_acc.insert(chrom.to_string(), acc.clone());
            chrom_to_acc.insert(format!("chr{}", chrom), acc.clone());
        }
        chrom_to_acc
    }
}

pub(crate) struct ReferenceReader {
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

impl ReferenceReader {
    fn new(path: impl AsRef<Path>) -> anyhow::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let index_path = format!("{}.fai", path.to_str().ok_or(anyhow!("Invalid path"))?);
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

    pub(crate) fn get(
        &self,
        ac: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> anyhow::Result<Option<Vec<u8>>> {
        if let Some(index_record) = self.index.get(ac) {
            fn seek_position(idx: &IndexRecord, start: u64) -> SeekFrom {
                assert!(start <= idx.length);

                let line_offset = start % idx.line_bases;
                let line_start = start / idx.line_bases * idx.line_bytes;
                let offset = SeekFrom::Start(idx.offset + line_start + line_offset);

                offset
            }

            let (start, end) = match (start, end) {
                (Some(start), Some(end)) => (start, end),
                (Some(start), None) => (start, index_record.length),
                (None, Some(end)) => (0, end),
                (None, None) => (0, index_record.length),
            };
            let length_bases = end - start;
            let start_line = start / index_record.line_bases;
            let end_line = end / index_record.line_bases;
            let num_lines = (end_line - start_line) + 1;
            let newline_bytes = index_record.line_bytes - index_record.line_bases;
            let num_newline_bytes = (num_lines - 1) * newline_bytes;

            let mut seq_with_newlines = vec![0; (length_bases + num_newline_bytes) as usize];
            let seek_from = seek_position(index_record, start);

            #[cfg(target_os = "windows")]
            {
                let mut reader = File::open("foo")?;
                reader.seek(seek_from)?;
                reader.read_exact(seq_with_newlines)?;
            }

            #[cfg(not(target_os = "windows"))]
            {
                let reader = &self.file;
                let position = match seek_from {
                    SeekFrom::Start(pos) => pos,
                    SeekFrom::Current(_) => unreachable!(),
                    SeekFrom::End(_) => unreachable!(),
                };
                reader.read_exact_at(&mut seq_with_newlines, position)?;
            }

            fn remove_newlines_and_uppercase(data: &mut Vec<u8>) {
                let mut new_idx = 0;
                let len = data.len();
                for old_idx in 0..len {
                    let c = data[old_idx] as char;

                    if !c.is_newline() {
                        data[new_idx] = c.to_ascii_uppercase() as u8;
                        new_idx += 1;
                    }
                }

                data.truncate(new_idx);
            }

            remove_newlines_and_uppercase(&mut seq_with_newlines);
            let seq = seq_with_newlines;
            assert_eq!(seq.len(), length_bases as usize);

            Ok(Some(seq))
        } else {
            Ok(None)
        }
    }
}

impl ProviderInterface for Provider {
    fn data_version(&self) -> &str {
        // TODO replace with proper data_version (see comment in `Provider::new`)
        &self.data_version
    }

    fn schema_version(&self) -> &str {
        // TODO replace with proper schema_version (see comment in `Provider::new`)
        &self.schema_version
    }

    fn get_assembly_map(
        &self,
        assembly: biocommons_bioutils::assemblies::Assembly,
    ) -> indexmap::IndexMap<String, String> {
        indexmap::IndexMap::from_iter(
            ASSEMBLY_INFOS[assembly]
                .sequences
                .iter()
                .map(|record| (record.refseq_ac.clone(), record.name.clone())),
        )
    }

    fn get_gene_info(&self, _hgnc: &str) -> Result<hgvs::data::interface::GeneInfoRecord, Error> {
        panic!("not implemented");
    }

    fn get_pro_ac_for_tx_ac(&self, tx_ac: &str) -> Result<Option<String>, Error> {
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx_idx = tx_idx as usize;
        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx];
        Ok(tx.protein.clone())
    }

    fn get_seq_part(
        &self,
        ac: &str,
        begin: Option<usize>,
        end: Option<usize>,
    ) -> Result<String, Error> {
        // In case the accession starts with "NC" or "NT" or "NW",
        // we need to look up the sequence in the reference FASTA mapping.
        let seq = if ac.starts_with("NC")
            || ac.starts_with("NT")
            || ac.starts_with("NW") && self.reference_available()
        {
            // Requires refseq_ac, i.e. refseq reference FASTA file
            let seq = self
                .reference_reader
                .as_ref()
                .unwrap()
                .get(ac, begin.map(|x| x as u64), end.map(|x| x as u64))
                .map_err(|_| Error::NoSequenceRecord(ac.to_string()))?
                .ok_or_else(|| Error::NoSequenceRecord(ac.to_string()))?;
            return String::from_utf8(seq).map_err(|_| Error::NoSequenceRecord(ac.to_string()));
        } else {
            // Otherwise, look up the sequence in the transcript database.
            let seq_idx = *self
                .seq_map
                .get(ac)
                .ok_or_else(|| Error::NoSequenceRecord(ac.to_string()))?;
            let seq_idx = seq_idx as usize;
            self.tx_seq_db.seq_db.as_ref().expect("no seq_db?").seqs[seq_idx].as_bytes()
        };

        let slice = match (begin, end) {
            (Some(begin), Some(end)) => {
                let begin = std::cmp::min(begin, seq.len());
                let end = std::cmp::min(end, seq.len());
                &seq[begin..end]
            }
            (Some(begin), None) => {
                let begin = std::cmp::min(begin, seq.len());
                &seq[begin..]
            }
            (None, Some(end)) => &seq[..end],
            (None, None) => seq,
        };
        Ok(String::from_utf8_lossy(slice).to_string())
    }

    fn get_acs_for_protein_seq(&self, seq: &str) -> Result<Vec<String>, Error> {
        Ok(vec![format!("MD5_{}", seq_md5(seq, true)?)])
    }

    fn get_similar_transcripts(
        &self,
        _tx_ac: &str,
    ) -> Result<Vec<hgvs::data::interface::TxSimilarityRecord>, Error> {
        panic!("not implemented");
    }

    fn get_tx_exons(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        _alt_aln_method: &str,
    ) -> Result<Vec<TxExonsRecord>, Error> {
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx_idx = tx_idx as usize;

        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx];
        for genome_alignment in &tx.genome_alignments {
            if genome_alignment.contig == alt_ac {
                return Ok(genome_alignment
                    .exons
                    .iter()
                    .map(|exon| TxExonsRecord {
                        hgnc: tx.gene_id.clone(),
                        tx_ac: tx_ac.to_string(),
                        alt_ac: alt_ac.to_string(),
                        alt_aln_method: ALT_ALN_METHOD.to_string(),
                        alt_strand: match Strand::try_from(genome_alignment.strand)
                            .expect("invalid strand")
                        {
                            Strand::Plus => 1,
                            Strand::Minus => -1,
                            _ => unreachable!("invalid strand {}", &genome_alignment.strand),
                        },
                        ord: exon.ord,
                        tx_start_i: exon.alt_cds_start_i.map(|val| val - 1).unwrap_or(-1),
                        tx_end_i: exon.alt_cds_end_i.unwrap_or(-1),
                        alt_start_i: exon.alt_start_i,
                        alt_end_i: exon.alt_end_i,
                        cigar: exon.cigar.clone(),
                        tx_aseq: None,
                        alt_aseq: None,
                        tx_exon_set_id: i32::MAX,
                        alt_exon_set_id: i32::MAX,
                        tx_exon_id: i32::MAX,
                        alt_exon_id: i32::MAX,
                        exon_aln_id: i32::MAX,
                    })
                    .collect());
            }
        }

        Err(Error::NoAlignmentFound(
            tx_ac.to_string(),
            alt_ac.to_string(),
        ))
    }

    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<TxInfoRecord>, Error> {
        let tx_acs = if let Some(tx_acs) = self.get_picked_transcripts(gene) {
            tx_acs
        } else {
            tracing::warn!("no transcripts found for gene: {}", gene);
            return Ok(Vec::default());
        };

        tx_acs
            .iter()
            .filter_map(|tx_ac| -> Option<Result<TxInfoRecord, Error>> {
                if let Some(tx_idx) = self.tx_map.get(tx_ac) {
                    let tx_idx = *tx_idx as usize;
                    let tx = &self
                        .tx_seq_db
                        .tx_db
                        .as_ref()
                        .expect("no tx_db?")
                        .transcripts[tx_idx];

                    if let Some(genome_alignment) = tx.genome_alignments.first() {
                        Some(Ok(TxInfoRecord {
                            hgnc: tx.gene_id.clone(),
                            cds_start_i: genome_alignment.cds_start,
                            cds_end_i: genome_alignment.cds_end,
                            tx_ac: tx.id.clone(),
                            alt_ac: genome_alignment.contig.to_string(),
                            alt_aln_method: "splign".into(),
                        }))
                    } else {
                        Some(Err(Error::NoAlignmentFound(
                            tx_ac.to_string(),
                            format!("{:?}", self.assembly()),
                        )))
                    }
                } else {
                    tracing::warn!("transcript ID not found {} for gene {}", tx_ac, gene);
                    None
                }
            })
            .collect::<Result<Vec<_>, _>>()
    }

    fn get_tx_for_region(
        &self,
        alt_ac: &str,
        _alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, Error> {
        let contig_idx = *self
            .tx_trees
            .contig_to_idx
            .get(alt_ac)
            .ok_or(Error::NoTranscriptFound(alt_ac.to_string()))?;
        let query = start_i..end_i;
        let tx_idxs = self.tx_trees.trees[contig_idx].find(query);

        Ok(tx_idxs
            .iter()
            .map(|entry| {
                let tx = &self
                    .tx_seq_db
                    .tx_db
                    .as_ref()
                    .expect("no tx_db?")
                    .transcripts[*entry.data() as usize];
                assert_eq!(
                    tx.genome_alignments.len(),
                    1,
                    "Can only have one alignment in Mehari"
                );
                let alt_strand = tx.genome_alignments.first().unwrap().strand;
                TxForRegionRecord {
                    tx_ac: tx.id.clone(),
                    alt_ac: alt_ac.to_string(),
                    alt_strand: match Strand::try_from(alt_strand).expect("invalid strand") {
                        Strand::Plus => 1,
                        Strand::Minus => -1,
                        _ => unreachable!("invalid strand {}", alt_strand),
                    },
                    alt_aln_method: ALT_ALN_METHOD.to_string(),
                    start_i,
                    end_i,
                }
            })
            .collect())
    }

    fn get_tx_identity_info(&self, tx_ac: &str) -> Result<TxIdentityInfo, Error> {
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx_idx = tx_idx as usize;
        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx];
        let is_selenoprotein = tx.tags.contains(&(TranscriptTag::Selenoprotein as i32));

        let hgnc = tx.gene_id.clone();

        let mut tmp = tx
            .genome_alignments
            .first()
            .unwrap()
            .exons
            .iter()
            .map(|exon| {
                (
                    exon.ord,
                    exon.alt_cds_end_i
                        .map(|alt_cds_end_i| alt_cds_end_i + 1 - exon.alt_cds_start_i.unwrap())
                        .unwrap_or_default(),
                )
            })
            .collect::<Vec<(i32, i32)>>();
        tmp.sort();

        let is_mitochondrial = MITOCHONDRIAL_ACCESSIONS
            .contains(&tx.genome_alignments.first().unwrap().contig.as_str());

        let lengths = tmp.into_iter().map(|(_, length)| length).collect();
        Ok(TxIdentityInfo {
            tx_ac: tx_ac.to_string(),
            alt_ac: tx_ac.to_string(), // sic(!)
            alt_aln_method: String::from("transcript"),
            cds_start_i: tx.start_codon.unwrap_or_default(),
            cds_end_i: tx.stop_codon.unwrap_or_default(),
            lengths,
            hgnc,
            translation_table: if is_mitochondrial {
                TranslationTable::VertebrateMitochondrial
            } else if is_selenoprotein {
                TranslationTable::Selenocysteine
            } else {
                TranslationTable::Standard
            },
        })
    }

    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        _alt_aln_method: &str,
    ) -> Result<TxInfoRecord, Error> {
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx_idx = tx_idx as usize;
        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx];

        for genome_alignment in &tx.genome_alignments {
            if genome_alignment.contig == alt_ac {
                return Ok(TxInfoRecord {
                    hgnc: tx.gene_id.clone(),
                    cds_start_i: genome_alignment.cds_start,
                    cds_end_i: genome_alignment.cds_end,
                    tx_ac: tx.id.clone(),
                    alt_ac: alt_ac.to_string(),
                    alt_aln_method: ALT_ALN_METHOD.to_string(),
                });
            }
        }

        Err(Error::NoAlignmentFound(
            tx_ac.to_string(),
            alt_ac.to_string(),
        ))
    }

    fn get_tx_mapping_options(&self, tx_ac: &str) -> Result<Vec<TxMappingOptionsRecord>, Error> {
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx_idx = tx_idx as usize;

        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx];

        let genome_alignment = tx.genome_alignments.first().unwrap();
        Ok(vec![TxMappingOptionsRecord {
            tx_ac: tx_ac.to_string(),
            alt_ac: genome_alignment.contig.clone(),
            alt_aln_method: NCBI_ALN_METHOD.to_string(),
        }])
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn test_sync() {
        fn is_sync<T: Sync>() {}
        is_sync::<super::Provider>();
    }
}
