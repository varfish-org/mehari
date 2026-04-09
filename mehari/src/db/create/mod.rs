//! Transcript database.

use crate::annotate::seqvars::ann::FeatureTag;
use crate::common::trace_rss_now;
use crate::pbs::txs::{SourceVersion, TxSeqDatabase};
use anyhow::{Error, anyhow};
use clap::Parser;
use derive_new::new;
use enumflags2::{BitFlag, BitFlags, bitflags};
use hgvs::data::cdot::json::models;
use hgvs::data::cdot::json::models::{BioType, Gene, GenomeAlignment, Tag, Transcript};
use hgvs::sequences::{TranslationTable, translate_cds};
use indexmap::IndexMap;
use itertools::Itertools;
use noodles::core::Region;
use noodles::gff::feature::record::Strand;
use noodles::gff::feature::record_buf::attributes::field::tag;
use nutype::nutype;
use once_cell::sync::Lazy;
use prost::Message;
use rayon::iter::Either;
use rayon::prelude::*;
use seqrepo::{AliasOrSeqId, Interface, SeqRepo};
use serde::Serialize;
use serde_json::json;
use serde_with::DisplayFromStr;
use serde_with::serde_as;
use std::cmp::{PartialEq, Reverse};
use std::collections::{HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::ops::{Not, RangeFull};
use std::path::Path;
use std::str::FromStr;
use std::{io::Write, path::PathBuf, time::Instant};
use strum::Display;
use thousands::Separable;

/// Mitochondrial accessions.
const MITOCHONDRIAL_ACCESSIONS: &[&str] = &["NC_012920.1", "NC_001807.4"];

/// Command line arguments for `db create txs` sub command.
#[derive(Parser, Debug)]
#[command(about = "Construct mehari transcripts and sequence database", long_about = None)]
pub struct Args {
    /// Targeted genome assembly to extract transcripts for.
    #[arg(long)]
    pub assembly: String,

    /// Version of the genome assembly, e.g. "GRCh37.p13".
    #[arg(long)]
    pub assembly_version: Option<String>,

    /// Paths to the transcript annotations to import (cdot JSON or arbitrary GFF3).
    #[arg(long, required = true)]
    pub annotation: Vec<PathBuf>,

    /// Version of annotation data (if applicable).
    #[arg(long)]
    pub annotation_version: Option<String>,

    /// Source of the transcripts. For example "RefSeq" or "Ensembl".
    #[arg(long)]
    pub transcript_source: String,

    /// Version of the transcript source. E.g. "112" for Ensembl.
    #[arg(long, required_if_eq("transcript_source", "ensembl"))]
    pub transcript_source_version: Option<String>,

    /// Path to the seqrepo instance directory to use.
    #[arg(long, required_unless_present = "transcript_sequences")]
    pub seqrepo: Option<PathBuf>,

    /// Path to FASTA file(s) containing transcript sequences (alternative to seqrepo).
    #[arg(long, required_unless_present = "seqrepo")]
    pub transcript_sequences: Option<PathBuf>,

    /// Path to TSV file for label transfer of transcripts.  Columns are
    /// transcript id (without version), (unused) gene symbol, and label.
    #[arg(long)]
    pub mane_transcripts: Option<PathBuf>,

    /// Disable rigorous filtering (useful for custom annotations).
    #[arg(long, default_value = "false")]
    pub disable_filters: bool,

    /// Number of threads to use for steps supporting parallel processing.
    #[arg(long, default_value = "1")]
    pub threads: usize,

    /// ZSTD compression level to use.
    #[arg(long, default_value = "19")]
    pub compression_level: i32,

    /// Path to output protobuf file to write to.
    #[arg(long)]
    pub output: PathBuf,
}

pub struct BasicIndexedFasta {
    fasta_path: PathBuf,
    index: noodles::fasta::fai::Index,
}

impl BasicIndexedFasta {
    pub fn new(fasta_path: PathBuf, index_path: PathBuf) -> Result<Self, Error> {
        let index = File::open(index_path)
            .map(BufReader::new)
            .map(noodles::fasta::fai::io::Reader::new)?
            .read_index()?;
        Ok(Self { fasta_path, index })
    }

    pub fn get_sequence(&self, id: &str) -> Result<Option<String>, Error> {
        let clean_id = id.strip_prefix("transcript:").unwrap_or(id);
        let clean_id_base = clean_id.split('.').next().unwrap_or(clean_id);

        let mut exact_name = None;
        for record in self.index.as_ref() {
            let rec_name = String::from_utf8_lossy(record.name());
            let rec_name_base = rec_name.split('.').next().unwrap_or(&rec_name);

            if rec_name == clean_id {
                exact_name = Some(rec_name.to_string());
                break;
            } else if exact_name.is_none() && rec_name_base == clean_id_base {
                exact_name = Some(rec_name.to_string());
            }
        }

        let target_name = match exact_name {
            Some(name) => name,
            None => {
                tracing::debug!(
                    "Sequence not found in FAI index for clean_id: '{}' (original id: '{}')",
                    clean_id,
                    id
                );
                return Ok(None);
            }
        };

        let mut reader = File::open(&self.fasta_path)
            .map(BufReader::new)
            .map(|r| noodles::fasta::io::IndexedReader::new(r, self.index.clone()))?;

        let region = Region::new::<String, RangeFull>(target_name, RangeFull);
        match reader.query(&region) {
            Ok(record) => {
                let seq = String::from_utf8(record.sequence().as_ref().to_vec())?;
                Ok(Some(seq))
            }
            Err(e) => {
                tracing::warn!("Failed to query region from indexed FASTA: {}", e);
                Ok(None)
            }
        }
    }
}

pub enum SequenceProvider {
    SeqRepo(SeqRepo),
    IndexedFasta(BasicIndexedFasta),
    FastaMap(HashMap<String, String>),
}

impl SequenceProvider {
    fn fetch_sequence(&self, alias: &AliasOrSeqId) -> Result<String, Error> {
        let id = match alias {
            AliasOrSeqId::Alias { value, .. } => value,
            AliasOrSeqId::SeqId(id) => id,
        };

        match self {
            SequenceProvider::SeqRepo(repo) => repo.fetch_sequence(alias).map_err(Into::into),
            SequenceProvider::FastaMap(map) => {
                let clean_id = id.strip_prefix("transcript:").unwrap_or(id);

                map.get(clean_id)
                    .or_else(|| map.get(clean_id.split('.').next().unwrap_or("")))
                    .cloned()
                    .ok_or_else(|| {
                        tracing::debug!(
                            "Sequence not found in FASTA map for clean_id: '{}' (original id: '{}')",
                            clean_id, id
                        );
                        anyhow!("Sequence not found in FASTA map for {}", id)
                    })
            }
            SequenceProvider::IndexedFasta(reader) => {
                if let Some(seq) = reader.get_sequence(id)? {
                    return Ok(seq);
                }

                Err(anyhow!("Sequence not found in indexed FASTA for {}", id))
            }
        }
    }
}

fn open_sequence_provider(args: &Args) -> Result<SequenceProvider, Error> {
    if let Some(seqrepo_path) = &args.seqrepo {
        return Ok(SequenceProvider::SeqRepo(open_seqrepo(seqrepo_path)?));
    }

    if let Some(fasta_path) = &args.transcript_sequences {
        let mut fai_path = fasta_path.clone();
        fai_path.as_mut_os_string().push(".fai");

        if fai_path.exists() {
            tracing::info!(
                "Opening indexed FASTA: {} (using index {})",
                fasta_path.display(),
                fai_path.display()
            );
            return Ok(SequenceProvider::IndexedFasta(BasicIndexedFasta::new(
                fasta_path.clone(),
                fai_path,
            )?));
        }

        tracing::info!(
            "Loading non-indexed FASTA into memory: {}",
            fasta_path.display()
        );
        let mut map = HashMap::new();
        let file = File::open(fasta_path)?;

        let is_gz = fasta_path
            .extension()
            .map(|ext| ext == "gz" || ext == "bgz")
            .unwrap_or(false);

        let buf_reader: Box<dyn std::io::BufRead> = if is_gz {
            Box::new(noodles_bgzf::io::Reader::new(file))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut reader = noodles::fasta::io::Reader::new(buf_reader);
        for result in reader.records() {
            let record = result?;
            let raw_name = String::from_utf8_lossy(record.name()).to_string();

            let id = raw_name
                .split_whitespace()
                .next()
                .unwrap_or(&raw_name)
                .to_string();
            let seq = String::from_utf8_lossy(record.sequence().as_ref()).to_string();

            if let Some((base, _version)) = id.split_once('.') {
                map.insert(base.to_string(), seq.clone());
            }

            map.insert(id, seq);
        }
        return Ok(SequenceProvider::FastaMap(map));
    }

    // clap guarantees this is unreachable, but we need to keep the compiler happy
    unreachable!("Either seqrepo or transcript_sequences must be set");
}

/// Helper struct for parsing the label TSV file.
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize)]
struct LabelEntry {
    /// Transcript identifier without version.
    transcript_id: String,
    /// Transcript version.
    transcript_version: usize,
    /// Gene symbol (unused).
    _gene_symbol: String,
    /// Label to transfer.
    label: String,
}
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub enum GeneId {
    Hgnc(usize),
    Fallback(String),
}

impl Display for GeneId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            GeneId::Hgnc(id) => write!(f, "HGNC:{}", id),
            GeneId::Fallback(id) => write!(f, "FALLBACK:{}", id),
        }
    }
}

impl FromStr for GeneId {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(num_str) = s.strip_prefix("HGNC:") {
            num_str.parse::<usize>().map(GeneId::Hgnc)
        } else if let Some(fallback_str) = s.strip_prefix("FALLBACK:") {
            Ok(GeneId::Fallback(fallback_str.to_string()))
        } else {
            s.parse::<usize>().map(GeneId::Hgnc)
        }
    }
}

#[nutype(
    sanitize(trim),
    validate(not_empty),
    derive(
        Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, AsRef, Deref, Borrow, Into,
        Display
    )
)]
struct TranscriptId(String);

impl TranscriptId {
    fn without_version(&self) -> Result<Self, Error> {
        Ok(Self::try_new(
            self.as_ref()
                .split('.')
                .next()
                .ok_or(anyhow!("No version to split off"))?,
        )?)
    }

    fn split_version(&self) -> (&str, u32) {
        let (ac, version) = self.rsplit_once('.').unwrap_or_else(|| {
            panic!("Invalid accession, expected format 'ac.version', got {self}")
        });
        (ac, version.parse::<u32>().expect("invalid version"))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Hash)]
#[serde(tag = "type", content = "value")]
enum Identifier {
    Gene(GeneId),
    Transcript(TranscriptId),
}
trait TranscriptExt {
    fn cds_length(&self) -> Option<u32>;

    fn protein_coding(&self) -> bool;

    fn is_on_contig(&self, contig: &str) -> bool;
}
impl TranscriptExt for Transcript {
    fn cds_length(&self) -> Option<u32> {
        self.start_codon
            .and_then(|start| self.stop_codon.map(|stop| stop.abs_diff(start)))
    }

    fn protein_coding(&self) -> bool {
        matches!(self.protein.as_ref().map(|p| !p.is_empty()), Some(true))
            || self
                .biotype
                .as_ref()
                .is_some_and(|bt| bt.iter().any(|b| b.is_protein_coding()))
    }

    fn is_on_contig(&self, contig: &str) -> bool {
        self.genome_builds.values().any(|gb| gb.contig == contig)
    }
}

trait BioTypeExt {
    fn is_protein_coding(&self) -> bool;
}

impl BioTypeExt for BioType {
    fn is_protein_coding(&self) -> bool {
        // see https://www.ensembl.org/info/genome/genebuild/biotypes.html
        // and https://www.ensembl.org/Help/Faq?id=468
        matches!(self, BioType::ProteinCoding | BioType::MRna)
    }
}

static DISCARD_BIOTYPES_TRANSCRIPTS: Lazy<HashSet<BioType>> = Lazy::new(|| {
    HashSet::from([
        // BioType::PseudogenicTranscript,
        BioType::IgCGene,
        BioType::IgDGene,
        BioType::IgJGene,
        BioType::IgVGene,
        BioType::TrCGene,
        BioType::TrDGene,
        BioType::TrJGene,
        BioType::TrVGene,
        BioType::AberrantProcessedTranscript,
        BioType::UnconfirmedTranscript,
        BioType::NmdTranscriptVariant,
    ])
});

// this used to contain `BioType::Pseudogene`,
// but we keep pseudogenes (as there are cases where a protein still is produced in certain individuals)
static DISCARD_BIOTYPES_GENES: Lazy<HashSet<BioType>> = Lazy::new(|| {
    HashSet::from([
        BioType::CGeneSegment,
        BioType::DGeneSegment,
        BioType::JGeneSegment,
        BioType::VGeneSegment,
    ])
});

#[derive(Debug, Clone, Default)]
struct TranscriptLoader {
    genome_release: String,
    annotation_version: String,
    transcript_id_to_transcript: HashMap<TranscriptId, Transcript>,
    gene_id_to_gene: HashMap<GeneId, Gene>,
    gene_id_to_transcript_ids: HashMap<GeneId, Vec<TranscriptId>>,
    discards: HashMap<Identifier, BitFlags<Reason>>,
    fixes: HashMap<Identifier, BitFlags<Fix>>,
    disable_filters: bool,
}

impl TranscriptLoader {
    fn new(genome_release: String, disable_filters: bool) -> Self {
        Self {
            genome_release,
            disable_filters,
            ..Default::default()
        }
    }

    fn merge(&mut self, other: &mut Self) -> &mut Self {
        assert_eq!(self.genome_release, other.genome_release);
        self.disable_filters |= other.disable_filters;

        // Use annotation_version from whichever provides it, preferring non-empty
        if self.annotation_version.is_empty() {
            self.annotation_version
                .clone_from(&other.annotation_version);
        }

        for (gene_id, gene) in other.gene_id_to_gene.drain() {
            if let Some(old_gene) = self.gene_id_to_gene.insert(gene_id.clone(), gene.clone()) {
                tracing::warn!(
                    "Overwriting gene: {}\n{:?},\n{:?}",
                    &gene_id,
                    &old_gene,
                    &gene
                );
            }
        }

        for (tx_id, tx) in other.transcript_id_to_transcript.drain() {
            if let Some(old_tx) = self.transcript_id_to_transcript.insert(tx_id, tx.clone()) {
                tracing::warn!(
                    "Overwriting transcript: {}\n{:?},\n{:?}",
                    old_tx.id,
                    &old_tx,
                    &tx
                );
            }
        }

        for (gene_id, mut tx_ids) in other.gene_id_to_transcript_ids.drain() {
            let ids = self.gene_id_to_transcript_ids.entry(gene_id).or_default();
            ids.append(&mut tx_ids);
            ids.sort_unstable();
            ids.dedup();
        }

        for (id, reason) in other.discards.drain() {
            *self.discards.entry(id).or_default() |= reason;
        }
        for (id, fix) in other.fixes.drain() {
            *self.fixes.entry(id).or_default() |= fix;
        }
        self
    }

    /// Load and extract from standard generic GFF3 using noodles::gff.
    fn load_gff3(&mut self, path: impl AsRef<Path>) -> Result<(), Error> {
        let file = File::open(path.as_ref())?;
        let reader: Box<dyn std::io::Read> = if path.as_ref().extension().is_some_and(|e| e == "gz")
        {
            Box::new(flate2::read::GzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let reader = BufReader::new(reader);
        let mut gff_reader = noodles::gff::io::Reader::new(reader);

        let mut tx_exons: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
        let mut tx_cds: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
        let mut tx_to_gene: HashMap<String, String> = HashMap::new();
        let mut tx_info: HashMap<String, (String, Strand)> = HashMap::new();
        let mut gene_symbols: HashMap<String, String> = HashMap::new();

        let mut raw_id_to_gene_id: HashMap<String, String> = HashMap::new();
        let mut raw_id_to_tx_id: HashMap<String, String> = HashMap::new();

        for result in gff_reader.record_bufs() {
            let record = result?;

            let contig = record.reference_sequence_name().to_string();
            let feature = record.ty().to_string();
            let strand = record.strand();

            let start = usize::from(record.start()) as i32 - 1;
            let end = usize::from(record.end()) as i32;

            let attrs = record.attributes();
            let get_attr = |key: &str| {
                attrs
                    .get(key.as_bytes())
                    .and_then(|v| v.as_string())
                    .map(|s| s.to_string())
            };

            let raw_id = get_attr(tag::ID);
            let raw_parent = get_attr(tag::PARENT);
            let name = get_attr(tag::NAME).or_else(|| get_attr("gene_name"));

            let resolve_id =
                |id: Option<String>, version: Option<String>, prefixes: &[&str]| -> String {
                    match (id, version) {
                        (Some(i), Some(v)) => format!("{i}.{v}"),
                        (Some(i), None) => i,
                        _ => {
                            let mut s = raw_id.clone().unwrap_or_default();
                            for prefix in prefixes {
                                s = s.replace(prefix, "");
                            }
                            s
                        }
                    }
                };

            match feature.as_str() {
                f if f.contains("gene") => {
                    let resolved_gene_id = resolve_id(
                        get_attr("gene_id"),
                        get_attr("version").or_else(|| get_attr("gene_version")),
                        &["gene:"],
                    );

                    if let Some(rid) = &raw_id {
                        raw_id_to_gene_id.insert(rid.clone(), resolved_gene_id.clone());
                    }

                    if !resolved_gene_id.is_empty() {
                        let fallback_id = GeneId::Fallback(resolved_gene_id.clone());
                        self.gene_id_to_gene.insert(
                            fallback_id.clone(),
                            Gene {
                                hgnc: Some(fallback_id.to_string()),
                                gene_symbol: name.clone(),
                                aliases: None,
                                biotype: None,
                                description: None,
                                map_location: None,
                                summary: None,
                                url: String::new(),
                            },
                        );
                        if let Some(n) = name {
                            gene_symbols.insert(resolved_gene_id, n);
                        }
                    }
                }
                f if f.contains("transcript") || f.contains("mRNA") || f.ends_with("RNA") => {
                    let resolved_tx_id = resolve_id(
                        get_attr("transcript_id"),
                        get_attr("version").or_else(|| get_attr("transcript_version")),
                        &["transcript:", "rna:", "rna-"],
                    );

                    if let Some(rid) = &raw_id {
                        raw_id_to_tx_id.insert(rid.clone(), resolved_tx_id.clone());
                    }

                    if !resolved_tx_id.is_empty() {
                        if let Some(p) = raw_parent {
                            let first_parent = p.split(',').next().unwrap();
                            let parent_gene = raw_id_to_gene_id
                                .get(first_parent)
                                .map(String::as_str) // Avoids cloning if present
                                .unwrap_or(first_parent)
                                .to_string();
                            tx_to_gene.insert(resolved_tx_id.clone(), parent_gene);
                        }
                        tx_info.insert(resolved_tx_id, (contig, strand));
                    }
                }
                "exon" | "CDS" => {
                    if let Some(p) = raw_parent {
                        let target_map = if feature == "exon" {
                            &mut tx_exons
                        } else {
                            &mut tx_cds
                        };

                        for parent_id in p.split(',') {
                            let clean_parent = raw_id_to_tx_id
                                .get(parent_id)
                                .map(String::as_str)
                                .unwrap_or(parent_id)
                                .to_string();

                            target_map
                                .entry(clean_parent)
                                .or_default()
                                .push((start, end));
                        }
                    }
                }
                _ => {}
            }
        }

        // Finalize transcripts by resolving genomic-to-transcript coordinates
        for (tx_id, (contig, gff_strand)) in tx_info {
            let mut exons = tx_exons.remove(&tx_id).unwrap_or_default();
            let cds_fragments = tx_cds.remove(&tx_id).unwrap_or_default();

            if exons.is_empty() {
                continue;
            }

            // Sort exons by genomic position
            exons.sort_by_key(|e| e.0);

            let is_reverse = matches!(gff_strand, Strand::Reverse);
            if is_reverse {
                exons.reverse();
            }

            let cds_start_genomic = cds_fragments.iter().map(|c| c.0).min();
            let cds_end_genomic = cds_fragments.iter().map(|c| c.1).max();

            let tx_strand = if is_reverse {
                models::Strand::Minus
            } else {
                models::Strand::Plus
            };

            let mut current_tx_pos = 0;
            let mut tx_cds_start = None;
            let mut tx_cds_end = None;

            let final_exons: Vec<_> = exons
                .into_iter()
                .enumerate()
                .map(|(i, (start, end))| {
                    let e_len = end - start;
                    let mut e_cds_start = -1;
                    let mut e_cds_end = -1;

                    if let (Some(cs), Some(ce)) = (cds_start_genomic, cds_end_genomic) {
                        let overlap_start = cs.max(start);
                        let overlap_end = ce.min(end);

                        if overlap_start < overlap_end {
                            e_cds_start = overlap_start;
                            e_cds_end = overlap_end;

                            // Calculate offset within this exon based on strand
                            let (offset_start, offset_end) = if !is_reverse {
                                (overlap_start - start, overlap_end - start)
                            } else {
                                (end - overlap_end, end - overlap_start)
                            };

                            if tx_cds_start.is_none() {
                                tx_cds_start = Some((current_tx_pos + offset_start) as u32);
                            }
                            tx_cds_end = Some((current_tx_pos + offset_end) as u32);
                        }
                    }

                    let exon_record = models::Exon {
                        alt_start_i: start,
                        alt_end_i: end,
                        ord: i as i32,
                        alt_cds_start_i: e_cds_start,
                        alt_cds_end_i: e_cds_end,
                        cigar: format!("{}M", e_len),
                    };

                    current_tx_pos += e_len;
                    exon_record
                })
                .collect();

            let alignment = GenomeAlignment {
                contig,
                strand: tx_strand,
                cds_start: cds_start_genomic,
                cds_end: cds_end_genomic,
                exons: final_exons,
                tag: None,
                note: None,
            };

            let gene_ref = tx_to_gene
                .get(&tx_id)
                .cloned()
                .unwrap_or_else(|| tx_id.clone());
            let gene_name = gene_symbols.get(&gene_ref).cloned();
            let fake_gene_id = GeneId::Fallback(gene_ref.clone());

            let transcript = Transcript {
                id: tx_id.clone(),
                hgnc: Some(fake_gene_id.to_string()),
                gene_name,
                gene_version: "".to_string(),
                biotype: None,
                protein: tx_cds_start.map(|_| "unspecified_protein".to_string()),
                start_codon: tx_cds_start.map(i32::try_from).transpose()?,
                stop_codon: tx_cds_end.map(i32::try_from).transpose()?,
                partial: None,
                genome_builds: IndexMap::from([(self.genome_release.clone(), alignment)]),
            };

            let t_id = TranscriptId::try_new(tx_id)?;
            self.transcript_id_to_transcript
                .insert(t_id.clone(), transcript);
            self.gene_id_to_transcript_ids
                .entry(fake_gene_id.clone())
                .or_default()
                .push(t_id);

            // Ensure the gene entry exists even if no explicit 'gene' feature was in GFF
            self.gene_id_to_gene
                .entry(fake_gene_id)
                .or_insert_with(|| Gene {
                    hgnc: Some(gene_ref),
                    gene_symbol: None,
                    aliases: None,
                    biotype: None,
                    description: None,
                    map_location: None,
                    summary: None,
                    url: "".into(),
                });
        }
        Ok(())
    }

    /// Load and extract from cdot JSON.
    fn load_cdot(&mut self, path: impl AsRef<Path>) -> Result<(), Error> {
        let models::Container {
            genes: cdot_genes,
            transcripts: cdot_transcripts,
            cdot_version,
            ..
        } = read_cdot_json(path.as_ref())?;
        let cdot_genes = cdot_genes.into_iter().collect::<HashMap<_, _>>();
        let cdot_transcripts = cdot_transcripts
            .into_iter()
            .map(|(txid, mut tx)| {
                let txid = if txid.starts_with("fake-rna-") && !txid.contains('.') {
                    format!("{txid}.0")
                } else {
                    txid
                };
                tx.id.clone_from(&txid);
                TranscriptId::try_new(txid).map(|t| (t, tx))
            })
            .collect::<Result<HashMap<_, _>, _>>()?;
        self.annotation_version = cdot_version;

        for (_, gene) in cdot_genes {
            if let Some(gene_id) = gene.hgnc.as_ref() {
                let gene_id: GeneId = gene_id.parse()?;
                self.gene_id_to_gene.insert(gene_id.clone(), gene);
                self.gene_id_to_transcript_ids.entry(gene_id).or_default();
            }
        }

        for (tx_id, mut tx) in cdot_transcripts {
            let gene_id = if let Some(hgnc_str) = tx.hgnc.as_ref() {
                hgnc_str.parse().ok()
            } else {
                None
            };

            let gene_id = match gene_id {
                Some(id) => id,
                None => {
                    // group by gene_name if available, otherwise fall back to transcript id
                    let grouping_key = if !tx.gene_version.is_empty() {
                        tx.gene_version.as_str()
                    } else {
                        tx.gene_name.as_deref().unwrap_or(tx.id.as_str())
                    };
                    let fake_id = GeneId::Fallback(grouping_key.to_string());

                    // only insert the dummy gene if we haven't seen this PseudoId before
                    self.gene_id_to_gene
                        .entry(fake_id.clone())
                        .or_insert_with(|| Gene {
                            hgnc: Some(fake_id.clone().to_string()),
                            gene_symbol: tx.gene_name.clone(),
                            aliases: None,
                            biotype: None,
                            description: None,
                            map_location: None,
                            summary: None,
                            url: "".into(),
                        });
                    tx.hgnc = Some(fake_id.clone().to_string());

                    fake_id
                }
            };

            self.gene_id_to_transcript_ids
                .entry(gene_id)
                .or_default()
                .push(tx_id.clone());
            self.transcript_id_to_transcript.insert(tx_id, tx);
        }

        Ok(())
    }

    fn filter_initial_hgnc_entries(&mut self) -> Result<(), Error> {
        if self.disable_filters {
            return Ok(());
        }
        tracing::info!("Filtering Hgnc entries …");
        for (gene_id, tx_ids) in &self.gene_id_to_transcript_ids {
            if tx_ids.is_empty() {
                *self
                    .discards
                    .entry(Identifier::Gene(gene_id.clone()))
                    .or_default() |= Reason::NoTranscripts;
            }
            if let Some(gene) = self.gene_id_to_gene.get(gene_id) {
                if gene
                    .biotype
                    .as_ref()
                    .is_some_and(|bt| bt.iter().any(|b| DISCARD_BIOTYPES_GENES.contains(b)))
                {
                    *self
                        .discards
                        .entry(Identifier::Gene(gene_id.clone()))
                        .or_default() |= Reason::Biotype;
                }
            } else {
                *self
                    .discards
                    .entry(Identifier::Gene(gene_id.clone()))
                    .or_default() |= Reason::MissingGene;
            }
        }
        // self.discard()?;
        tracing::info!("… done filter Hgnc entries");
        Ok(())
    }

    fn apply_fixes(
        &mut self,
        transcript_id_to_tags: &Option<HashMap<TranscriptId, Vec<FeatureTag>>>,
    ) {
        self.fix_transcript_genome_builds();
        if let Some(txid_to_label) = transcript_id_to_tags {
            self.update_transcript_tags(txid_to_label);
        }
        self.fix_cds();
    }

    fn filter_genes(&mut self) -> Result<(), Error> {
        if self.disable_filters {
            return Ok(());
        }
        tracing::info!("Filtering genes …");
        // Check whether the gene has a gene symbol.
        let missing_symbol = |_gene_id: GeneId, gene: &Gene, _txs: &[&Transcript]| -> bool {
            gene.gene_symbol.is_none() || gene.gene_symbol.as_ref().unwrap().is_empty()
        };
        // Check whether the gene has a biotype that should be discarded.
        let biotype = |_gene_id: GeneId, gene: &Gene, _txs: &[&Transcript]| -> bool {
            gene.biotype
                .as_ref()
                .map(|bt| bt.iter().any(|bt| DISCARD_BIOTYPES_GENES.contains(bt)))
                .unwrap_or(false)
        };
        // Check whether all transcripts for a gene are predicted.
        let predicted_only = |_gene_id: GeneId, _gene: &Gene, txs: &[&Transcript]| -> bool {
            !txs.is_empty() && txs.iter().all(|tx| tx.id.starts_with('X'))
        };
        type Filter = fn(GeneId, &Gene, &[&Transcript]) -> bool;
        let filters: [(Filter, Reason); 3] = [
            (missing_symbol, Reason::MissingGeneSymbol),
            (biotype, Reason::Biotype),
            (predicted_only, Reason::PredictedTranscriptsOnly),
        ];

        let discarded_genes = self
            .gene_id_to_gene
            .iter()
            .filter_map(|(gene_id, gene)| {
                let txs = self
                    .gene_id_to_transcript_ids
                    .get(gene_id)
                    .map(|tx_ids| {
                        tx_ids
                            .iter()
                            .filter_map(|tx_id| self.transcript_id_to_transcript.get(tx_id))
                            .collect_vec()
                    })
                    .unwrap_or_default();
                let reason = filters
                    .iter()
                    .filter_map(|(f, r)| f(gene_id.clone(), gene, &txs).then_some(*r))
                    .fold(BitFlags::<Reason>::default(), |a, b| a | b);
                let empty = reason.is_empty();
                (!empty).then_some((Identifier::Gene(gene_id.clone()), reason))
            })
            .collect_vec();

        for (id, reason) in discarded_genes {
            self.mark_discarded(&id, reason)?;
        }
        // self.discard()?;
        tracing::info!("… done filtering genes");
        Ok(())
    }

    /// Filter transcripts for gene.
    ///
    /// We employ the following rules:
    ///
    /// - Remove redundant transcripts with the same identifier and pick only the
    ///   transcripts that have the highest version number for one assembly.
    /// - Do not pick any `XM_`/`XR_` (NCBI predicted only) transcripts.
    /// - Do not pick any `NR_` transcripts when there are coding `NM_` transcripts.
    fn filter_transcripts(&mut self) -> Result<(), Error> {
        if self.disable_filters {
            return Ok(());
        }
        tracing::info!("Filtering transcripts …");

        /// Params used as args for filter functions defined below
        #[derive(new)]
        struct Params<'a> {
            tx: &'a Transcript,
            tx_id: &'a str,
        }

        // Check whether the transcript has an HGNC id.
        let missing_hgnc = |p: &Params| -> bool {
            p.tx.hgnc.is_none()
                || p.tx.hgnc.as_ref().unwrap().is_empty()
                || p.tx.hgnc.as_ref().unwrap().starts_with("FALLBACK:")
        };

        // Check whether the transcript has any genome builds.
        let empty_genome_builds = |p: &Params| -> bool { p.tx.genome_builds.is_empty() };

        // Check whether the transcript has tags cds_start_NF or cds_end_NF.
        let cds_start_end_nf = |p: &Params| -> bool {
            p.tx.genome_builds.values().any(|gb| {
                if let Some(tag) = &gb.tag {
                    tag.contains(&Tag::Other("cds_start_NF".into()))
                        || tag.contains(&Tag::Other("cds_end_NF".into()))
                } else {
                    false
                }
            })
        };

        // Check whether the transcript is marked as partial.
        // Ignore "partial" tag for NR transcripts.
        // For protein coding transcripts which are marked as partial,
        // check whether there's at least 1bp of 5' and 3' UTR
        // (and then discard the 'partial' tag).
        let partial = |p: &Params| -> bool {
            if !matches!(p.tx.partial, Some(1)) {
                return false;
            }
            if p.tx_id.starts_with("NR") {
                false
            } else if p.tx.protein_coding() {
                assert_eq!(p.tx.genome_builds.len(), 1);
                let alignment = p.tx.genome_builds.values().next().unwrap();
                let exons = alignment.exons.clone();
                if exons.is_empty() {
                    return true;
                }
                let five_prime_trunc = {
                    let cds_start = alignment.cds_start;
                    let first = exons.iter().min_by_key(|e| e.alt_start_i).unwrap();
                    Some(first.alt_start_i) == cds_start
                };
                let three_prime_trunc = {
                    let cds_end = alignment.cds_end;
                    let last = exons.iter().max_by_key(|e| e.alt_end_i).unwrap();
                    Some(last.alt_end_i) == cds_end
                };
                five_prime_trunc && three_prime_trunc
            } else {
                true
            }
        };

        // Check whether the transcript is predicted.
        let predicted = |p: &Params| -> bool { p.tx_id.starts_with('X') };

        // Check whether the transcript has a biotype that should be discarded.
        let biotype = |p: &Params| -> bool {
            p.tx.biotype
                .as_ref()
                .map(|bt| {
                    bt.iter()
                        .any(|bt| DISCARD_BIOTYPES_TRANSCRIPTS.contains(bt))
                })
                .unwrap_or(false)
        };

        // We have two sets of filters here:
        // `tx_filters`: for every transcript known to us at this point (so every transcript from cdot)
        // `hgnc_grouped_tx_filters`: all transcripts associated to an hgnc id
        // The latter set relies on the grouping of transcripts by hgnc id,
        // e.g. to discard transcripts with an older version within the same hgnc group
        type Filter = fn(&Params) -> bool;
        let tx_filters: [(Filter, Reason); 6] = [
            (missing_hgnc, Reason::MissingHgncId),
            (empty_genome_builds, Reason::EmptyGenomeBuilds),
            (partial, Reason::OnlyPartialAlignmentInRefSeq),
            (predicted, Reason::PredictedTranscript),
            (biotype, Reason::Biotype),
            (cds_start_end_nf, Reason::CdsStartOrEndNotConfirmed),
        ];

        #[derive(new)]
        struct GroupParams<'a> {
            idx: usize,
            tx_ids: &'a [TranscriptId],
            txs: &'a [&'a Transcript],
        }
        type GroupFilter = fn(&GroupParams) -> bool;

        // Check whether we have a newer version of the same transcript accession available.
        let old_version = |p: &GroupParams| -> bool {
            let versioned = Self::by_release_and_version(p.txs);
            let (tx_ac, tx_version) = p.tx_ids[p.idx].split_version();
            for release in p.txs[p.idx].genome_builds.keys() {
                if let Some(other) = versioned.get(&(release.into(), tx_ac.to_string()))
                    && other.iter().any(|(version, _tx)| *version > tx_version)
                {
                    return true;
                }
            }
            false
        };

        // Check whether we have already seen an NM transcript for the gene.
        let nm_instead_of_nr = |p: &GroupParams| -> bool {
            p.tx_ids[p.idx].starts_with("NR_")
                && p.tx_ids.iter().any(|tx_id| tx_id.starts_with("NM_"))
        };

        let hgnc_grouped_tx_filters: [(GroupFilter, Reason); 2] = [
            (old_version, Reason::OldVersion),
            (nm_instead_of_nr, Reason::UseNmTranscriptInsteadOfNr),
        ];

        let start = Instant::now();

        // Apply first set of filters (which do not depend on hgnc grouping)
        let discarded = self
            .transcript_id_to_transcript
            .iter()
            .filter_map(|(tx_id, tx)| {
                let p = Params::new(tx, tx_id);
                let reason = tx_filters
                    .iter()
                    .filter_map(|(f, r)| f(&p).then_some(*r))
                    .fold(BitFlags::<Reason>::default(), |a, b| a | b);
                let empty = reason.is_empty();
                (!empty).then(|| (Identifier::Transcript(tx_id.clone()), reason))
            })
            .collect_vec();

        for (id, reason) in discarded {
            self.mark_discarded(&id, reason)?;
        }

        // Apply second set of filters (which depend on hgnc grouping)
        let discarded = self
            .gene_id_to_transcript_ids
            .iter()
            .filter_map(|(_gene_id, tx_ids)| {
                let txs = tx_ids
                    .iter()
                    .map(|tx_id| self.transcript_id_to_transcript.get(tx_id).unwrap())
                    .collect_vec();
                let tx_reasons = tx_ids
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, tx)| {
                        let p = GroupParams::new(idx, tx_ids, &txs);
                        let reason = hgnc_grouped_tx_filters
                            .iter()
                            .filter_map(|(f, r)| f(&p).then_some(*r))
                            .fold(BitFlags::<Reason>::default(), |a, b| a | b);
                        let empty = reason.is_empty();
                        if !empty {
                            Some((Identifier::Transcript(tx.clone()), reason))
                        } else {
                            None
                        }
                    })
                    .collect_vec();
                if tx_reasons.is_empty() {
                    None
                } else {
                    Some(tx_reasons)
                }
            })
            .flatten()
            .collect_vec();

        for (id, reason) in discarded {
            self.mark_discarded(&id, reason)?;
        }
        // self.discard()?;

        tracing::info!("… done filtering transcripts in {:?}", start.elapsed());
        Ok(())
    }

    fn filter_transcripts_with_sequence(
        &mut self,
        seq_provider: &mut SequenceProvider,
    ) -> Result<HashMap<TranscriptId, String>, Error> {
        tracing::info!("Filtering transcripts with sequences …");
        let start = Instant::now();

        let five_prime_truncated = |tx: &Transcript| -> bool {
            let is_mt = MITOCHONDRIAL_ACCESSIONS.iter().any(|a| tx.is_on_contig(a));
            if tx.protein_coding() && !is_mt {
                tx.genome_builds.iter().any(|(_release, alignment)| {
                    let cds_start = alignment.cds_start;
                    let mut exons = alignment.exons.clone();
                    exons.is_empty() || {
                        exons.sort_unstable_by_key(|e| e.alt_start_i);
                        let first = exons.first().unwrap();
                        Some(first.alt_start_i) == cds_start
                    }
                })
            } else {
                false
            }
        };
        let three_prime_truncated = |tx: &Transcript| -> bool {
            let is_mt = MITOCHONDRIAL_ACCESSIONS.iter().any(|a| tx.is_on_contig(a));
            if tx.protein_coding() && !is_mt {
                tx.genome_builds.iter().any(|(_release, alignment)| {
                    let cds_end = alignment.cds_end;
                    let mut exons = alignment.exons.clone();
                    exons.is_empty() || {
                        exons.sort_unstable_by_key(|e| e.alt_end_i);
                        let last = exons.last().unwrap();
                        Some(last.alt_end_i) == cds_end
                    }
                })
            } else {
                false
            }
        };

        // Check transcript's CDS length for being a multiple of 3.
        //
        // Note that the chrMT transcripts have been fixed earlier already to
        // accommodate for how they are fixed by poly-A tailing.
        let invalid_cds_length = |tx: &Transcript| -> bool {
            if tx.protein_coding() {
                tx.cds_length().is_none_or(|l| l % 3 != 0)
            } else {
                false
            }
        };

        let (discards, keeps_with_reasons): (Vec<_>, Vec<_>) = self
            .transcript_id_to_transcript
            .par_iter()
            .partition_map(|(tx_id, tx)| {
                if let Some(d) = self.discards.get(&Identifier::Transcript(tx_id.clone()))
                    && d.intersects(Reason::hard())
                {
                    return Either::Left((Identifier::Transcript(tx_id.clone()), *d));
                }

                let has_invalid_cds_length = if self.disable_filters {
                    false
                } else {
                    invalid_cds_length(tx)
                };

                let namespace: String = if tx_id.starts_with("ENST") {
                    String::from("Ensembl")
                } else {
                    String::from("NCBI")
                };
                let seq = seq_provider.fetch_sequence(&AliasOrSeqId::Alias {
                    value: (*tx_id).to_string(),
                    namespace: Some(namespace.clone()),
                });

                if let Ok(seq) = seq {
                    let is_mt = MITOCHONDRIAL_ACCESSIONS.iter().any(|a| tx.is_on_contig(a));
                    let is_seleno = tx
                        .biotype
                        .as_ref()
                        .map(|bt| bt.contains(&BioType::Selenoprotein))
                        .unwrap_or(false);
                    if seq.is_empty() {
                        return Either::Left((
                            Identifier::Transcript(tx_id.clone()),
                            Reason::MissingSequence.into(),
                        ));
                    };
                    // Skip transcript if it is coding and the translated CDS does not have a stop codon.
                    let cds = tx.start_codon.and_then(|start| {
                        tx.stop_codon.map(|stop| (start as usize, stop as usize))
                    });
                    let cds = if tx.protein_coding() { cds } else { None };
                    if let Some((cds_start, cds_end)) = cds {
                        let cds_length = cds_end.saturating_sub(cds_start);
                        let delta = (3 - (cds_length % 3)).max(cds_end.saturating_sub(seq.len()));
                        let seq = append_poly_a(seq, delta);

                        let safe_end = cds_end.min(seq.len());
                        let safe_start = cds_start.min(safe_end);

                        let tx_seq_to_translate = &seq[safe_start..safe_end];

                        let aa_sequence_result = translate_cds(
                            tx_seq_to_translate,
                            true,
                            "*",
                            if is_mt {
                                TranslationTable::VertebrateMitochondrial
                            } else if is_seleno {
                                TranslationTable::Selenocysteine
                            } else {
                                TranslationTable::Standard
                            },
                        );

                        let mut reason = Reason::empty();

                        if let Ok(aa_sequence) = aa_sequence_result {
                            let has_missing_stop_codon = (!is_mt && !aa_sequence.ends_with('*'))
                                || (is_mt && !aa_sequence.contains('*'));
                            if has_missing_stop_codon && !self.disable_filters {
                                reason |= Reason::MissingStopCodon;
                                if five_prime_truncated(tx) {
                                    reason |= Reason::FivePrimeEndTruncated;
                                }
                                if three_prime_truncated(tx) {
                                    reason |= Reason::ThreePrimeEndTruncated;
                                }
                            }
                        }

                        if has_invalid_cds_length
                            || self
                                .fixes
                                .get(&Identifier::Transcript(tx_id.clone()))
                                .is_some_and(|f| f.contains(Fix::Cds))
                        {
                            reason |= Reason::InvalidCdsLength;
                        }

                        if reason.intersects(Reason::hard()) {
                            Either::Left((Identifier::Transcript(tx_id.clone()), reason))
                        } else {
                            Either::Right((tx_id.clone(), seq, reason))
                        }
                    } else {
                        Either::Right((tx_id.clone(), seq, Reason::empty()))
                    }
                } else {
                    Either::Left((
                        Identifier::Transcript(tx_id.clone()),
                        Reason::MissingSequence.into(),
                    ))
                }
            });

        for (id, reason) in discards {
            self.mark_discarded(&id, reason)?;
        }
        // self.discard()?;

        let mut sequence_map = HashMap::new();
        for (tx_id, seq, reason) in keeps_with_reasons {
            if !reason.is_empty() {
                self.mark_discarded(&Identifier::Transcript(tx_id.clone()), reason)?;
            }
            sequence_map.insert(tx_id, seq);
        }

        tracing::info!(
            "… done filtering transcripts with sequence provider in {:?}",
            start.elapsed()
        );

        Ok(sequence_map)
    }

    /// Remove all hgnc ids for which there are no transcripts left.
    fn filter_empty_hgnc_mappings(&mut self) -> Result<(), Error> {
        tracing::info!("Removing empty HGNC mappings …");
        let gene_ids_gene: HashSet<_> = self.gene_id_to_gene.keys().collect();
        let gene_ids_tx: HashSet<_> = self.gene_id_to_transcript_ids.keys().collect();
        let gene_id_but_no_tx_id = &gene_ids_gene - &gene_ids_tx;
        let tx_id_but_no_gene_id = &gene_ids_tx - &gene_ids_gene;
        for gene_id in gene_id_but_no_tx_id {
            *self
                .discards
                .entry(Identifier::Gene(gene_id.clone()))
                .or_default() |= Reason::NoTranscripts;
        }
        for gene_id in tx_id_but_no_gene_id {
            *self
                .discards
                .entry(Identifier::Gene(gene_id.clone()))
                .or_default() |= Reason::MissingGene;
        }

        for (gene_id, txs) in self.gene_id_to_transcript_ids.iter() {
            if txs.is_empty() {
                *self
                    .discards
                    .entry(Identifier::Gene(gene_id.clone()))
                    .or_default() |= Reason::NoTranscriptLeft;
            }
        }
        // self.discard()?;
        tracing::info!("… done removing empty HGNC mappings");
        Ok(())
    }

    fn update_pseudogene_status(&mut self) -> Result<(), Error> {
        let empty = Reason::empty();
        for (gene_id, gene) in &self.gene_id_to_gene {
            if !self
                .discards
                .get(&Identifier::Gene(gene_id.clone()))
                .unwrap_or(&empty)
                .contains(Reason::NoTranscriptLeft)
            {
                continue;
            }
            if gene
                .biotype
                .as_ref()
                .is_some_and(|bt| bt.contains(&BioType::Pseudogene))
            {
                *self
                    .discards
                    .entry(Identifier::Gene(gene_id.clone()))
                    .or_default() |= Reason::Pseudogene;
            }
        }
        Ok(())
    }

    fn gather_transcript_stats(&self) -> Result<(usize, usize, usize), Error> {
        let (mut n_mane_select, mut n_mane_plus_clinical, mut n_mt) = (0, 0, 0);
        for tx in self.transcript_id_to_transcript.values() {
            let mut is_mane_select = false;
            let mut is_mane_plus_clinical = false;
            for gb in tx.genome_builds.values() {
                if MITOCHONDRIAL_ACCESSIONS.contains(&gb.contig.as_str()) {
                    n_mt += 1;
                }
                if let Some(tag) = &gb.tag {
                    if tag.contains(&Tag::ManeSelect) {
                        is_mane_select = true;
                    }
                    if tag.contains(&Tag::ManePlusClinical) {
                        is_mane_plus_clinical = true;
                    }
                }
            }
            if is_mane_select {
                n_mane_select += 1;
            }
            if is_mane_plus_clinical {
                n_mane_plus_clinical += 1;
            }
        }
        Ok((n_mt, n_mane_select, n_mane_plus_clinical))
    }

    fn update_transcript_tags(
        &mut self,
        transcript_id_to_tags: &HashMap<TranscriptId, Vec<FeatureTag>>,
    ) {
        self.transcript_id_to_transcript
            .iter_mut()
            .for_each(|(tx_id, tx)| {
                // transfer MANE-related labels from TSV file
                let tx_id_no_version = tx_id.without_version().unwrap();
                if let Some(tags) = transcript_id_to_tags.get(&tx_id_no_version) {
                    let mut any_change = false;
                    tx.genome_builds.iter_mut().for_each(|(_, alignment)| {
                        if let Some(alignment_tag) = &mut alignment.tag {
                            let tags_to_add: Vec<Tag> = tags
                                .iter()
                                .filter_map(|t| {
                                    let standard = Tag::from(t.clone());
                                    if alignment_tag.contains(&standard) {
                                        None
                                    } else {
                                        let backported = Tag::from(t.to_backported());
                                        if alignment_tag.contains(&backported) {
                                            None
                                        } else {
                                            Some(backported)
                                        }
                                    }
                                })
                                .collect();

                            if !tags_to_add.is_empty() {
                                alignment_tag.extend(tags_to_add);
                                alignment_tag.sort();
                                alignment_tag.dedup();
                                any_change = true;
                            }
                        }
                    });

                    if any_change {
                        self.fixes
                            .entry(Identifier::Transcript(tx_id.clone()))
                            .or_default()
                            .insert(Fix::Tags);
                    }
                }
            });
    }

    fn fix_transcript_genome_builds(&mut self) {
        let genome_release = self.genome_release.clone();
        self.transcript_id_to_transcript
            .values_mut()
            .for_each(|tx| {
                let n = tx.genome_builds.len();
                tx.genome_builds
                    .retain(|key, _| key.eq_ignore_ascii_case(&genome_release));
                if n != tx.genome_builds.len() {
                    self.fixes
                        .entry(Identifier::Transcript(
                            TranscriptId::try_new(&tx.id).unwrap(),
                        ))
                        .or_default()
                        .insert(Fix::GenomeBuild);
                }
            });
    }

    fn fix_cds(&mut self) {
        self.transcript_id_to_transcript
            .values_mut()
            .filter(|tx| tx.protein_coding())
            .filter_map(|tx| {
                tx.start_codon
                    .and_then(|start| tx.stop_codon.map(|stop| (start, stop)))
                    .and_then(|(cds_start, cds_end)| {
                        let cds_len = cds_end - cds_start;
                        (cds_len % 3 != 0).then_some((tx, cds_start, cds_end, cds_len))
                    })
            })
            .for_each(|(tx, _cds_start, cds_end, cds_len)| {
                assert_eq!(tx.genome_builds.len(), 1);
                for gb in tx.genome_builds.values_mut() {
                    let delta = 3 - (cds_len % 3);
                    if delta == 0 {
                        continue;
                    };
                    tx.stop_codon = Some(cds_end + delta);
                    let exon = gb
                        .exons
                        .iter_mut()
                        .max_by_key(|g| g.alt_cds_end_i)
                        .expect("No exons found during fix_cds");
                    exon.alt_cds_end_i += delta;
                    exon.cigar.push_str(&format!("{}I", delta));
                    self.fixes
                        .entry(Identifier::Transcript(
                            TranscriptId::try_new(&tx.id).unwrap(),
                        ))
                        .or_default()
                        .insert(Fix::Cds);
                }
            });
    }

    fn gene_name(&self, id: &Identifier) -> Option<String> {
        match id {
            Identifier::Gene(gene_id) => self
                .gene_id_to_gene
                .get(gene_id)
                .and_then(|gene| gene.gene_symbol.clone()),
            Identifier::Transcript(tx_id) => self
                .transcript_id_to_transcript
                .get(tx_id)
                .and_then(|tx| tx.gene_name.clone()),
        }
    }

    fn tags(&self, id: &Identifier) -> Option<Vec<Tag>> {
        match id {
            Identifier::Transcript(tx_id) => {
                self.transcript_id_to_transcript.get(tx_id).map(|tx| {
                    tx.genome_builds
                        .values()
                        .flat_map(|gb| gb.tag.clone())
                        .flatten()
                        .sorted_unstable()
                        .dedup()
                        .collect()
                })
            }
            Identifier::Gene(gene_id) => {
                self.gene_id_to_transcript_ids.get(gene_id).map(|tx_ids| {
                    tx_ids
                        .iter()
                        .flat_map(|tx_id| {
                            self.tags(&Identifier::Transcript(tx_id.clone()))
                                .unwrap_or_default()
                        })
                        .sorted_unstable()
                        .dedup()
                        .collect()
                })
            }
        }
    }

    /// Group transcripts by release and version (sorted descending).
    fn by_release_and_version<'a>(
        txs: &[&'a Transcript],
    ) -> HashMap<(String, String), Vec<(u32, &'a Transcript)>> {
        let mut versioned = txs
            .iter()
            .flat_map(|tx| {
                let releases = tx.genome_builds.keys().cloned().collect_vec();
                let tx_id = TranscriptId::try_new(&tx.id).unwrap();
                let (ac, version) = tx_id.split_version();
                let ac = ac.to_string();
                releases
                    .into_iter()
                    .map(move |r| ((r, ac.clone()), (version, *tx)))
            })
            .into_group_map();
        versioned.values_mut().for_each(|txs| {
            txs.sort_unstable_by_key(|(version, _)| Reverse(*version));
        });
        versioned
    }

    fn discard(&mut self, remove: bool) -> Result<(), Error> {
        let (n_transcripts_pre, n_hgnc_ids_pre) = (
            self.transcript_id_to_transcript.len(),
            self.gene_id_to_transcript_ids.len(),
        );
        type Counter = HashMap<Reason, usize>;
        #[derive(Debug, Clone, Copy, Display, PartialEq, Eq, Hash)]
        enum Kind {
            Transcript,
            GeneId,
        }
        let mut count_groups = HashMap::<Kind, Counter>::new();
        let kinds = [Kind::Transcript, Kind::GeneId];
        for kind in &kinds {
            let mut counter = Counter::new();
            for r in Reason::all().iter() {
                counter.entry(r).or_default();
            }
            count_groups.insert(*kind, counter);
        }

        for (id, reason) in self.discards.clone() {
            if reason.is_empty() {
                tracing::warn!("Empty discard reason for {:#?}", id);
                continue;
            }
            self._discard_id(&id, reason, remove)?;
        }

        for (id, reason) in &self.discards {
            match id {
                Identifier::Transcript(_) => {
                    for r in reason.iter() {
                        *count_groups
                            .entry(Kind::Transcript)
                            .or_default()
                            .entry(r)
                            .or_default() += 1;
                    }
                }
                Identifier::Gene(_) => {
                    for r in reason.iter() {
                        *count_groups
                            .entry(Kind::GeneId)
                            .or_default()
                            .entry(r)
                            .or_default() += 1;
                    }
                }
            }
        }
        let (n_transcripts_post, n_hgnc_ids_post) = (
            self.transcript_id_to_transcript.len(),
            self.gene_id_to_transcript_ids.len(),
        );
        tracing::info!(
            "Discarded {} transcripts and {} HGNC IDs",
            n_transcripts_pre.abs_diff(n_transcripts_post),
            n_hgnc_ids_pre.abs_diff(n_hgnc_ids_post)
        );
        for kind in kinds {
            tracing::info!(
                "{}: {:?}",
                kind,
                Reason::all()
                    .iter()
                    .filter(|r| count_groups[&kind][r] > 0)
                    .map(|r| format!("{:#?}: {}", r, count_groups[&kind][&r]))
                    .join(", ")
            );
        }
        Ok(())
    }

    fn _discard_id(
        &mut self,
        id: &Identifier,
        reason: BitFlags<Reason>,
        remove: bool,
    ) -> Result<(), Error> {
        *self.discards.entry(id.clone()).or_default() |= reason;
        match id {
            Identifier::Gene(gene_id) => {
                let txs = if remove {
                    self.gene_id_to_transcript_ids.remove(gene_id)
                } else {
                    self.gene_id_to_transcript_ids.get(gene_id).cloned()
                };
                if remove {
                    self.gene_id_to_gene.remove(gene_id);
                }
                for tx_id in txs.unwrap_or_default() {
                    *self
                        .discards
                        .entry(Identifier::Transcript(tx_id.clone()))
                        .or_default() |= reason;
                    if remove {
                        let _transcript = self.transcript_id_to_transcript.remove(&tx_id);
                    }
                }
            }
            Identifier::Transcript(tx_id) => {
                let transcript = if remove {
                    self.transcript_id_to_transcript.remove(tx_id)
                } else {
                    self.transcript_id_to_transcript.get(tx_id).cloned()
                };

                if let Some(hgnc) = transcript.and_then(|t| t.hgnc) {
                    let gene_id: GeneId = hgnc.parse().expect("Invalid GeneId");

                    if remove {
                        let txs = self.gene_id_to_transcript_ids.get_mut(&gene_id);
                        if let Some(txs) = txs {
                            txs.retain(|tx| tx != tx_id);
                            if txs.is_empty() {
                                *self
                                    .discards
                                    .entry(Identifier::Gene(gene_id.clone()))
                                    .or_default() |= Reason::NoTranscriptLeft;
                                self.gene_id_to_transcript_ids.remove(&gene_id);
                                self.gene_id_to_gene.remove(&gene_id);
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn mark_discarded(&mut self, id: &Identifier, reason: BitFlags<Reason>) -> Result<(), Error> {
        if reason.is_empty() {
            panic!("Empty discard reason for {:#?}", id);
        }
        *self.discards.entry(id.clone()).or_default() |= reason;
        Ok(())
    }

    /// For each transcript that has been discarded for whatever reason, propagate the reason to its
    /// parent hgnc entry *if* the hgnc entry already exists and has a non-empty reason
    /// (to avoid erroneously discarding hgnc entries for a non-important reason).
    fn propagate_discard_reasons(&mut self, _raw: &Self) -> Result<(), Error> {
        // First check whether all transcripts of a gene have been marked as discarded.
        for (gene_id, _) in self.gene_id_to_gene.iter() {
            let tx_ids = self
                .gene_id_to_transcript_ids
                .get(gene_id)
                .map(|v| v.as_slice())
                .unwrap_or_default();
            if !tx_ids.is_empty()
                && tx_ids.iter().all(|tx_id| {
                    self.discards
                        .get(&Identifier::Transcript(tx_id.clone()))
                        .is_some_and(|d| d.intersects(Reason::hard()))
                })
            {
                *self
                    .discards
                    .entry(Identifier::Gene(gene_id.clone()))
                    .or_default() |= Reason::NoTranscriptLeft;
            }
        }

        Ok(())
    }

    /// Perform protobuf file construction.
    ///
    /// This can be done by simply converting the models from ``hgvs-rs`` to the prost generated data structures.
    fn build_protobuf(
        &mut self,
        sequence_map: &mut HashMap<TranscriptId, String>,
        source_version: SourceVersion,
    ) -> Result<TxSeqDatabase, Error> {
        tracing::info!("Constructing protobuf data structures …");
        let start = Instant::now();

        let (gene_ids, tx_ids): (Vec<GeneId>, Vec<TranscriptId>) = {
            // ensure we use all hgnc_ids that are available to us;
            // this includes hgnc_ids for which we only have transcripts but no genes or vice versa
            let gene_ids = self
                .gene_id_to_gene
                .keys()
                .cloned()
                .chain(self.gene_id_to_transcript_ids.keys().cloned())
                .sorted_unstable()
                .dedup()
                .collect_vec();

            // We do, however, discard transcripts that have no sequence or hgnc id.
            let tx_ids = self
                .gene_id_to_transcript_ids
                .values()
                .flatten()
                .sorted_unstable()
                .dedup()
                .filter(|&tx_id| {
                    self.discards
                        .get(&Identifier::Transcript(tx_id.clone()))
                        .is_none_or(|reason| !reason.intersects(Reason::hard()))
                })
                .cloned()
                .collect();
            (gene_ids, tx_ids)
        };

        // Construct sequence database.
        tracing::info!("  Constructing sequence database …");
        let seq_db = {
            // filter out transcripts that have no sequence, enumerate, collect; build SequenceDb
            let (aliases, aliases_idx, seqs) = tx_ids
                .iter()
                .filter_map(|tx_id| {
                    sequence_map
                        .remove(tx_id)
                        .map(|seq| ((*tx_id).to_string(), seq))
                })
                .enumerate()
                .map(|(idx, (tx_id, seq))| (tx_id, idx as u32, seq))
                .multiunzip();

            crate::pbs::txs::SequenceDb {
                aliases,
                aliases_idx,
                seqs,
            }
        };

        tracing::info!("  Creating transcript records for each gene…");
        let empty_txs = vec![];
        let data_transcripts = gene_ids
            .iter()
            .flat_map(|gene_id| {
                let tx_ids = self
                    .gene_id_to_transcript_ids
                    .get(gene_id)
                    .unwrap_or(&empty_txs);
                tx_ids
                    .iter()
                    .sorted_unstable()
                    .map(|tx_id| self.protobuf_transcript(gene_id, tx_id))
            })
            .collect();

        tracing::info!(" … done creating transcripts in {:#?}", start.elapsed());

        // Build mapping of gene HGNC symbol to transcript IDs.
        tracing::info!("  Build gene symbol to transcript ID mapping …");
        let empty = vec![];
        let gene_to_tx = gene_ids
            .iter()
            .map(|gene_id| {
                let tx_ids = self
                    .gene_id_to_transcript_ids
                    .get(gene_id)
                    .unwrap_or(&empty);
                let hgnc_reason = self.discards.get(&Identifier::Gene(gene_id.clone()));
                let filtered = hgnc_reason.is_some_and(|reason| reason.intersects(Reason::hard()));
                crate::pbs::txs::GeneToTxId {
                    gene_id: gene_id.to_string(),
                    tx_ids: tx_ids
                        .iter()
                        .sorted_unstable()
                        .map(|tx_id| tx_id.to_string())
                        .collect(),
                    filtered: Some(filtered),
                    filter_reason: hgnc_reason.map(|r| r.bits()),
                }
            })
            .collect::<Vec<_>>();
        tracing::info!(" … done building gene symbol to transcript ID mapping");

        // Compose transcript database from transcripts and gene to transcript mapping.
        tracing::info!("  Composing transcript seq database …");
        let tx_db = crate::pbs::txs::TranscriptDb {
            transcripts: data_transcripts,
            gene_to_tx,
        };

        // Compose the final transcript and sequence database.
        let tx_seq_db = TxSeqDatabase {
            tx_db: Some(tx_db),
            seq_db: Some(seq_db),
            version: Some(crate::common::version().to_string()),
            source_version: vec![source_version],
        };

        tracing::info!(" … done composing transcript seq database");

        Ok(tx_seq_db)
    }

    fn protobuf_transcript(
        &self,
        gene_id: &GeneId,
        tx_id: &TranscriptId,
    ) -> crate::pbs::txs::Transcript {
        let tx = self
            .transcript_id_to_transcript
            .get(tx_id)
            .unwrap_or_else(|| panic!("No transcript for id {:?}", tx_id));

        // Combine discard reasons for transcript and gene level.
        let tx_reason = self
            .discards
            .get(&Identifier::Transcript(tx_id.clone()))
            .copied()
            .unwrap_or_default();
        let gene_reason = self
            .discards
            .get(&Identifier::Gene(gene_id.clone()))
            .copied()
            .unwrap_or_default();
        let combined_reason = tx_reason | gene_reason;

        // Only filter out the transcript if it intersects with a "hard" reason.
        // Soft-filtered transcripts will have `filtered = false`, but their `filter_reason`
        // will still contain the bitmask for downstream inspection.
        let filtered = combined_reason.intersects(Reason::hard());

        // ... build genome alignment for selected:
        let (genome_alignments, tags): (Vec<_>, Vec<_>) = tx
            .genome_builds
            .iter()
            .map(|(genome_build, alignment)| {
                Self::protobuf_genome_alignment(genome_build, alignment)
            })
            .unzip();

        let seleno_tag = tx
            .biotype
            .as_ref()
            .is_some_and(|bt| bt.contains(&BioType::Selenoprotein))
            .then(|| crate::pbs::txs::TranscriptTag::Selenoprotein.into());

        let tags = tags
            .into_iter()
            .flatten()
            .chain(seleno_tag)
            .sorted_unstable()
            .dedup()
            .collect();

        // Now, just obtain the basic properties and create a new `data::Transcript`.
        let Gene { gene_symbol, .. } =
            self.gene_id_to_gene
                .get(gene_id)
                .cloned()
                .unwrap_or_else(|| Gene {
                    aliases: None,
                    biotype: None,
                    description: None,
                    gene_symbol: None,
                    hgnc: Some(gene_id.to_string()),
                    map_location: None,
                    summary: None,
                    url: "".to_string(),
                });

        let Transcript {
            protein,
            start_codon,
            stop_codon,
            ..
        } = tx.clone();
        let biotype = if tx.protein_coding() {
            crate::pbs::txs::TranscriptBiotype::Coding.into()
        } else {
            crate::pbs::txs::TranscriptBiotype::NonCoding.into()
        };

        crate::pbs::txs::Transcript {
            id: (*tx_id).to_string(),
            gene_symbol: gene_symbol.unwrap_or("<MISSING>".to_string()),
            gene_id: gene_id.to_string(),
            biotype,
            tags,
            protein,
            start_codon,
            stop_codon,
            genome_alignments,
            filtered: Some(filtered),
            filter_reason: combined_reason
                .is_empty()
                .not()
                .then(|| combined_reason.bits()),
        }
    }

    fn protobuf_genome_alignment(
        genome_build: &str,
        alignment: &GenomeAlignment,
    ) -> (crate::pbs::txs::GenomeAlignment, HashSet<i32>) {
        let mut tags = HashSet::new();
        let genome_build = genome_build.to_lowercase();

        let GenomeAlignment {
            contig,
            cds_start,
            cds_end,
            ..
        } = alignment.clone();
        let strand = match alignment.strand {
            models::Strand::Plus => crate::pbs::txs::Strand::Plus,
            models::Strand::Minus => crate::pbs::txs::Strand::Minus,
        };
        if let Some(tag) = alignment.tag.as_ref() {
            for t in tag {
                let elem = match t {
                    Tag::Basic => crate::pbs::txs::TranscriptTag::Basic.into(),
                    Tag::EnsemblCanonical => {
                        crate::pbs::txs::TranscriptTag::EnsemblCanonical.into()
                    }
                    Tag::ManeSelect => crate::pbs::txs::TranscriptTag::ManeSelect.into(),
                    Tag::ManePlusClinical => {
                        crate::pbs::txs::TranscriptTag::ManePlusClinical.into()
                    }
                    Tag::RefSeqSelect => crate::pbs::txs::TranscriptTag::RefSeqSelect.into(),
                    Tag::GencodePrimary => crate::pbs::txs::TranscriptTag::GencodePrimary.into(),
                    Tag::Other(v) => match v.as_str() {
                        "EnsemblGraft" => crate::pbs::txs::TranscriptTag::EnsemblGraft.into(),
                        "basic-backport" => crate::pbs::txs::TranscriptTag::BasicBackport.into(),
                        "ensembl_canonical-backport" => {
                            crate::pbs::txs::TranscriptTag::EnsemblCanonicalBackport.into()
                        }
                        "mane_select-backport" => {
                            crate::pbs::txs::TranscriptTag::ManeSelectBackport.into()
                        }
                        "mane_plus_clinical-backport" => {
                            crate::pbs::txs::TranscriptTag::ManePlusClinicalBackport.into()
                        }
                        "ref_seq_select-backport" => {
                            crate::pbs::txs::TranscriptTag::RefSeqSelectBackport.into()
                        }
                        "selenoprotein-backport" => {
                            crate::pbs::txs::TranscriptTag::SelenoproteinBackport.into()
                        }
                        "gencode_primary-backport" => {
                            crate::pbs::txs::TranscriptTag::GencodePrimaryBackport.into()
                        }
                        x if x.ends_with("-backport") => {
                            crate::pbs::txs::TranscriptTag::OtherBackport.into()
                        }
                        _ => crate::pbs::txs::TranscriptTag::Other.into(),
                    },
                };
                tags.insert(elem);
            }
        }
        // Look into any "note" string for a selenoprotein marker and
        // add this as a tag.
        if let Some(note) = alignment.note.as_ref() {
            let needle = "UGA stop codon recoded as selenocysteine";
            if note.contains(needle) {
                tags.insert(crate::pbs::txs::TranscriptTag::Selenoprotein.into());
            }
        }
        // and construct vector of all exons
        let exons: Vec<_> = alignment
            .exons
            .iter()
            .map(|exon| {
                let models::Exon {
                    alt_start_i,
                    alt_end_i,
                    ord,
                    alt_cds_start_i,
                    alt_cds_end_i,
                    cigar,
                } = exon.clone();
                crate::pbs::txs::ExonAlignment {
                    alt_start_i,
                    alt_end_i,
                    ord,
                    alt_cds_start_i: if alt_cds_start_i == -1 {
                        None
                    } else {
                        Some(alt_cds_start_i)
                    },
                    alt_cds_end_i: if alt_cds_end_i == -1 {
                        None
                    } else {
                        Some(alt_cds_end_i)
                    },
                    cigar,
                }
            })
            .collect();

        let genome_alignment = crate::pbs::txs::GenomeAlignment {
            genome_build,
            contig,
            cds_start,
            cds_end,
            strand: strand.into(),
            exons,
        };
        (genome_alignment, tags)
    }
}

fn append_poly_a(seq: String, length: usize) -> String {
    // Append poly-A for chrMT transcripts (which are from ENSEMBL).
    // This also potentially fixes the stop codon.
    let mut seq = seq.into_bytes();
    seq.extend_from_slice(b"A".repeat(length).as_slice());
    String::from_utf8(seq).expect("must be valid UTF-8")
}

#[bitflags]
#[repr(u32)]
#[derive(Debug, Clone, Copy, Serialize, Hash, PartialEq, Eq)]
pub(crate) enum Reason {
    Biotype,
    CdsEndAfterSequenceEnd,
    DeselectedGene,
    EmptyGenomeBuilds,
    FivePrimeEndTruncated,
    InvalidCdsLength,
    MissingGene,
    MissingGeneSymbol,
    MissingHgncId,
    MissingSequence,
    MissingStopCodon,
    NoTranscriptLeft,
    NoTranscripts,
    OldVersion,
    OnlyPartialAlignmentInRefSeq,
    PredictedTranscript,
    PredictedTranscriptsOnly,
    Pseudogene,
    ThreePrimeEndTruncated,
    TranscriptPriority,
    UseNmTranscriptInsteadOfNr,
    CdsStartOrEndNotConfirmed,
}

impl Reason {
    /// Reasons that make a transcript completely unusable (Hard Filter).
    pub fn hard() -> BitFlags<Reason> {
        Reason::MissingSequence
            | Reason::EmptyGenomeBuilds
            | Reason::DeselectedGene
            | Reason::NoTranscripts
            | Reason::NoTranscriptLeft
    }

    #[allow(dead_code)]
    /// Reasons that mean the transcript is flawed but salvageable (Soft Filter).
    pub fn soft() -> BitFlags<Reason> {
        BitFlags::all() ^ Self::hard()
    }
}

#[bitflags]
#[repr(u8)]
#[derive(Debug, Clone, Copy, Serialize, Hash, PartialEq, Eq)]
enum Fix {
    Cds,
    GenomeBuild,
    Tags,
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "type", content = "value")]
enum ReportEntry {
    Discard(Discard),
    SoftFilter(SoftFilter),
    Fix(LogFix),
    Log(serde_json::Value),
}
#[serde_as]
#[derive(Debug, Clone, Serialize)]
struct Discard {
    source: String,
    #[serde_as(as = "DisplayFromStr")]
    reason: BitFlags<Reason>,
    id: Identifier,
    gene_name: Option<String>,
    tags: Option<Vec<Tag>>,
}

#[serde_as]
#[derive(Debug, Clone, Serialize)]
struct SoftFilter {
    source: String,
    #[serde_as(as = "DisplayFromStr")]
    reason: BitFlags<Reason>,
    id: Identifier,
    gene_name: Option<String>,
    tags: Option<Vec<Tag>>,
}

#[serde_as]
#[derive(Debug, Clone, Serialize)]
struct LogFix {
    source: String,
    #[serde_as(as = "DisplayFromStr")]
    fix: BitFlags<Fix>,
    id: Identifier,
    gene_name: Option<String>,
    tags: Option<Vec<Tag>>,
}

fn txid_to_label(
    label_tsv_path: impl AsRef<Path>,
) -> Result<HashMap<TranscriptId, Vec<FeatureTag>>, Error> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_path(label_tsv_path.as_ref())?;

    rdr.deserialize()
        .map(|result| {
            result
                .map_err(anyhow::Error::from)
                .and_then(|entry: LabelEntry| {
                    TranscriptId::try_new(entry.transcript_id)
                        .map(|txid| {
                            (
                                txid,
                                entry
                                    .label
                                    .split(',')
                                    .map(|s| {
                                        let cdot_tag = models::str_to_tag(s);
                                        FeatureTag::from(cdot_tag)
                                    })
                                    .collect::<Vec<_>>(),
                            )
                        })
                        .map_err(anyhow::Error::from)
                })
        })
        .collect()
}
pub(crate) fn read_cdot_json(path: impl AsRef<Path>) -> Result<models::Container, Error> {
    Ok(if path.as_ref().extension().unwrap_or_default() == "gz" {
        tracing::info!("(from gzip compressed file)");
        serde_json::from_reader(std::io::BufReader::new(flate2::read::GzDecoder::new(
            File::open(path)?,
        )))?
    } else {
        tracing::info!("(from uncompressed file)");
        serde_json::from_reader(std::io::BufReader::new(File::open(path)?))?
    })
}

/// Create file-backed `SeqRepo`.
fn open_seqrepo(path: impl AsRef<Path>) -> Result<SeqRepo, Error> {
    tracing::info!("Opening seqrepo…");
    let start = Instant::now();
    let seqrepo = PathBuf::from(path.as_ref());
    let p = path.as_ref().to_str();
    let path = seqrepo
        .parent()
        .ok_or(anyhow::anyhow!("Could not get parent from {:?}", &p))?
        .to_str()
        .unwrap()
        .to_string();
    let instance = seqrepo
        .file_name()
        .ok_or(anyhow::anyhow!("Could not get basename from {:?}", &p))?
        .to_str()
        .unwrap()
        .to_string();
    let seqrepo = SeqRepo::new(path, &instance)?;
    tracing::info!("… seqrepo opened in {:?}", start.elapsed());
    Ok(seqrepo)
}

/// Load the annotations (JSON or GFF3).
fn load_annotations(args: &Args) -> Result<TranscriptLoader, Error> {
    tracing::info!("Loading annotations …");
    let start = Instant::now();
    let labels = args
        .mane_transcripts
        .as_ref()
        .map(txid_to_label)
        .transpose()?;
    let loaders = args
        .annotation
        .iter()
        .map(|path| {
            let mut loader = TranscriptLoader::new(args.assembly.clone(), args.disable_filters);

            let ext = path.extension().unwrap_or_default().to_string_lossy();
            let file_stem = path.file_stem().unwrap_or_default().to_string_lossy();
            let is_gff3 = ext.ends_with("gff")
                || ext.ends_with("gff3")
                || (ext == "gz" && (file_stem.ends_with("gff3") || file_stem.ends_with("gff")));

            if is_gff3 {
                loader.load_gff3(path).map(|_| loader)
            } else {
                loader.load_cdot(path).map(|_| loader)
            }
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let mut merged = loaders
        .into_iter()
        .reduce(|mut a, mut b| {
            a.merge(&mut b);
            a
        })
        .unwrap();
    merged.apply_fixes(&labels);

    tracing::info!(
        "… done loading annotations in {:?} -- #transcripts = {}, #hgnc_ids = {}",
        start.elapsed(),
        merged
            .transcript_id_to_transcript
            .len()
            .separate_with_underscores(),
        merged
            .gene_id_to_transcript_ids
            .len()
            .separate_with_underscores()
    );

    Ok(merged)
}

/// Main entry point for `db create txs` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), Error> {
    fn _run(common: &crate::common::Args, args: &Args) -> Result<(), Error> {
        let mut report_file =
            File::create(format!("{}.report.jsonl", args.output.display())).map(BufWriter::new)?;
        let mut report = |r: ReportEntry| -> Result<(), Error> {
            writeln!(report_file, "{}", serde_json::to_string(&r)?)?;
            Ok(())
        };
        tracing::info!(
            "Building transcript and sequence database file\ncommon args: {:#?}\nargs: {:#?}",
            common,
            args
        );

        let mut tx_data = load_annotations(args)?;
        for (id, fix) in tx_data.fixes.iter() {
            report(ReportEntry::Fix(LogFix {
                source: "annotations".into(),
                fix: *fix,
                id: id.clone(),
                gene_name: tx_data.gene_name(id),
                tags: tx_data.tags(id),
            }))?;
        }

        let raw_tx_data = tx_data.clone();
        trace_rss_now();
        report(ReportEntry::Log(json!({
            "source": "annotations",
            "total_transcripts": tx_data.transcript_id_to_transcript.len(),
            "total_hgnc_ids": tx_data.gene_id_to_transcript_ids.len()
        })))?;

        // … then filter hgnc entries with no transcripts to boot …
        tx_data.filter_initial_hgnc_entries()?;
        // … then filter genes (missing hgnc id and/or symbol) …
        tx_data.filter_genes()?;
        // … then filter transcripts …
        tx_data.filter_transcripts()?;
        // … ensure there are no hgnc keys without associated transcripts left …
        tx_data.filter_empty_hgnc_mappings()?;

        // Open seqrepo / FASTA …
        let mut seq_provider = open_sequence_provider(args)?;
        // … and filter transcripts based on their sequences,
        // e.g. checking whether their translation contains a stop codon …
        let mut sequence_map = tx_data.filter_transcripts_with_sequence(&mut seq_provider)?;
        tx_data.filter_empty_hgnc_mappings()?;
        // … if there are genes with no transcripts left, check whether they are pseudogenes …
        tx_data.update_pseudogene_status()?;

        // … trigger the discard routine, but do not remove anything, just make sure they are consistent,
        // and report the stats.
        let remove = false;
        tx_data.discard(remove)?;

        // … and update all discard annotations …
        tx_data.propagate_discard_reasons(&raw_tx_data)?;

        trace_rss_now();

        report(ReportEntry::Log(json!({
            "source": "annotations_filtered",
            "total_transcripts": tx_data.transcript_id_to_transcript.len(),
            "total_hgnc_ids": tx_data.gene_id_to_transcript_ids.len()
        })))?;

        // For some stats, count number of chrMt, MANE Select and MANE Plus Clinical transcripts.
        let (n_mt, n_mane_select, n_mane_plus_clinical) = tx_data.gather_transcript_stats()?;

        report(ReportEntry::Log(json!({
            "source": "annotations_filtered",
            "n_mt": n_mt,
            "n_mane_select": n_mane_select,
            "n_mane_plus_clinical": n_mane_plus_clinical
        })))?;

        let source_version = bundle_source_version_information(args);

        // … and finally construct protobuf txdb data structures.
        let tx_db = tx_data.build_protobuf(&mut sequence_map, source_version)?;

        // List all discarded transcripts and genes.
        for (id, reason) in tx_data.discards.into_iter().sorted_unstable() {
            if reason.intersects(Reason::hard()) {
                // if it has at least one hard reason → actually discarded.
                report(ReportEntry::Discard(Discard {
                    source: "protobuf".into(),
                    reason,
                    id: id.clone(),
                    gene_name: raw_tx_data.gene_name(&id),
                    tags: raw_tx_data.tags(&id),
                }))?;
            } else {
                // only soft reasons → kept but flagged
                report(ReportEntry::SoftFilter(SoftFilter {
                    source: "protobuf".into(),
                    reason,
                    id: id.clone(),
                    gene_name: raw_tx_data.gene_name(&id),
                    tags: raw_tx_data.tags(&id),
                }))?;
            }
        }
        trace_rss_now();

        write_tx_db(tx_db, &args.output, args.compression_level)?;

        tracing::info!("Done building transcript and sequence database file");
        Ok(())
    }

    fn bundle_source_version_information(args: &Args) -> SourceVersion {
        let assembly = args.assembly.clone();

        let assembly_version = args.assembly_version.clone();

        let source_name = args.transcript_source.clone();

        let source_version = args.transcript_source_version.clone().unwrap_or("".into());
        let annotation_version = args.annotation_version.clone().unwrap_or("".into());
        let annotation_name = args
            .annotation
            .iter()
            .map(|p| p.file_stem().unwrap().to_string_lossy())
            .collect::<Vec<_>>()
            .join(",");

        SourceVersion {
            mehari_version: crate::common::version().to_string(),
            assembly: assembly.to_string(),
            assembly_version,
            source_name: source_name.to_string(),
            source_version,
            annotation_version,
            annotation_name,
        }
    }

    let threadpool = rayon::ThreadPoolBuilder::default()
        .num_threads(args.threads)
        .build()?;

    threadpool.install(|| _run(common, args))
}

pub(crate) fn write_tx_db(
    tx_db: TxSeqDatabase,
    path: impl AsRef<Path>,
    compression_level: i32,
) -> Result<(), Error> {
    tracing::info!("Writing out final database …");
    let path = path.as_ref();
    let mut buf = prost::bytes::BytesMut::with_capacity(tx_db.encoded_len());
    tx_db
        .encode(&mut buf)
        .map_err(|e| anyhow!("failed to encode: {}", e))?;
    tracing::info!("  … done constructing final tx and seq database");

    // Write out the final transcript and sequence database.
    tracing::info!("  Writing out final database …");
    // Open file and if necessary, wrap in a decompressor.
    let file = std::fs::File::create(path)
        .map_err(|e| anyhow!("failed to create file {}: {}", path.display(), e))?;
    let ext = path.extension().map(|s| s.to_str());
    let mut writer: Box<dyn Write> = if ext == Some(Some("zst")) {
        Box::new(
            zstd::Encoder::new(file, compression_level)
                .map_err(|e| anyhow!("failed to open zstd encoder for {}: {}", path.display(), e))?
                .auto_finish(),
        )
    } else {
        Box::new(file)
    };
    writer
        .write_all(&buf)
        .map_err(|e| anyhow!("failed to write to {}: {}", path.display(), e))?;
    tracing::info!("  … done writing out final database");

    Ok(())
}

#[cfg(test)]
pub mod test {
    use std::path::{Path, PathBuf};

    use clap_verbosity_flag::Verbosity;
    use itertools::Itertools;
    use rstest::rstest;
    use temp_testdir::TempDir;

    use crate::common::Args as CommonArgs;
    use crate::db::create::{GeneId, TranscriptLoader};
    use crate::db::dump;

    use super::{Args, run};

    #[test]
    fn filter_transcripts_brca1() -> Result<(), anyhow::Error> {
        let path_tsv = Path::new("tests/data/db/create/txs/txs_main.tsv");
        let mut tx_data = TranscriptLoader::new("GRCh37".to_string(), false);
        let labels = super::txid_to_label(path_tsv)?;
        tx_data.load_cdot(Path::new(
            "tests/data/db/create/txs/cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json",
        ))?;
        tx_data.apply_fixes(&Some(labels));

        eprintln!("{:#?}", &tx_data.gene_id_to_transcript_ids);
        insta::assert_yaml_snapshot!(
            tx_data
                .gene_id_to_transcript_ids
                .get(&GeneId::Hgnc(1100))
                .unwrap()
                .iter()
                .map(|s| s.as_str())
                .sorted_unstable()
                .collect::<Vec<_>>()
        );

        tx_data.filter_transcripts()?;
        insta::assert_yaml_snapshot!(
            tx_data
                .gene_id_to_transcript_ids
                .get(&GeneId::Hgnc(1100))
                .unwrap()
                .iter()
                .map(|s| s.as_str())
                .sorted_unstable()
                .collect::<Vec<_>>()
        );

        insta::assert_snapshot!(&tx_data.annotation_version);

        Ok(())
    }

    #[rstest]
    #[case("grch37")]
    #[case("grch38")]
    fn run_smoke_brca1_opa1(#[case] assembly: String) -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            output: tmp_dir.join("out.bin.zst"),
            annotation: vec![PathBuf::from(
                "tests/data/db/create/txs/cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json",
            )],
            mane_transcripts: Some(PathBuf::from("tests/data/db/create/txs/txs_main.tsv")),
            seqrepo: Some(PathBuf::from("tests/data/db/create/txs/latest")),
            transcript_sequences: None,
            assembly: assembly.clone(),
            assembly_version: None,
            transcript_source: "refseq".to_string(),
            transcript_source_version: None,
            disable_filters: false,
            threads: 1,
            annotation_version: Some("0.2.22".to_string()),
            compression_level: 19,
        };

        run(&common_args, &args)?;

        let mut buf: Vec<u8> = Vec::new();
        dump::run_with_write(
            &Default::default(),
            &dump::Args {
                path_db: tmp_dir.join("out.bin.zst"),
            },
            &mut buf,
        )?;
        crate::common::set_snapshot_suffix!("{}", assembly.to_lowercase());
        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }

    #[test]
    fn run_smoke_selenoproteins() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            output: tmp_dir.join("out.bin.zst"),
            annotation: vec![PathBuf::from(
                "tests/data/db/create/seleonoproteins/cdot-0.2.22.refseq.grch38.selenon.json",
            )],
            mane_transcripts: Some(PathBuf::from("tests/data/db/create/txs/txs_main.tsv")),
            seqrepo: Some(PathBuf::from("tests/data/db/create/seleonoproteins/latest")),
            transcript_sequences: None,
            assembly: "grch38".to_string(),
            assembly_version: None,
            transcript_source: "refseq".to_string(),
            transcript_source_version: None,
            disable_filters: false,
            threads: 1,
            annotation_version: Some("0.2.22".to_string()),
            compression_level: 19,
        };

        run(&common_args, &args)?;

        let mut buf: Vec<u8> = Vec::new();
        dump::run_with_write(
            &Default::default(),
            &dump::Args {
                path_db: tmp_dir.join("out.bin.zst"),
            },
            &mut buf,
        )?;
        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[test]
    fn run_smoke_mitochondrial() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(5, 0),
        };
        let args = Args {
            output: tmp_dir.join("out.bin.zst"),
            annotation: vec![PathBuf::from(
                "tests/data/db/create/mitochondrial/cdot-0.2.23.ensembl.chrMT.grch37.gff3.json",
            )],
            mane_transcripts: None,
            seqrepo: Some(PathBuf::from("tests/data/db/create/mitochondrial/latest")),
            transcript_sequences: None,
            assembly: "grch37".to_string(),
            assembly_version: None,
            transcript_source: "ensembl".to_string(),
            transcript_source_version: Some("98".into()),
            disable_filters: false,
            threads: 1,
            annotation_version: Some("0.2.23".to_string()),
            compression_level: 19,
        };

        run(&common_args, &args)?;

        let mut buf: Vec<u8> = Vec::new();
        dump::run_with_write(
            &Default::default(),
            &dump::Args {
                path_db: tmp_dir.join("out.bin.zst"),
            },
            &mut buf,
        )?;
        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }
}
