//! Data structures for storing transcripts after deserializing from flatbuffers.
//!
//! The structure is taken 1:1 from the flatbuffers IDL so in the future it may become more
//! Rust-y if we switch to bincode.

use enumset::{EnumSet, EnumSetType};
use serde::{Deserialize, Serialize};

// Enumeration for `Transcript::biotype`.
#[derive(Debug, Serialize, Deserialize)]
pub enum TranscriptBiotype {
    Coding,
    NonCoding,
}

/// Transcript tag enumeration.
#[derive(EnumSetType, Debug, Serialize, Deserialize)]
pub enum TranscriptTag {
    Basic,
    EnsemblCanonical,
    ManeSelect,
    ManePlusClinical,
    RefSeqSelect,
}

/// Enumeration for the known genome builds.
#[derive(Debug, Serialize, Deserialize)]
enum GenomeBuild {
    Grch37,
    Grch38,
}

/// Enumeration for the two strands of the genome.
#[derive(Debug, Serialize, Deserialize)]
enum Strand {
    Plus,
    Minus,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ExonAlignment {
    /// Start position on reference.
    alt_start_i: i32,
    /// End position on reference.
    alt_end_i: i32,
    /// Exon number.
    ord: i32,
    /// CDS start coordinate.
    alt_cds_start_i: i32,
    /// CDS end coordinate.
    alt_cds_end_i: i32,
    /// CIGAR string of alignment, empty indicates full matches.
    cigar: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct GenomeAlignment {
    /// The genome build identifier.
    genome_build: GenomeBuild,
    /// Accession of the contig sequence.
    contig: String,
    /// CDS end position, `-1` to indicate `None`.
    cds_start: i32,
    /// CDS end position, `-1` to indicate `None`.
    cds_end: i32,
    /// The strand.
    strand: Strand,
    /// Exons of the alignment.
    exons: Vec<ExonAlignment>,
}
// Store information about a transcript.
#[derive(Debug, Serialize, Deserialize)]
pub struct Transcript {
    // Transcript accession with version, e.g., `"NM_007294.3"` or `"ENST00000461574.1"` for BRCA1.
    pub id: String,
    // HGNC symbol, e.g., `"BRCA1"`
    pub gene_name: String,
    // HGNC gene identifier, e.g., `"1100"` for BRCA1.
    pub gene_id: String,
    // Transcript biotype.
    pub biotype: TranscriptBiotype,
    // Transcript flags, values from `TranscriptTag`, stored as OR-ed ubyte values.
    pub tags: EnumSet<TranscriptTag>,
    // Identifier of the corresponding protein.
    pub protein: String,
    // CDS start codon.
    pub start_codon: i32,
    // CDS stop codon.
    pub stop_codon: i32,
    // Alignments on the different genome builds.
    pub genome_alignments: Vec<GenomeAlignment>,
}

// Mapping from gene to transcript ID.
#[derive(Debug, Serialize, Deserialize)]
pub struct GeneToTxId {
    /// Gene HGNC symbol, serves as gene identifier.
    pub gene_name: String,
    /// Vector of all transcript IDs.
    pub tx_ids: Vec<String>,
}

/// Container for the transcript-related database.
#[derive(Debug, Serialize, Deserialize)]
pub struct TranscriptDb {
    /// Vector of all transcripts.
    pub transcripts: Vec<Transcript>,
    /// Mapping from gene ID to vector of all transcript IDs.
    pub gene_to_tx: Vec<GeneToTxId>,
}

/// Stores long array of sequences with an "index" of sequence names to their
/// index.
///
/// The fields `aliases` and `aliases_idx` have the same length and `aliases_idx[i]`
/// stores the index into `seqs` for the sequence `aliases[i]`.  In other words.
/// `seqs[aliases_idx[i]]` stores the sequence for `aliases[i]`.
#[derive(Debug, Serialize, Deserialize)]
pub struct SequenceDb {
    /// The sequence aliases, cf. `aliases_idx`.
    pub aliases: Vec<String>,
    /// The corresponding index in `seqs`, cf. `aliases`.
    pub aliases_idx: Vec<u32>,
    /// The corresponding sequences.
    pub seqs: Vec<String>,
}

//// Database of transcripts with sequences.
#[derive(Debug, Serialize, Deserialize)]
struct TxSeqDatabase {
    /// Store transcripts with their aliases.
    pub tx_db: TranscriptDb,
    /// Store sequence with their aliases.
    pub seq_db: SequenceDb,
}
