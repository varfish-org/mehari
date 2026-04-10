use crate::db::create::DisplayFromStr;
use anyhow::{Error, anyhow};
use enumflags2::BitFlags;
use enumflags2::bitflags;
use hgvs::data::cdot::json::models::Tag;
use nutype::nutype;
use serde::Serialize;
use serde_with::serde_as;
use std::fmt::{Display, Formatter};
use std::str::FromStr;

/// Helper struct for parsing the label TSV file.
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize)]
pub struct LabelEntry {
    /// Transcript identifier without version.
    pub(crate) transcript_id: String,
    /// Transcript version.
    transcript_version: usize,
    /// Gene symbol (unused).
    _gene_symbol: String,
    /// Label to transfer.
    pub(crate) label: String,
}

#[nutype(
    sanitize(trim),
    validate(not_empty),
    derive(
        Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, AsRef, Deref, Borrow, Into,
        Display
    )
)]
pub struct TranscriptId(String);

impl TranscriptId {
    pub(crate) fn without_version(&self) -> Result<Self, Error> {
        Ok(Self::try_new(
            self.as_ref()
                .split('.')
                .next()
                .ok_or(anyhow!("No version to split off"))?,
        )?)
    }

    pub(crate) fn split_version(&self) -> (&str, u32) {
        let (ac, version) = self.rsplit_once('.').unwrap_or_else(|| {
            panic!("Invalid accession, expected format 'ac.version', got {self}")
        });
        (ac, version.parse::<u32>().expect("invalid version"))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Hash)]
#[serde(tag = "type", content = "value")]
pub enum Identifier {
    Gene(GeneId),
    Transcript(TranscriptId),
}

#[bitflags]
#[repr(u32)]
#[derive(Debug, Clone, Copy, Serialize, Hash, PartialEq, Eq)]
pub enum Reason {
    Biotype,
    CdsEndAfterSequenceEnd,
    DeselectedGene,
    EmptyGenomeBuilds,
    FivePrimeEndTruncated,
    InvalidCdsLength,
    MissingGene,
    MissingGeneSymbol,
    MissingGeneId,
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
pub enum Fix {
    Cds,
    GenomeBuild,
    Tags,
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "type", content = "value")]
pub enum ReportEntry {
    Discard(Discard),
    SoftFilter(SoftFilter),
    Fix(LogFix),
    Log(serde_json::Value),
}

#[serde_as]
#[derive(Debug, Clone, Serialize)]
pub struct Discard {
    pub source: String,
    #[serde_as(as = "DisplayFromStr")]
    pub reason: BitFlags<Reason>,
    pub id: Identifier,
    pub gene_name: Option<String>,
    pub tags: Option<Vec<Tag>>,
}

#[serde_as]
#[derive(Debug, Clone, Serialize)]
pub struct SoftFilter {
    pub source: String,
    #[serde_as(as = "DisplayFromStr")]
    pub reason: BitFlags<Reason>,
    pub id: Identifier,
    pub gene_name: Option<String>,
    pub tags: Option<Vec<Tag>>,
}

#[serde_as]
#[derive(Debug, Clone, Serialize)]
pub struct LogFix {
    pub source: String,
    #[serde_as(as = "DisplayFromStr")]
    pub fix: BitFlags<Fix>,
    pub id: Identifier,
    pub gene_name: Option<String>,
    pub tags: Option<Vec<Tag>>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub enum GeneId {
    Hgnc(usize),
    Gene(String),
}

impl Display for GeneId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            GeneId::Hgnc(id) => write!(f, "HGNC:{}", id),
            GeneId::Gene(id) => write!(f, "GENE:{}", id),
        }
    }
}

impl FromStr for GeneId {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(num_str) = s.strip_prefix("HGNC:") {
            num_str.parse::<usize>().map(GeneId::Hgnc)
        } else if let Some(gene_str) = s.strip_prefix("GENE:") {
            Ok(GeneId::Gene(gene_str.to_string()))
        } else {
            s.parse::<usize>().map(GeneId::Hgnc)
        }
    }
}
