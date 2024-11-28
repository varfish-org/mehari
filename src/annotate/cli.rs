use clap::Args as ClapArgs;
use strum::{Display, VariantArray};

#[derive(Debug, ClapArgs)]
#[group(required = true, multiple = true)]
pub struct Sources {
    #[arg(long)]
    pub transcripts: Option<Vec<String>>,
    #[arg(long)]
    pub frequencies: Option<String>,
    #[arg(long)]
    pub clinvar: Option<String>,
}

#[derive(Debug, ClapArgs)]
pub struct TranscriptSettings {
    /// The transcript source.
    #[arg(long, value_enum, default_value_t = TranscriptSource::Both)]
    pub transcript_source: TranscriptSource,

    /// Whether to report only the worst consequence for each picked transcript.
    #[arg(long)]
    pub report_most_severe_consequence_by: Option<ConsequenceBy>,

    /// Which kind of transcript to pick / restrict to. Default is not to pick at all.
    ///
    /// Depending on `--pick-transcript-mode`, if multiple transcripts match the selection,
    /// either the first one is kept or all are kept.
    #[arg(long)]
    pub pick_transcript: Vec<TranscriptPickType>,

    /// Determines how to handle multiple transcripts. Default is to keep all.
    ///
    /// When transcript picking is enabled via `--pick-transcript`,
    /// either keep the first one found or keep all that match.
    #[arg(long, default_value = "all")]
    pub pick_transcript_mode: TranscriptPickMode,
}

#[derive(
    Debug,
    Copy,
    Clone,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Display,
    clap::ValueEnum,
    VariantArray,
    parse_display::FromStr,
)]
pub enum ConsequenceBy {
    Gene,
    Transcript,
    // or "Variant"?
    Allele,
}

#[derive(
    Debug,
    Copy,
    Clone,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Display,
    clap::ValueEnum,
    VariantArray,
    parse_display::FromStr,
)]
pub enum TranscriptPickType {
    ManeSelect,
    ManePlusClinical,
    Length,
    EnsemblCanonical,
    RefSeqSelect,
    GencodePrimary,
    Basic,
}

#[derive(Debug, Copy, Clone, Display, clap::ValueEnum, Default)]
pub enum TranscriptPickMode {
    #[default]
    First,
    All,
}

/// Enum that allows to select the transcript source.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Default,
    serde::Deserialize,
    serde::Serialize,
    clap::ValueEnum,
)]
pub enum TranscriptSource {
    /// ENSEMBL
    Ensembl,
    /// RefSeq
    RefSeq,
    /// Both
    #[default]
    Both,
}
