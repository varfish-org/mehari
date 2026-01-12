use clap::Args as ClapArgs;
use strum::{Display, VariantArray};

#[derive(Debug, ClapArgs)]
#[group(required = true, multiple = true)]
pub struct Sources {
    /// Transcript database containing the transcript information.
    ///
    /// Pre-built databases are available at https://github.com/varfish-org/mehari-data-tx/releases
    #[arg(long)]
    pub transcripts: Option<Vec<String>>,

    /// Frequency database.
    ///
    /// The frequency database contains gnomAD frequencies for the variants.
    /// Pre-built databases are available at TODO
    #[arg(long)]
    pub frequencies: Option<Vec<String>>,

    /// ClinVar database.
    ///
    /// The ClinVar database contains clinical significance information for the variants.
    /// Pre-built databases are available at https://github.com/varfish-org/annonars-data-clinvar/releases
    #[arg(long)]
    pub clinvar: Option<Vec<String>>,
}

#[derive(Debug, ClapArgs, Default, Clone)]
pub struct TranscriptSettings {
    /// The transcript source.
    #[arg(long, value_enum, default_value_t = TranscriptSource::Both)]
    pub transcript_source: TranscriptSource,

    /// Whether to report only the most severe consequence, grouped by gene, transcript, or allele.
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

    /// Whether to keep intergenic variants.
    #[arg(long, default_value = "false")]
    pub keep_intergenic: bool,

    /// Whether to report splice variants in UTRs.
    #[arg(long, default_value = "false")]
    pub discard_utr_splice_variants: bool,

    /// Whether to do hgvs shifting for hgvs.g like vep does
    #[arg(long, default_value = "false", hide = true)]
    pub vep_hgvs_shift: bool,
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
    ManeSelectBackport,
    ManePlusClinical,
    ManePlusClinicalBackport,
    Length,
    EnsemblCanonical,
    EnsemblCanonicalBackport,
    RefSeqSelect,
    RefSeqSelectBackport,
    GencodePrimary,
    GencodePrimaryBackport,
    Basic,
    BasicBackport,
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
    Hash,
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
