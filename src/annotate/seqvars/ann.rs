//! Code for annotating variants based on molecular consequence.
use parse_display::{Display, FromStr};

/// Putative impact level.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Display, FromStr)]
#[display(style = "UPPERCASE")]
pub enum PutativeImpact {
    High,
    Moderate,
    Low,
    Modifier,
}

/// Putative impact.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Display, FromStr)]
#[display(style = "snake_case")]
pub enum Consequence {
    // high impact
    ChromosomeNumberVariation,
    ExonLossVariant,
    FrameshiftVariant,
    RareAminoAcidVariant,
    SpliceAcceptorVariant,
    SpliceDonorVariant,
    StartLost,
    StopGained,
    StopLost,
    TranscriptAblation,
    // moderate impact
    #[display("3_prime_UTR_truncation")]
    ThreePrimeUtrTruncation,
    #[display("5_prime_UTR_truncation")]
    FivePrimeUtrTruncaction,
    CodingSequenceVairant,
    ConservativeInframeDeletion,
    ConservativeInframeInsertion,
    DisruptiveInframeDeletion,
    DisruptiveInframeInsertion,
    MissenseVariant,
    RegulatoryRegionAblation,
    SpliceRegionVariant,
    #[display("TFBS_ablation")]
    TbfsAblation,
    // low impact
    #[display("5_prime_UTR_premature_start_codon_gain_variant")]
    FivePrimeUtrPrematureStartCodonGainVariant,
    InitiatorCodonVariant,
    StartRetained,
    StopRetainedVariant,
    SynonymousVariant,
    // modifier
    #[display("3_prime_UTR_variant")]
    ThreePrimeUtrVariant,
    #[display("5_prime_UTR_variant")]
    FivePrimeUtrVariant,
    CodingSequenceVariant,
    ConservedIntergenicVariant,
    ConservedIntronVariant,
    DownstreamGeneVariant,
    ExonVariant,
    FeatureElongation,
    FeatureTruncation,
    GeneVariant,
    IntergenicVegion,
    IntragenicVariant,
    IntronVariant,
    #[display("mature_miRNA_variant")]
    MatureMirnaVariant,
    #[display("miRNA")]
    Mirna,
    #[display("NMD_transcript_variant")]
    NmdTranscriptVariant,
    NonCodingTranscriptExonVariant,
    NonCodingTranscriptVariant,
    RegulatoryRegionAmplification,
    RegulatoryRegionVariant,
    #[display("TF_binding_site_variant")]
    TfBindingSiteVariant,
    #[display("TFBS_amplification")]
    TfbsAmplification,
    TranscriptAmplification,
    TranscriptVariant,
    UpstreamGeneVariant,
}

impl Consequence {
    pub fn to_level(&self) -> PutativeImpact {
        match self {
            Consequence::ChromosomeNumberVariation
            | Consequence::ExonLossVariant
            | Consequence::FrameshiftVariant
            | Consequence::RareAminoAcidVariant
            | Consequence::SpliceAcceptorVariant
            | Consequence::SpliceDonorVariant
            | Consequence::StartLost
            | Consequence::StopGained
            | Consequence::StopLost
            | Consequence::TranscriptAblation => PutativeImpact::High,
            Consequence::ThreePrimeUtrTruncation
            | Consequence::FivePrimeUtrTruncaction
            | Consequence::CodingSequenceVairant
            | Consequence::ConservativeInframeDeletion
            | Consequence::ConservativeInframeInsertion
            | Consequence::DisruptiveInframeDeletion
            | Consequence::DisruptiveInframeInsertion
            | Consequence::MissenseVariant
            | Consequence::RegulatoryRegionAblation
            | Consequence::SpliceRegionVariant
            | Consequence::TbfsAblation => PutativeImpact::Moderate,
            Consequence::FivePrimeUtrPrematureStartCodonGainVariant
            | Consequence::InitiatorCodonVariant
            | Consequence::StartRetained
            | Consequence::StopRetainedVariant
            | Consequence::SynonymousVariant => PutativeImpact::Low,
            Consequence::ThreePrimeUtrVariant
            | Consequence::FivePrimeUtrVariant
            | Consequence::CodingSequenceVariant
            | Consequence::ConservedIntergenicVariant
            | Consequence::ConservedIntronVariant
            | Consequence::DownstreamGeneVariant
            | Consequence::ExonVariant
            | Consequence::FeatureElongation
            | Consequence::FeatureTruncation
            | Consequence::GeneVariant
            | Consequence::IntergenicVegion
            | Consequence::IntragenicVariant
            | Consequence::IntronVariant
            | Consequence::MatureMirnaVariant
            | Consequence::Mirna
            | Consequence::NmdTranscriptVariant
            | Consequence::NonCodingTranscriptExonVariant
            | Consequence::NonCodingTranscriptVariant
            | Consequence::RegulatoryRegionAmplification
            | Consequence::RegulatoryRegionVariant
            | Consequence::TfBindingSiteVariant
            | Consequence::TfbsAmplification
            | Consequence::TranscriptAmplification
            | Consequence::TranscriptVariant
            | Consequence::UpstreamGeneVariant => PutativeImpact::Modifier,
        }
    }
}

/// Enumeration for `AnnField::allele`.
pub enum Allele {
    /// A simple value for the allele.
    Alt { alternative: String },
    /// Allele encoding in the case of cancer samples (when an ALT is the reference).
    ///
    /// For example, if `T` is the major allele, the donor carriers `C` and there is
    /// a somatic mutation to `G`: `G-C`.
    AltRef {
        alternative: String,
        reference: String,
    },
    /// Compound annotation of other variant in the annotation.
    Compound {
        alternative: String,
        other_chrom: String,
        other_pos: u32,
        other_ref: String,
        other_alt: String,
    },
}

/// Sequence ontology feature.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Display, FromStr)]
#[display(style = "snake_case")]
pub enum SoFeature {
    Transcript,
}

/// Enum for `AnnField::feature_type`.
pub enum FeatureType {
    SoTerm { term: SoFeature },
    Custom { value: String },
}

/// Encode feature biotype.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Display, FromStr)]
pub enum FeatureBiotype {
    Coding,
    Noncoding,
}

/// Encode exon/intron rank.
pub struct Rank {
    pub ord: u32,
    pub total: u32,
}

/// Position, optionally with total length.
pub struct Pos {
    pub ord: u32,
    pub total: Option<u32>,
}

/// A message to be used in `AnnField::messages`.
pub enum Message {
    ErrorChromosomeNotFound,
    ErrorOutOfChromosomeRange,
    WarningRefDoesNotMatchGenome,
    WarningSequnenceNotAvailable,
    WarningTranscriptIncomplete,
    WarningTranscriptMultipleStopCodons,
    WarningTranscriptsNoStartCodon,
    InfoRealignThreePrime,
    InfoCompoundAnnotation,
    InfoNonReferenceAnnotation,
}

/// Representation of an `ANN` field.
pub struct AnnField {
    /// The alternative allele that this annotation refers to.
    pub allele: Allele,
    /// The consequences of the allele.
    pub consequences: Vec<Consequence>,
    /// The putative impact.
    pub putative_impact: PutativeImpact,
    /// The gene identifier.
    pub gene_id: String,
    /// The feature type.
    pub feature_type: FeatureType,
    /// The feature identifier.
    pub feature_id: String,
    /// The feature biotype.
    pub feature_biotype: FeatureBiotype,
    /// The exon / intron rank.
    pub rank: Rank,
    /// HGVS c. notation.
    pub hgvs_c: String,
    /// HGVS p. notation.
    pub hgvs_p: String,
    /// cDNA position.
    pub cdna_pos: Pos,
    /// CDS position.
    pub cds_pos: Pos,
    /// Protein position.
    pub protein_pos: Pos,
    /// Distance to feature.
    pub distance: Option<i32>,
    /// Optional list of warnings and error messages.
    pub messages: Option<Vec<Message>>,
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use super::*;

    #[test]
    fn putative_impact_display() {
        assert_eq!(format!("{}", PutativeImpact::High), "HIGH");
    }
}
