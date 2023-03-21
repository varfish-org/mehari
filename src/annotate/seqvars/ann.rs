//! Code for annotating variants based on molecular consequence.
use std::str::FromStr;

use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alphanumeric1, digit1},
    combinator::{all_consuming, map},
    sequence::tuple,
    IResult,
};
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
    pub fn to_impact(&self) -> PutativeImpact {
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

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Display)]
/// Enumeration for `AnnField::allele`.
pub enum Allele {
    /// A simple value for the allele.
    #[display("{alternative}")]
    Alt { alternative: String },
    /// Allele encoding in the case of cancer samples (when an ALT is the reference).
    ///
    /// For example, if `T` is the major allele, the donor carriers `C` and there is
    /// a somatic mutation to `G`: `G-C`.
    #[display("{alternative}-{reference}")]
    AltRef {
        alternative: String,
        reference: String,
    },
    /// Compound annotation of other variant in the annotation.
    #[display("{alternative}-{other_chrom}:{other_pos}_{other_ref}>{other_alt}")]
    Compound {
        alternative: String,
        other_chrom: String,
        other_pos: u32,
        other_ref: String,
        other_alt: String,
    },
}

mod parse {
    use nom::bytes::complete::take_while1;

    pub static NA_IUPAC: &str = "ACGTURYMKWSBDHVNacgturymkwsbdhvn";

    pub fn na1(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        take_while1(|c: char| NA_IUPAC.contains(c))(input)
    }
}

impl Allele {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        all_consuming(alt((
            Self::parse_compound,
            Self::parse_alt_ref,
            Self::parse_alt,
        )))(input)
    }

    fn parse_compound(input: &str) -> IResult<&str, Self> {
        map(
            tuple((
                parse::na1,
                tag("-"),
                alphanumeric1,
                tag(":"),
                digit1,
                tag("_"),
                parse::na1,
                tag(">"),
                parse::na1,
            )),
            |(alternative, _, other_chrom, _, other_pos, _, other_ref, _, other_alt)| {
                Allele::Compound {
                    alternative: alternative.to_string(),
                    other_chrom: other_chrom.to_string(),
                    other_pos: other_pos.parse::<u32>().unwrap(),
                    other_ref: other_ref.to_string(),
                    other_alt: other_alt.to_string(),
                }
            },
        )(input)
    }

    fn parse_alt_ref(input: &str) -> IResult<&str, Self> {
        map(
            tuple((parse::na1, tag("-"), parse::na1)),
            |(alternative, _, reference)| Allele::AltRef {
                alternative: alternative.to_string(),
                reference: reference.to_string(),
            },
        )(input)
    }

    fn parse_alt(input: &str) -> IResult<&str, Self> {
        map(parse::na1, |alternative| Allele::Alt {
            alternative: alternative.to_string(),
        })(input)
    }
}

impl FromStr for Allele {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
            .map(|(_, value)| value)
            .map_err(|e| anyhow::anyhow!("{}", e))
    }
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
    use std::str::FromStr;

    use pretty_assertions::assert_eq;

    use super::*;

    #[test]
    fn putative_impact_display() {
        assert_eq!(format!("{}", PutativeImpact::High), "HIGH");
        assert_eq!(format!("{}", PutativeImpact::Moderate), "MODERATE");
        assert_eq!(format!("{}", PutativeImpact::Low), "LOW");
        assert_eq!(format!("{}", PutativeImpact::Modifier), "MODIFIER");
    }

    #[test]
    fn putative_impact_from_str() -> Result<(), anyhow::Error> {
        assert_eq!(PutativeImpact::from_str("HIGH")?, PutativeImpact::High);
        assert_eq!(
            PutativeImpact::from_str("MODERATE")?,
            PutativeImpact::Moderate
        );
        assert_eq!(PutativeImpact::from_str("LOW")?, PutativeImpact::Low);
        assert_eq!(
            PutativeImpact::from_str("MODIFIER")?,
            PutativeImpact::Modifier
        );

        Ok(())
    }

    #[test]
    fn consequence_display() {
        assert_eq!(
            format!("{}", Consequence::ChromosomeNumberVariation),
            "chromosome_number_variation"
        );
        assert_eq!(
            format!("{}", Consequence::ThreePrimeUtrTruncation),
            "3_prime_UTR_truncation"
        );
        assert_eq!(
            format!("{}", Consequence::FivePrimeUtrTruncaction),
            "5_prime_UTR_truncation"
        );
        assert_eq!(format!("{}", Consequence::TbfsAblation), "TFBS_ablation");
        assert_eq!(
            format!(
                "{}",
                Consequence::FivePrimeUtrPrematureStartCodonGainVariant
            ),
            "5_prime_UTR_premature_start_codon_gain_variant"
        );
        assert_eq!(
            format!("{}", Consequence::ThreePrimeUtrVariant),
            "3_prime_UTR_variant"
        );
        assert_eq!(
            format!("{}", Consequence::FivePrimeUtrVariant),
            "5_prime_UTR_variant"
        );
        assert_eq!(
            format!("{}", Consequence::MatureMirnaVariant),
            "mature_miRNA_variant"
        );
        assert_eq!(format!("{}", Consequence::Mirna), "miRNA");
        assert_eq!(
            format!("{}", Consequence::NmdTranscriptVariant),
            "NMD_transcript_variant"
        );
        assert_eq!(
            format!("{}", Consequence::TfBindingSiteVariant),
            "TF_binding_site_variant"
        );
        assert_eq!(
            format!("{}", Consequence::TfbsAmplification),
            "TFBS_amplification"
        );
    }

    #[test]
    fn consequence_from_str() -> Result<(), anyhow::Error> {
        assert_eq!(
            Consequence::from_str("chromosome_number_variation")?,
            Consequence::ChromosomeNumberVariation
        );
        assert_eq!(
            Consequence::from_str("3_prime_UTR_truncation")?,
            Consequence::ThreePrimeUtrTruncation,
        );
        assert_eq!(
            Consequence::from_str("5_prime_UTR_truncation")?,
            Consequence::FivePrimeUtrTruncaction,
        );
        assert_eq!(
            Consequence::from_str("TFBS_ablation")?,
            Consequence::TbfsAblation,
        );
        assert_eq!(
            Consequence::from_str("5_prime_UTR_premature_start_codon_gain_variant")?,
            Consequence::FivePrimeUtrPrematureStartCodonGainVariant,
        );
        assert_eq!(
            Consequence::from_str("3_prime_UTR_variant")?,
            Consequence::ThreePrimeUtrVariant,
        );
        assert_eq!(
            Consequence::from_str("5_prime_UTR_variant")?,
            Consequence::FivePrimeUtrVariant,
        );
        assert_eq!(
            Consequence::from_str("mature_miRNA_variant")?,
            Consequence::MatureMirnaVariant,
        );
        assert_eq!(Consequence::from_str("miRNA")?, Consequence::Mirna,);
        assert_eq!(
            Consequence::from_str("NMD_transcript_variant")?,
            Consequence::NmdTranscriptVariant,
        );
        assert_eq!(
            Consequence::from_str("TF_binding_site_variant")?,
            Consequence::TfBindingSiteVariant,
        );
        assert_eq!(
            Consequence::from_str("TFBS_amplification")?,
            Consequence::TfbsAmplification,
        );

        Ok(())
    }

    #[test]
    fn consequence_to_impact() {
        assert_eq!(
            Consequence::ChromosomeNumberVariation.to_impact(),
            PutativeImpact::High,
        );
        assert_eq!(
            Consequence::ThreePrimeUtrTruncation.to_impact(),
            PutativeImpact::Moderate,
        );
        assert_eq!(
            Consequence::FivePrimeUtrPrematureStartCodonGainVariant.to_impact(),
            PutativeImpact::Low,
        );
        assert_eq!(
            Consequence::ThreePrimeUtrVariant.to_impact(),
            PutativeImpact::Modifier,
        );
    }

    #[test]
    fn allele_display() {
        assert_eq!(
            format!(
                "{}",
                Allele::Alt {
                    alternative: String::from("A")
                }
            ),
            "A",
        );
        assert_eq!(
            format!(
                "{}",
                Allele::AltRef {
                    alternative: String::from("A"),
                    reference: String::from("C"),
                }
            ),
            "A-C",
        );
        assert_eq!(
            format!(
                "{}",
                Allele::Compound {
                    alternative: String::from("A"),
                    other_chrom: String::from("chr1"),
                    other_pos: 123,
                    other_ref: String::from("C"),
                    other_alt: String::from("T"),
                }
            ),
            "A-chr1:123_C>T",
        );
    }

    #[test]
    fn allele_from_str() -> Result<(), anyhow::Error> {
        assert_eq!(
            Allele::from_str("A")?,
            Allele::Alt {
                alternative: String::from("A")
            },
        );
        assert_eq!(
            Allele::from_str("A-C")?,
            Allele::AltRef {
                alternative: String::from("A"),
                reference: String::from("C"),
            }
        );
        assert_eq!(
            Allele::from_str("A-chr1:123_C>T")?,
            Allele::Compound {
                alternative: String::from("A"),
                other_chrom: String::from("chr1"),
                other_pos: 123,
                other_ref: String::from("C"),
                other_alt: String::from("T"),
            },
        );

        Ok(())
    }
}
