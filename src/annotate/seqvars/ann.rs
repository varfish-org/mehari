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
use strum::IntoEnumIterator;

/// Putative impact level.
#[derive(
    Debug,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Clone,
    Copy,
    Display,
    FromStr,
    serde::Deserialize,
    serde::Serialize,
    strum::EnumIter,
)]
#[display(style = "UPPERCASE")]
#[serde(rename_all = "UPPERCASE")]
pub enum PutativeImpact {
    High,
    Moderate,
    Low,
    Modifier,
}

/// Putative impact.
#[derive(
    Debug,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Clone,
    Copy,
    Display,
    FromStr,
    serde::Deserialize,
    serde::Serialize,
    strum::EnumIter,
)]
#[display(style = "snake_case")]
#[serde(rename_all = "snake_case")]
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
    #[serde(rename = "3_prime_UTR_truncation")]
    ThreePrimeUtrTruncation,
    #[display("5_prime_UTR_truncation")]
    #[serde(rename = "5_prime_UTR_truncation")]
    FivePrimeUtrTruncaction,
    ConservativeInframeDeletion,
    ConservativeInframeInsertion,
    DisruptiveInframeDeletion,
    DisruptiveInframeInsertion,
    MissenseVariant,
    RegulatoryRegionAblation,
    SpliceRegionVariant,
    #[display("TFBS_ablation")]
    #[serde(rename = "TFBS_ablation")]
    TbfsAblation,
    // low impact
    #[display("5_prime_UTR_premature_start_codon_gain_variant")]
    #[serde(rename = "5_prime_UTR_premature_start_codon_gain_variant")]
    FivePrimeUtrPrematureStartCodonGainVariant,
    InitiatorCodonVariant,
    StartRetained,
    StopRetainedVariant,
    SynonymousVariant,
    // modifier
    #[display("3_prime_UTR_variant")]
    #[serde(rename = "3_prime_UTR_variant")]
    ThreePrimeUtrVariant,
    #[display("5_prime_UTR_variant")]
    #[serde(rename = "5_prime_UTR_variant")]
    FivePrimeUtrVariant,
    CodingSequenceVariant,
    ConservedIntergenicVariant,
    ConservedIntronVariant,
    DownstreamGeneVariant,
    ExonVariant,
    FeatureElongation,
    FeatureTruncation,
    GeneVariant,
    IntergenicVariant,
    IntronVariant,
    #[display("mature_miRNA_variant")]
    #[serde(rename = "mature_miRNA_variant")]
    MatureMirnaVariant,
    #[display("miRNA")]
    #[serde(rename = "miRNA")]
    Mirna,
    #[display("NMD_transcript_variant")]
    #[serde(rename = "NMD_transcript_variant")]
    NmdTranscriptVariant,
    NonCodingTranscriptExonVariant,
    NonCodingTranscriptIntronVariant,
    RegulatoryRegionAmplification,
    RegulatoryRegionVariant,
    #[display("TF_binding_site_variant")]
    #[serde(rename = "TF_binding_site_variant")]
    TfBindingSiteVariant,
    #[display("TFBS_amplification")]
    #[serde(rename = "TFBS_amplification")]
    TfbsAmplification,
    TranscriptAmplification,
    TranscriptVariant,
    UpstreamGeneVariant,
}

impl From<Consequence> for PutativeImpact {
    fn from(val: Consequence) -> Self {
        match val {
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
            | Consequence::IntergenicVariant
            | Consequence::IntronVariant
            | Consequence::MatureMirnaVariant
            | Consequence::Mirna
            | Consequence::NmdTranscriptVariant
            | Consequence::NonCodingTranscriptExonVariant
            | Consequence::NonCodingTranscriptIntronVariant
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

impl Consequence {
    /// Return vector of all values of `Consequence`.
    pub fn all() -> Vec<Self> {
        Self::iter().collect()
    }

    pub fn impact(&self) -> PutativeImpact {
        PutativeImpact::from(*self)
    }
}

#[derive(
    Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Display, serde::Deserialize, serde::Serialize,
)]
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
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Display,
    FromStr,
    serde::Deserialize,
    serde::Serialize,
)]
#[display(style = "snake_case")]
pub enum SoFeature {
    Transcript,
}

/// Enum for `AnnField::feature_type`.
#[derive(
    Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Display, serde::Deserialize, serde::Serialize,
)]
pub enum FeatureType {
    #[display("{term}")]
    SoTerm { term: SoFeature },
    #[display("{value}")]
    Custom { value: String },
}

impl FromStr for FeatureType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        SoFeature::from_str(s)
            .map(|term| FeatureType::SoTerm { term })
            .or_else(|_| {
                Ok(FeatureType::Custom {
                    value: s.to_string(),
                })
            })
    }
}

/// Encode feature biotype.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Display,
    FromStr,
    serde::Deserialize,
    serde::Serialize,
    strum::EnumIter,
)]
pub enum FeatureBiotype {
    Coding,
    Noncoding,
}

impl FeatureBiotype {
    pub fn is_coding(&self) -> bool {
        match self {
            FeatureBiotype::Coding => true,
            FeatureBiotype::Noncoding => false,
        }
    }
}

/// Encode exon/intron rank.
#[derive(
    Clone,
    Debug,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Display,
    FromStr,
    Default,
    serde::Deserialize,
    serde::Serialize,
)]
#[display("{ord}/{total}")]
pub struct Rank {
    pub ord: i32,
    pub total: i32,
}

/// Position, optionally with total length.
#[derive(
    Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, serde::Deserialize, serde::Serialize,
)]
pub struct Pos {
    pub ord: i32,
    pub total: Option<i32>,
}

impl std::fmt::Display for Pos {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(total) = self.total {
            write!(f, "{}/{}", self.ord, total)
        } else {
            write!(f, "{}", self.ord)
        }
    }
}

impl Pos {
    fn parse_number_neg(input: &str) -> IResult<&str, i32> {
        map(tuple((tag("-"), digit1::<&str, _>)), |(sign, num)| {
            let num = num.parse::<i32>().unwrap();
            if sign == "-" {
                -num
            } else {
                num
            }
        })(input)
    }

    fn parse_number_nosign(input: &str) -> IResult<&str, i32> {
        map(digit1::<&str, _>, |num| num.parse::<i32>().unwrap())(input)
    }

    fn parse_number(input: &str) -> IResult<&str, i32> {
        alt((Self::parse_number_neg, Self::parse_number_nosign))(input)
    }

    fn parse_with_total(input: &str) -> IResult<&str, Self> {
        map(
            tuple((Self::parse_number, tag("/"), digit1)),
            |(num, _, total)| Pos {
                ord: num,
                total: Some(total.parse::<i32>().unwrap()),
            },
        )(input)
    }

    fn parse_no_total(input: &str) -> IResult<&str, Self> {
        map(Self::parse_number, |num| Pos {
            ord: num,
            total: None,
        })(input)
    }
}

impl FromStr for Pos {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        all_consuming(alt((Self::parse_with_total, Self::parse_no_total)))(s)
            .map(|(_, value)| value)
            .map_err(|e| anyhow::anyhow!("{}", e))
    }
}

/// A message to be used in `AnnField::messages`.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Display,
    FromStr,
    serde::Deserialize,
    serde::Serialize,
)]
#[display(style = "SNAKE_CASE")]
pub enum Message {
    ErrorChromosomeNotFound,
    ErrorOutOfChromosomeRange,
    WarningRefDoesNotMatchGenome,
    WarningSequenceNotAvailable,
    WarningTranscriptIncomplete,
    WarningTranscriptMultipleStopCodons,
    WarningTranscriptsNoStartCodon,
    InfoRealignThreePrime,
    InfoCompoundAnnotation,
    InfoNonReferenceAnnotation,
}

/// Representation of an `ANN` field.
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct AnnField {
    /// The alternative allele that this annotation refers to.
    pub allele: Allele,
    /// The consequences of the allele.
    pub consequences: Vec<Consequence>,
    /// The putative impact.
    pub putative_impact: PutativeImpact,
    /// The gene symbol.
    pub gene_symbol: String,
    /// The gene identifier.
    pub gene_id: String,
    /// The feature type.
    pub feature_type: FeatureType,
    /// The feature identifier.
    pub feature_id: String,
    /// The feature biotype.
    pub feature_biotype: FeatureBiotype,
    /// The exon / intron rank.
    pub rank: Option<Rank>,
    /// HGVS c. notation.
    pub hgvs_t: Option<String>,
    /// HGVS p. notation.
    pub hgvs_p: Option<String>,
    /// cDNA position.
    pub tx_pos: Option<Pos>,
    /// CDS position.
    pub cds_pos: Option<Pos>,
    /// Protein position.
    pub protein_pos: Option<Pos>,
    /// Distance to feature.
    pub distance: Option<i32>,
    /// Optional list of warnings and error messages.
    pub messages: Option<Vec<Message>>,
}

impl FromStr for AnnField {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split('|');
        let allele = fields.next().unwrap().parse()?;
        let consequences = fields
            .next()
            .unwrap()
            .split('&')
            .map(|s| s.parse())
            .collect::<Result<Vec<_>, _>>()?;
        let putative_impact = fields.next().unwrap().parse()?;
        let gene_symbol = fields.next().unwrap().to_string();
        let gene_id = fields.next().unwrap().to_string();
        let feature_type = fields.next().unwrap().parse()?;
        let feature_id = fields.next().unwrap().to_string();
        let feature_biotype = fields.next().unwrap().parse()?;
        let rank = fields.next().unwrap();
        let rank = if rank.is_empty() {
            None
        } else {
            Some(rank.parse()?)
        };
        let hgvs_t = fields.next().unwrap();
        let hgvs_t = if hgvs_t.is_empty() {
            None
        } else {
            Some(hgvs_t.to_string())
        };
        let hgvs_p = fields.next().unwrap();
        let hgvs_p = if hgvs_p.is_empty() {
            None
        } else {
            Some(hgvs_p.to_string())
        };
        let tx_pos = fields.next().unwrap();
        let tx_pos = if tx_pos.is_empty() {
            None
        } else {
            Some(tx_pos.parse()?)
        };
        let cds_pos = fields.next().unwrap();
        let cds_pos = if cds_pos.is_empty() {
            None
        } else {
            Some(cds_pos.parse()?)
        };
        let protein_pos = fields.next().unwrap();
        let protein_pos = if protein_pos.is_empty() {
            None
        } else {
            Some(protein_pos.parse()?)
        };
        let distance = fields.next().unwrap();
        let distance = if distance.is_empty() {
            None
        } else {
            Some(distance.parse()?)
        };
        let messages = fields.next().unwrap();
        let messages = if messages.is_empty() {
            None
        } else {
            Some(
                messages
                    .split('&')
                    .map(|s| s.parse())
                    .collect::<Result<Vec<_>, _>>()?,
            )
        };
        Ok(AnnField {
            allele,
            consequences,
            putative_impact,
            gene_symbol,
            gene_id,
            feature_type,
            feature_id,
            feature_biotype,
            rank,
            hgvs_t,
            hgvs_p,
            tx_pos,
            cds_pos,
            protein_pos,
            distance,
            messages,
        })
    }
}

impl std::fmt::Display for AnnField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.allele)?;
        write!(f, "|")?;
        for (i, csq) in self.consequences.iter().enumerate() {
            if i > 0 {
                write!(f, "&{}", csq)?;
            } else {
                write!(f, "{}", csq)?;
            }
        }
        write!(f, "|")?;
        write!(f, "{}", self.putative_impact)?;
        write!(f, "|")?;
        write!(f, "{}", self.gene_symbol)?;
        write!(f, "|")?;
        write!(f, "{}", self.gene_id)?;
        write!(f, "|")?;
        write!(f, "{}", self.feature_type)?;
        write!(f, "|")?;
        write!(f, "{}", self.feature_id)?;
        write!(f, "|")?;
        write!(f, "{}", self.feature_biotype)?;
        write!(f, "|")?;
        if let Some(rank) = &self.rank {
            write!(f, "{}", rank)?;
        }
        write!(f, "|")?;
        if let Some(hgvs_c) = &self.hgvs_t {
            write!(f, "{}", hgvs_c)?;
        }
        write!(f, "|")?;
        if let Some(hgvs_p) = &self.hgvs_p {
            write!(f, "{}", hgvs_p)?;
        }
        write!(f, "|")?;
        if let Some(cdna_pos) = &self.tx_pos {
            write!(f, "{}", cdna_pos)?;
        }
        write!(f, "|")?;
        if let Some(cds_pos) = &self.cds_pos {
            write!(f, "{}", cds_pos)?;
        }
        write!(f, "|")?;
        if let Some(protein_pos) = &self.protein_pos {
            write!(f, "{}", protein_pos)?;
        }
        write!(f, "|")?;
        if let Some(distance) = self.distance {
            write!(f, "{}", distance)?;
        }
        write!(f, "|")?;
        if let Some(messages) = &self.messages {
            for (i, csq) in messages.iter().enumerate() {
                if i > 0 {
                    write!(f, "&{}", csq)?;
                } else {
                    write!(f, "{}", csq)?;
                }
            }
        }

        Ok(())
    }
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
        {
            let p: PutativeImpact = Consequence::ChromosomeNumberVariation.into();
            assert_eq!(p, PutativeImpact::High,);
        }
        {
            let p: PutativeImpact = Consequence::MissenseVariant.into();
            assert_eq!(p, PutativeImpact::Moderate,);
        }
        {
            let p: PutativeImpact = Consequence::SynonymousVariant.into();
            assert_eq!(p, PutativeImpact::Low,);
        }
        {
            let p: PutativeImpact = Consequence::UpstreamGeneVariant.into();
            assert_eq!(p, PutativeImpact::Modifier,);
        }
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

    #[test]
    fn so_feature_display_from_str() -> Result<(), anyhow::Error> {
        assert_eq!(format!("{}", SoFeature::Transcript), "transcript",);
        assert_eq!(SoFeature::from_str("transcript")?, SoFeature::Transcript,);

        Ok(())
    }

    #[test]
    fn feature_type_display() -> Result<(), anyhow::Error> {
        assert_eq!(
            format!(
                "{}",
                FeatureType::SoTerm {
                    term: SoFeature::Transcript
                }
            ),
            "transcript",
        );
        assert_eq!(
            format!(
                "{}",
                FeatureType::Custom {
                    value: String::from("foo")
                }
            ),
            "foo",
        );

        Ok(())
    }

    #[test]
    fn so_feature_from_str() -> Result<(), anyhow::Error> {
        assert_eq!(
            FeatureType::from_str("transcript")?,
            FeatureType::SoTerm {
                term: SoFeature::Transcript
            },
        );
        assert_eq!(
            FeatureType::from_str("foo")?,
            FeatureType::Custom {
                value: String::from("foo")
            },
        );

        Ok(())
    }

    #[test]
    fn feature_biotype() -> Result<(), anyhow::Error> {
        assert_eq!(FeatureBiotype::from_str("Coding")?, FeatureBiotype::Coding);
        assert_eq!(
            FeatureBiotype::from_str("Noncoding")?,
            FeatureBiotype::Noncoding
        );

        assert_eq!(format!("{}", FeatureBiotype::Coding), "Coding");
        assert_eq!(format!("{}", FeatureBiotype::Noncoding), "Noncoding");

        Ok(())
    }

    #[test]
    fn rank_display() -> Result<(), anyhow::Error> {
        assert_eq!(format!("{}", Rank { ord: 1, total: 2 }), "1/2",);

        Ok(())
    }

    #[test]
    fn rank_from_str() -> Result<(), anyhow::Error> {
        assert_eq!(Rank::from_str("1/2")?, Rank { ord: 1, total: 2 },);

        Ok(())
    }

    #[test]
    fn pos_display() -> Result<(), anyhow::Error> {
        assert_eq!(
            format!(
                "{}",
                Pos {
                    ord: 1,
                    total: None
                }
            ),
            "1",
        );
        assert_eq!(
            format!(
                "{}",
                Pos {
                    ord: 1,
                    total: Some(2)
                }
            ),
            "1/2",
        );

        Ok(())
    }

    #[test]
    fn pos_from_str() -> Result<(), anyhow::Error> {
        assert_eq!(
            Pos::from_str("1")?,
            Pos {
                ord: 1,
                total: None
            },
        );
        assert_eq!(
            Pos::from_str("1/2")?,
            Pos {
                ord: 1,
                total: Some(2)
            },
        );

        Ok(())
    }

    #[test]
    fn message() -> Result<(), anyhow::Error> {
        assert_eq!(
            format!("{}", Message::ErrorChromosomeNotFound),
            "ERROR_CHROMOSOME_NOT_FOUND",
        );
        assert_eq!(
            Message::from_str("ERROR_CHROMOSOME_NOT_FOUND")?,
            Message::ErrorChromosomeNotFound,
        );

        Ok(())
    }

    #[test]
    fn ann_field_display() {
        let value = AnnField {
            allele: Allele::Alt {
                alternative: String::from("A"),
            },
            consequences: vec![Consequence::MissenseVariant],
            putative_impact: PutativeImpact::Moderate,
            gene_symbol: String::from("GENE"),
            gene_id: String::from("HGNC:gene_id"),
            feature_type: FeatureType::SoTerm {
                term: SoFeature::Transcript,
            },
            feature_id: String::from("feature_id"),
            feature_biotype: FeatureBiotype::Coding,
            rank: Some(Rank { ord: 1, total: 2 }),
            hgvs_t: Some(String::from("HGVS.c")),
            hgvs_p: Some(String::from("HGVS.p")),
            tx_pos: Some(Pos {
                ord: 1,
                total: None,
            }),
            cds_pos: Some(Pos {
                ord: 1,
                total: Some(2),
            }),
            protein_pos: Some(Pos {
                ord: 1,
                total: None,
            }),
            distance: Some(1),
            messages: Some(vec![Message::ErrorChromosomeNotFound]),
        };

        assert_eq!(
            format!("{}", &value),
            "A|missense_variant|MODERATE|GENE|HGNC:gene_id|transcript|feature_id|Coding|1/2|HGVS.c\
            |HGVS.p|1|1/2|1|1|ERROR_CHROMOSOME_NOT_FOUND"
        );
    }

    #[test]
    fn ann_field_from_str() -> Result<(), anyhow::Error> {
        let value = "A|missense_variant|MODERATE|GENE|HGNC:gene_id|transcript|feature_id|\
        Coding|1/2|HGVS.c|HGVS.p|1|1/2|1|1|ERROR_CHROMOSOME_NOT_FOUND";

        let field = AnnField::from_str(value)?;
        assert_eq!(format!("{}", &field), value);

        Ok(())
    }
}
