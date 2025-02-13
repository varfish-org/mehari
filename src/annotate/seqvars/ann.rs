//! Code for annotating variants based on molecular consequence.
use enumflags2::bitflags;
use nom::Parser;
use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{alphanumeric1, digit1},
    combinator::{all_consuming, map},
    IResult,
};
use parse_display::{Display, FromStr};
use std::str::FromStr;
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
    utoipa::ToSchema,
)]
#[display(style = "UPPERCASE")]
#[serde(rename_all = "snake_case")]
pub enum PutativeImpact {
    High,
    Moderate,
    Low,
    Modifier,
}

/// Putative impact.
#[bitflags]
#[repr(u64)]
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
    utoipa::ToSchema,
)]
#[display(style = "snake_case")]
#[serde(rename_all = "snake_case")]
pub enum Consequence {
    // high impact
    /// "A feature ablation whereby the deleted region includes a transcript feature."
    /// SO:transcript_ablation, VEP:transcript_ablation
    TranscriptAblation,

    /// "A sequence variant whereby an exon is lost from the transcript."
    /// SO:exon_loss_variant, VEP:transcript_ablation
    ExonLossVariant,

    /// "A splice variant that changes the 2 base region at the 3' end of an intron."
    /// SO:splice_acceptor_variant, VEP:splice_acceptor_variant
    SpliceAcceptorVariant,

    /// "A splice variant that changes the 2 base region at the 5' end of an intron."
    /// SO:splice_donor_variant, VEP:splice_donor_variant
    SpliceDonorVariant,

    /// "A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript."
    /// SO:stop_gained, VEP:stop_gained
    StopGained,

    /// "A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three."
    /// SO:frameshift_variant, VEP:frameshift_variant
    FrameshiftVariant,

    /// "A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript."
    /// SO:stop_lost, VEP:stop_lost
    StopLost,

    /// "A codon variant that changes at least one base of the canonical start codon."
    /// SO:start_lost, VEP:start_lost
    StartLost,

    /// "A feature amplification of a region containing a transcript."
    /// SO:transcript_amplification, VEP:transcript_amplification
    TranscriptAmplification,

    // Currently never written out (because hgvs::parser::ProteinEdit::Ext not produced)
    /// "A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence."
    /// SO:feature_elongation, VEP:feature_elongation
    FeatureElongation,

    /// "A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence."
    /// SO:feature_truncation, VEP:feature_truncation
    FeatureTruncation,

    // moderate impact
    /// "An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon."
    /// SO:disruptive_inframe_insertion, VEP:inframe_insertion
    DisruptiveInframeInsertion,

    /// "An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon."
    /// SO:disruptive_inframe_deletion, VEP:inframe_deletion
    DisruptiveInframeDeletion,

    /// "An inframe increase in cds length that inserts one or more codons into the coding sequence between existing codons."
    /// SO:conservative_inframe_insertion, VEP:inframe_insertion
    ConservativeInframeInsertion,

    /// "An inframe decrease in cds length that deletes one or more entire codons from the coding sequence but does not change any remaining codons."
    /// SO:conservative_inframe_deletion, VEP:inframe_deletion
    ConservativeInframeDeletion,

    /// "A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved."
    /// SO:missense_variant, VEP:missense_variant
    MissenseVariant,

    /// "A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid."
    /// SO:rare_amino_acid_variant
    RareAminoAcidVariant,

    // Not used by mehari, but by VEP (we're usually more specific)
    // /// "A sequence_variant which is predicted to change the protein encoded in the coding sequence."
    // /// SO:protein_altering_variant, VEP:missense_variant
    // ProteinAlteringVariant,

    // low impact
    /// "A sequence variant that causes a change at the 5th base pair after the start of the intron in the orientation of the transcript."
    /// SO:splice_donor_5th_base_variant, VEP:splice_donor_5th_base_variant
    #[display("splice_donor_5th_base_variant")]
    #[serde(rename = "splice_donor_5th_base_variant")]
    SpliceDonorFifthBaseVariant,

    /// "A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron."
    /// SO:splice_region_variant, VEP:splice_region_variant
    SpliceRegionVariant,

    /// "A sequence variant that falls in the region between the 3rd and 6th base after splice junction (5' end of intron)."
    /// SO:splice_donor_region_variant, VEP:splice_donor_region_variant
    SpliceDonorRegionVariant,

    /// "A sequence variant that falls in the polypyrimidine tract at 3' end of intron between 17 and 3 bases from the end (acceptor -3 to acceptor -17)."
    /// SO:splice_polypyrimidine_tract_variant, VEP:splice_polypyrimidine_tract_variant
    SplicePolypyrimidineTractVariant,

    // Not used by mehari, but by VEP
    // /// "A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed."
    // /// SO:incomplete_terminal_codon_variant, VEP:incomplete_terminal_codon_variant
    // IncompleteTerminalCodonVariant
    /// "A sequence variant where at least one base in the start codon is changed, but the start remains."
    /// SO:start_retained_variant, VEP:start_retained_variant
    StartRetainedVariant,

    /// "A sequence variant where at least one base in the terminator codon is changed, but the terminator remains."
    /// SO:stop_retained_variant, VEP:stop_retained_variant
    StopRetainedVariant,

    /// "A sequence variant where there is no resulting change to the encoded amino acid."
    /// SO:synonymous_variant, VEP:synonymous_variant
    SynonymousVariant,

    // modifier
    /// "A sequence variant that changes the coding sequence."
    /// SO:coding_sequence_variant, VEP:coding_sequence_variant
    CodingSequenceVariant,

    // Not yet implemented
    /// "A transcript variant located with the sequence of the mature miRNA."
    /// SO:mature_miRNA_variant, VEP:mature_miRNA_variant
    #[display("mature_miRNA_variant")]
    #[serde(rename = "mature_miRNA_variant")]
    MatureMirnaVariant,

    /// "A UTR variant of exonic sequence of the 5' UTR."
    /// SO:5_prime_UTR_exon_variant, VEP:5_prime_UTR_variant
    #[display("5_prime_UTR_exon_variant")]
    #[serde(rename = "5_prime_UTR_exon_variant")]
    FivePrimeUtrExonVariant,

    /// "A UTR variant of intronic sequence of the 5' UTR."
    /// SO:5_prime_UTR_intron_variant, VEP:5_prime_UTR_variant
    #[display("5_prime_UTR_intron_variant")]
    #[serde(rename = "5_prime_UTR_intron_variant")]
    FivePrimeUtrIntronVariant,

    /// "A UTR variant of exonic sequence of the 3' UTR."
    /// SO:3_prime_UTR_exon_variant, VEP:3_prime_UTR_variant
    #[display("3_prime_UTR_exon_variant")]
    #[serde(rename = "3_prime_UTR_exon_variant")]
    ThreePrimeUtrExonVariant,

    /// "A UTR variant of intronic sequence of the 3' UTR."
    /// SO:3_prime_UTR_intron_variant, VEP:3_prime_UTR_variant
    #[display("3_prime_UTR_intron_variant")]
    #[serde(rename = "3_prime_UTR_intron_variant")]
    ThreePrimeUtrIntronVariant,

    /// "A sequence variant that changes non-coding exon sequence in a non-coding transcript."
    /// SO:non_coding_transcript_exon_variant, VEP:non_coding_transcript_variant
    NonCodingTranscriptExonVariant,

    /// "A sequence variant that changes non-coding intron sequence in a non-coding transcript."
    /// SO:non_coding_transcript_intron_variant, VEP:non_coding_transcript_variant
    NonCodingTranscriptIntronVariant,

    // Not used by mehari, but by VEP
    // /// "A transcript variant of a protein coding gene."
    // /// SO:coding_transcript_variant, VEP:coding_transcript_variant
    // CodingTranscriptVariant,
    /// "A sequence variant located 5' of a gene."
    /// SO:upstream_gene_variant, VEP:upstream_gene_variant
    UpstreamGeneVariant,

    /// "A sequence variant located 3' of a gene."
    /// SO:downstream_gene_variant, VEP:downstream_gene_variant
    DownstreamGeneVariant,

    /// "A feature ablation whereby the deleted region includes a transcription factor binding site."
    /// SO:TFBS_ablation, VEP:TFBS_ablation
    #[display("TFBS_ablation")]
    #[serde(rename = "TFBS_ablation")]
    TfbsAblation,

    /// "A feature amplification of a region containing a transcription factor binding site."
    /// SO:TFBS_amplification, VEP:TFBS_amplification
    #[display("TFBS_amplification")]
    #[serde(rename = "TFBS_amplification")]
    TfbsAmplification,

    /// "A sequence variant located within a transcription factor binding site."
    /// SO:TF_binding_site_variant, VEP:TF_binding_site_variant
    #[display("TF_binding_site_variant")]
    #[serde(rename = "TF_binding_site_variant")]
    TfBindingSiteVariant,

    /// "A feature ablation whereby the deleted region includes a regulatory region."
    /// SO:regulatory_region_ablation, VEP:regulatory_region_ablation
    RegulatoryRegionAblation,

    /// "A feature amplification of a region containing a regulatory region."
    /// SO:regulatory_region_amplification, VEP:regulatory_region_amplification
    RegulatoryRegionAmplification,

    /// "A sequence variant located within a regulatory region."
    /// SO:regulatory_region_variant, VEP:regulatory_region_variant
    RegulatoryRegionVariant,

    /// "A sequence variant located in the intergenic region, between genes."
    /// SO:intergenic_variant, VEP:intergenic_variant
    IntergenicVariant,

    // Not used by mehari, but by VEP
    // /// "A sequence_variant is a non exact copy of a sequence_feature or genome exhibiting one or more sequence_alteration."
    // /// SO:sequence_variant, VEP:sequence_variant
    // SequenceVariant,
    /// "A transcript variant occurring within an intron."
    /// SO:intron_variant, VEP:intron_variant
    IntronVariant,

    /// "A sequence variant where the structure of the gene is changed."
    /// SO:gene_variant
    GeneVariant,
}

impl From<Consequence> for PutativeImpact {
    fn from(val: Consequence) -> Self {
        use Consequence::*;
        match val {
            TranscriptAblation
            | ExonLossVariant
            | SpliceAcceptorVariant
            | SpliceDonorVariant
            | StopGained
            | FrameshiftVariant
            | StopLost
            | StartLost
            | TranscriptAmplification
            | FeatureElongation
            | FeatureTruncation => PutativeImpact::High,
            DisruptiveInframeInsertion
            | DisruptiveInframeDeletion
            | ConservativeInframeInsertion
            | ConservativeInframeDeletion
            | MissenseVariant
            | RareAminoAcidVariant => PutativeImpact::Moderate,
            SpliceDonorFifthBaseVariant
            | SpliceRegionVariant
            | SpliceDonorRegionVariant
            | SplicePolypyrimidineTractVariant
            | StartRetainedVariant
            | StopRetainedVariant
            | SynonymousVariant => PutativeImpact::Low,
            CodingSequenceVariant
            | MatureMirnaVariant
            | FivePrimeUtrExonVariant
            | FivePrimeUtrIntronVariant
            | ThreePrimeUtrExonVariant
            | ThreePrimeUtrIntronVariant
            | NonCodingTranscriptExonVariant
            | NonCodingTranscriptIntronVariant
            | UpstreamGeneVariant
            | DownstreamGeneVariant
            | TfbsAblation
            | TfbsAmplification
            | TfBindingSiteVariant
            | RegulatoryRegionAblation
            | RegulatoryRegionAmplification
            | RegulatoryRegionVariant
            | IntergenicVariant
            | IntronVariant
            | GeneVariant => PutativeImpact::Modifier,
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
    use nom::Parser;

    pub static NA_IUPAC: &str = "ACGTURYMKWSBDHVNacgturymkwsbdhvn";

    pub fn na1(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        take_while1(|c: char| NA_IUPAC.contains(c)).parse(input)
    }
}

impl Allele {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        all_consuming(alt((
            Self::parse_compound,
            Self::parse_alt_ref,
            Self::parse_alt,
        )))
        .parse(input)
    }

    fn parse_compound(input: &str) -> IResult<&str, Self> {
        map(
            (
                parse::na1,
                tag("-"),
                alphanumeric1,
                tag(":"),
                digit1,
                tag("_"),
                parse::na1,
                tag(">"),
                parse::na1,
            ),
            |(alternative, _, other_chrom, _, other_pos, _, other_ref, _, other_alt)| {
                Allele::Compound {
                    alternative: alternative.to_string(),
                    other_chrom: other_chrom.to_string(),
                    other_pos: other_pos.parse::<u32>().unwrap(),
                    other_ref: other_ref.to_string(),
                    other_alt: other_alt.to_string(),
                }
            },
        )
        .parse(input)
    }

    fn parse_alt_ref(input: &str) -> IResult<&str, Self> {
        map(
            (parse::na1, tag("-"), parse::na1),
            |(alternative, _, reference)| Allele::AltRef {
                alternative: alternative.to_string(),
                reference: reference.to_string(),
            },
        )
        .parse(input)
    }

    fn parse_alt(input: &str) -> IResult<&str, Self> {
        map(parse::na1, |alternative| Allele::Alt {
            alternative: alternative.to_string(),
        })
        .parse(input)
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
    utoipa::ToSchema,
)]
#[display(style = "snake_case")]
pub enum SoFeature {
    Transcript,
}

/// Enum for `AnnField::feature_type`.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Display,
    serde::Deserialize,
    serde::Serialize,
    utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
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
    utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum FeatureBiotype {
    /// Is coding transcript.
    Coding,
    /// Is non-coding transcript.
    Noncoding,
    /// Is in MANE Select set.
    ManeSelect,
    /// Is in MANE Plus Clinical set.
    ManePlusClinical,
}

impl FeatureBiotype {
    pub fn is_coding(&self) -> bool {
        matches!(self, FeatureBiotype::Coding)
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
    utoipa::ToSchema,
)]
#[display("{ord}/{total}")]
pub struct Rank {
    pub ord: i32,
    pub total: i32,
}

impl Rank {
    #[inline]
    pub fn is_first(&self) -> bool {
        self.ord == 1
    }

    #[inline]
    pub fn is_last(&self) -> bool {
        self.ord == self.total
    }
}

/// Position, optionally with total length.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Default,
    serde::Deserialize,
    serde::Serialize,
    utoipa::ToSchema,
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
        map((tag("-"), digit1::<&str, _>), |(sign, num)| {
            let num = num.parse::<i32>().unwrap();
            if sign == "-" {
                -num
            } else {
                num
            }
        })
        .parse(input)
    }

    fn parse_number_nosign(input: &str) -> IResult<&str, i32> {
        map(digit1::<&str, _>, |num| num.parse::<i32>().unwrap()).parse(input)
    }

    fn parse_number(input: &str) -> IResult<&str, i32> {
        alt((Self::parse_number_neg, Self::parse_number_nosign)).parse(input)
    }

    fn parse_with_total(input: &str) -> IResult<&str, Self> {
        map((Self::parse_number, tag("/"), digit1), |(num, _, total)| {
            Pos {
                ord: num,
                total: Some(total.parse::<i32>().unwrap()),
            }
        })
        .parse(input)
    }

    fn parse_no_total(input: &str) -> IResult<&str, Self> {
        map(Self::parse_number, |num| Pos {
            ord: num,
            total: None,
        })
        .parse(input)
    }
}

impl FromStr for Pos {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        all_consuming(alt((Self::parse_with_total, Self::parse_no_total)))
            .parse(s)
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
    utoipa::ToSchema,
)]
#[display(style = "SNAKE_CASE")]
#[serde(rename_all = "snake_case")]
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
    pub feature_biotype: Vec<FeatureBiotype>,
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
    /// Strand of the alignment
    pub strand: i32,
    /// Optional list of warnings and error messages.
    pub messages: Option<Vec<Message>>,
}

impl Default for AnnField {
    fn default() -> Self {
        Self {
            allele: Allele::Alt {
                alternative: Default::default(),
            },
            consequences: vec![],
            putative_impact: PutativeImpact::Modifier,
            gene_symbol: Default::default(),
            gene_id: Default::default(),
            feature_type: FeatureType::SoTerm {
                term: SoFeature::Transcript,
            },
            feature_id: Default::default(),
            feature_biotype: vec![FeatureBiotype::Coding],
            rank: Default::default(),
            hgvs_t: Default::default(),
            hgvs_p: Default::default(),
            tx_pos: Default::default(),
            cds_pos: Default::default(),
            protein_pos: Default::default(),
            distance: Default::default(),
            strand: Default::default(),
            messages: Default::default(),
        }
    }
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
        let feature_biotype = fields
            .next()
            .unwrap()
            .split('&')
            .map(|s| s.parse())
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| anyhow::anyhow!("could not parse feature biotype: {}", e))?;
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
        let strand = fields.next().unwrap();
        let strand = strand.parse()?;
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
            strand,
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
        write!(
            f,
            "{}",
            self.feature_biotype
                .iter()
                .map(|t| format!("{}", t))
                .collect::<Vec<_>>()
                .join("&")
        )?;
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
        write!(f, "{}", self.strand)?;
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
        assert_eq!(format!("{}", Consequence::TfbsAblation), "TFBS_ablation");
        assert_eq!(
            format!("{}", Consequence::ThreePrimeUtrExonVariant),
            "3_prime_UTR_exon_variant"
        );
        assert_eq!(
            format!("{}", Consequence::FivePrimeUtrIntronVariant),
            "5_prime_UTR_intron_variant"
        );
        assert_eq!(
            format!("{}", Consequence::MatureMirnaVariant),
            "mature_miRNA_variant"
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
            Consequence::from_str("TFBS_ablation")?,
            Consequence::TfbsAblation,
        );
        assert_eq!(
            Consequence::from_str("3_prime_UTR_exon_variant")?,
            Consequence::ThreePrimeUtrExonVariant,
        );
        assert_eq!(
            Consequence::from_str("5_prime_UTR_intron_variant")?,
            Consequence::FivePrimeUtrIntronVariant,
        );
        assert_eq!(
            Consequence::from_str("mature_miRNA_variant")?,
            Consequence::MatureMirnaVariant,
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
            feature_biotype: vec![FeatureBiotype::Coding],
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
            strand: 0,
            messages: Some(vec![Message::ErrorChromosomeNotFound]),
        };

        assert_eq!(
            format!("{}", &value),
            "A|missense_variant|MODERATE|GENE|HGNC:gene_id|transcript|feature_id|Coding|1/2|HGVS.c\
            |HGVS.p|1|1/2|1|1|0|ERROR_CHROMOSOME_NOT_FOUND"
        );
    }

    #[test]
    fn ann_field_from_str() -> Result<(), anyhow::Error> {
        let value = "A|missense_variant|MODERATE|GENE|HGNC:gene_id|transcript|feature_id|\
        Coding|1/2|HGVS.c|HGVS.p|1|1/2|1|1|0|ERROR_CHROMOSOME_NOT_FOUND";

        let field = AnnField::from_str(value)?;
        assert_eq!(format!("{}", &field), value);

        Ok(())
    }
}
