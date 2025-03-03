/// Contains (placeholders for) structs and traits for the plugin interface.
use derive_builder::Builder;
use derive_new::new;
use getset::{CopyGetters, Getters};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, new)]
pub struct Annotation {
    pub annotation: String
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct HeaderInfo {
    pub tag: String,
    pub description: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Getters, Builder)]
#[getset(get = "pub")]
pub struct TranscriptVariationAllele {
    transcript: Transcript,
    transcript_variation: TranscriptVariation,
    variation_feature: VariationFeature,
    overlap_consequences: Vec<OverlapConsequence>,
}

#[derive(
    Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Getters, CopyGetters, Builder, Default,
)]
pub struct Transcript {
    #[getset(get = "pub")]
    introns: Vec<Intron>,
    #[getset(get = "pub")]
    exons: Vec<Exon>,
    #[getset(get_copy = "pub")]
    strand: isize,
}

impl Transcript {
    pub fn has_introns(&self) -> bool {
        !self.introns.is_empty()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Getters, new)]
#[getset(get = "pub")]
pub struct TranscriptVariation {
    cds: CodingSequence,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, CopyGetters, new)]
#[getset(get_copy = "pub")]
pub struct CodingSequence {
    start: u64,
    end: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, CopyGetters, new)]
#[getset(get_copy = "pub")]
pub struct Intron {
    start: u64,
    end: u64,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, CopyGetters, new)]
#[getset(get_copy = "pub")]
pub struct Exon {
    start: u64,
    end: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Getters, new, Default)]
#[getset(get = "pub")]
pub struct OverlapConsequence {
    so_term: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Getters, new)]
#[getset(get = "pub")]
pub struct VariationFeature {
    seq_region: SeqRegion,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, CopyGetters, new)]
#[getset(get_copy = "pub")]
pub struct SeqRegion {
    start: u64,
    end: u64,
}
