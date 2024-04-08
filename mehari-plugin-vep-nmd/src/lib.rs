use extism_convert::{Json, ToBytes};
use extism_pdk::*;
use serde::Serialize;

use mehari_plugins::*;

#[derive(ToBytes, Serialize)]
#[encoding(Json)]
struct HeaderInfo {
    tag: String,
    description: String,
}

#[plugin_fn]
pub fn process(record: String) -> FnResult<String> {
    let mut record: Dummy = serde_json::from_str(&record).expect("Failed deserializing record");
    record.change_something("after");
    let record = serde_json::to_string(&record).expect("Failed serializing record");
    Ok(record)
}

#[plugin_fn]
pub fn feature_type() -> FnResult<String> {
    Ok("Transcript".into())
}

#[plugin_fn]
pub fn header_info() -> FnResult<HeaderInfo> {
    Ok(HeaderInfo {
        tag: "NMD".into(),
        description: "Nonsense-mediated mRNA decay escaping variants prediction".into(),
    })
}

// Included SO terms
const INCLUDE_SO: &[&str] = &[
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
];

struct TranscriptVariationAllele {
    transcript: Transcript,
    transcript_variation: TranscriptVariation,
    variation_feature: VariationFeature,
}

impl TranscriptVariationAllele {
    fn overlap_consequences(&self) -> &Vec<OverlapConsequence> {
        todo!()
    }
}

struct Transcript;

impl Transcript {
    fn introns(&self) -> &Vec<Intron> {
        todo!()
    }

    fn exons(&self) -> &Vec<Exon> {
        todo!()
    }

    fn has_introns(&self) -> bool {
        todo!()
    }

    fn strand(&self) -> isize {
        todo!()
    }
}

struct TranscriptVariation;

impl TranscriptVariation {
    fn cds_start(&self) -> Option<u64> {
        todo!()
    }
    fn cds_end(&self) -> Option<u64> {
        todo!()
    }
}

struct Intron;
struct Exon;

impl Exon {
    fn start(&self) -> u64 {
        todo!()
    }

    fn end(&self) -> u64 {
        todo!()
    }
}
struct OverlapConsequence;

impl OverlapConsequence {
    fn so_term(&self) -> &'static str {
        todo!()
    }
}

struct VariationFeature;

impl VariationFeature {
    fn seq_region_end(&self) -> u64 {
        todo!()
    }
}

// Core logic for NMD prediction
fn run(tva: &TranscriptVariationAllele) -> Option<&'static str> {
    // Condition for inclusion
    if !tva
        .overlap_consequences()
        .iter()
        .any(|oc| INCLUDE_SO.contains(&oc.so_term()))
    {
        return None;
    }

    // Access transcript and variation information
    let tr = &tva.transcript;
    let tv = &tva.transcript_variation;

    // Rules for NMD prediction
    if let Some(variant_coding_region) = tv.cds_end() {
        if variant_coding_region <= 101 || variant_exon_check(tva) || !tr.has_introns() {
            return Some("NMD_escaping_variant");
        }
    }
    None
}

// Check variant location within exons
fn variant_exon_check(tva: &TranscriptVariationAllele) -> bool {
    let tr = &tva.transcript;
    let vf_end = tva.variation_feature.seq_region_end();
    let exons = tr.exons();

    // Check last exon
    let last_exon = exons.last().unwrap();
    if vf_end >= last_exon.start() && vf_end <= last_exon.end() {
        return true;
    }

    // Check second-to-last exon
    if let Some(second_last_exon) = exons.get(exons.len() - 2) {
        let coding_region_end = if tr.strand() == -1 {
            second_last_exon.start().saturating_add(51)
        } else {
            second_last_exon.end().saturating_sub(51)
        };
        if vf_end >= coding_region_end && vf_end <= second_last_exon.end() {
            return true;
        }
    }

    false
}
