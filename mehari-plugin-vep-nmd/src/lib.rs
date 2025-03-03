use extism_convert::Json;
use extism_pdk::*;

use mehari_plugins::*;

#[plugin_fn]
pub fn process(Json(tva): Json<TranscriptVariationAllele>) -> FnResult<Option<Json<Annotation>>> {
    Ok(transcript_escapes_nmd(&tva).then_some(Json(Annotation::new("NMD".into()))))
}

#[plugin_fn]
pub fn feature_type() -> FnResult<String> {
    Ok("Transcript".into())
}

#[plugin_fn]
pub fn header_info() -> FnResult<Json<HeaderInfo>> {
    Ok(Json(HeaderInfo {
        tag: "NMD".into(),
        description: "Nonsense-mediated mRNA decay escaping variants prediction".into(),
    }))
}

// Core logic for NMD prediction
fn transcript_escapes_nmd(tva: &TranscriptVariationAllele) -> bool {
    // To qualify for NMD, at least one of the following consequences is required:
    //   "stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant",
    if !tva
        .overlap_consequences()
        .iter()
        .any(|oc| INCLUDE_SO.contains(&oc.so_term().as_ref()))
    {
        return false;
    }

    let transcript = tva.transcript();
    let exons = transcript.exons();
    let strand = transcript.strand();
    let variant_feature_end_position = tva.variation_feature().seq_region().end();

    // Rules for NMD prediction
    variant_within_last_exon(variant_feature_end_position, exons)
        || variant_upstream_of_penultimate_exon(variant_feature_end_position, strand, exons)
        || variant_within_first_100_coding_bases(tva.transcript_variation())
        || variant_is_intronless(transcript)
}

// Included SO terms
const INCLUDE_SO: [&str; 4] = [
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
];

/// Checks whether the variant location falls within the last exon of the transcript
fn variant_within_last_exon(variant_feature_end_position: u64, exons: &[Exon]) -> bool {
    let vf_end = variant_feature_end_position;
    exons
        .last()
        .map(|exon| vf_end >= exon.start() && vf_end <= exon.end())
        .unwrap_or(false)
}

fn variant_upstream_of_penultimate_exon(
    variant_feature_end_position: u64,
    strand: isize,
    exons: &[Exon],
) -> bool {
    let vf_end = variant_feature_end_position;
    exons
        .get(exons.len() - 2)
        .map(|second_last_exon| {
            let coding_region_end = if strand == -1 {
                second_last_exon.start().saturating_add(51)
            } else {
                second_last_exon.end().saturating_sub(51)
            };
            vf_end >= coding_region_end && vf_end <= second_last_exon.end()
        })
        .unwrap_or(false)
}

fn variant_within_first_100_coding_bases(tv: &TranscriptVariation) -> bool {
    tv.cds().end() <= 101
}

fn variant_is_intronless(transcript: &Transcript) -> bool {
    !transcript.has_introns()
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_has_nmd() {
        panic!()
    }
}
