//! Reference plugin: mehari's built-in `nmd` plugin, ported to the plugin interface in
//! `mehari/wit/plugin.wit`.

wit_bindgen::generate!({
    path: "../../mehari/wit",
    world: "plugin",
});

/// Premature stop codons starting within this many coding bases escape NMD by re-initiation.
const START_PROXIMAL_BASES: i32 = 100;
/// Premature stop codons starting at most this many bases before the last exon-exon junction
/// escape NMD.
const LAST_JUNCTION_BASES: i32 = 50;

struct Nmd;

impl Guest for Nmd {
    fn columns() -> Vec<String> {
        vec!["nmd".into()]
    }

    fn annotate(
        _variants: Vec<VcfVariant>,
        annotations: Vec<Annotation>,
        contexts: &Contexts,
    ) -> Vec<Vec<Option<String>>> {
        (0..)
            .zip(&annotations)
            .map(|(index, ann)| vec![predict(contexts, index, ann)])
            .collect()
    }
}

/// One value for every stop-gained or frameshift annotation, see `NmdPrediction` in mehari.
fn predict(contexts: &Contexts, index: u32, ann: &Annotation) -> Option<String> {
    let truncating = ann
        .consequences
        .iter()
        .any(|csq| csq == "stop_gained" || csq == "frameshift_variant");
    truncating.then(|| classify(contexts, index, ann).unwrap_or("unknown").into())
}

/// `None` if the premature stop codon (PTC) or the exon structure is unknown.
fn classify(contexts: &Contexts, index: u32, ann: &Annotation) -> Option<&'static str> {
    let transcript = contexts.transcript(index)?;
    let reference = contexts.reference(index)?;
    let alternative = contexts.alternative(index)?;
    let Some(stop) = alternative.protein.find('*') else {
        return Some("no_stop_codon");
    };

    // 1-based positions on the mutant transcript.
    let ptc = reference.cds_start + 3 * stop as i32;
    let last_exon_start = transcript.exons.last()?.tx_start;
    // Insertions and deletions before the last junction move it.
    let last_junction = if ann.cdna_pos? < last_exon_start {
        last_exon_start + alternative.transcript.len() as i32 - reference.transcript.len() as i32
    } else {
        last_exon_start
    };

    let ptc_c = ptc - reference.cds_start + 1;
    Some(
        if last_junction - ptc <= LAST_JUNCTION_BASES || ptc_c <= START_PROXIMAL_BASES {
            "escaping"
        } else {
            "triggering"
        },
    )
}

export!(Nmd);
