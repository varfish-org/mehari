//! Reference plugin: mehari's built-in `nmd` plugin, ported to the plugin interface in
//! `mehari/wit/plugin.wit`.

wit_bindgen::generate!({
    path: "../../mehari/wit",
    world: "plugin",
});

use exports::mehari::plugin::transcript_plugin::Guest;
use mehari::plugin::types::{Annotation, Prediction, Setup};

/// Premature stop codons starting within this many coding bases escape NMD by re-initiation.
const START_PROXIMAL_BASES: i32 = 100;
/// Premature stop codons starting at most this many bases before the last exon-exon junction
/// escape NMD.
const LAST_JUNCTION_BASES: i32 = 50;

struct Nmd;

impl Guest for Nmd {
    fn init(setup: &Setup) -> Result<(), String> {
        if let Some((key, _)) = setup.options().first() {
            return Err(format!("unknown option {key:?}"));
        }
        setup.add_string_column("nmd");
        Ok(())
    }

    fn annotate(prediction: &Prediction) -> Result<(), String> {
        for annotation in prediction.annotations() {
            if let Some(value) = predict(&annotation) {
                annotation.set_string("nmd", value);
            }
        }
        Ok(())
    }
}

/// One value for every stop-gained or frameshift annotation, see `NmdPrediction` in mehari.
fn predict(annotation: &Annotation) -> Option<&'static str> {
    let truncating = annotation
        .consequences()
        .iter()
        .any(|csq| csq == "stop_gained" || csq == "frameshift_variant");
    truncating.then(|| classify(annotation).unwrap_or("unknown"))
}

/// `None` if the premature stop codon (PTC) or the exon structure is unknown.
fn classify(annotation: &Annotation) -> Option<&'static str> {
    let exons = annotation.exons();
    let cds_start = annotation.reference_cds_start()?;
    let reference_len = annotation.reference_transcript()?.len() as i32;
    let alternative_len = annotation.alternative_transcript()?.len() as i32;
    let Some(stop) = annotation.alternative_protein()?.find('*') else {
        return Some("no_stop_codon");
    };

    // 1-based positions on the mutant transcript.
    let ptc = cds_start + 3 * stop as i32;
    let last_exon_start = exons.last()?.tx_start;
    // Insertions and deletions before the last junction move it.
    let last_junction = if annotation.cdna_pos()?.ord < last_exon_start {
        last_exon_start + alternative_len - reference_len
    } else {
        last_exon_start
    };

    let ptc_c = ptc - cds_start + 1;
    Some(
        if last_junction - ptc <= LAST_JUNCTION_BASES || ptc_c <= START_PROXIMAL_BASES {
            "escaping"
        } else {
            "triggering"
        },
    )
}

export!(Nmd);
