use hgvs::parser::HgvsVariant;
use itertools::Itertools;

use super::{TranscriptPlugin, TxAnnotation, TxContext};
use crate::annotate::seqvars::consequence::VcfVariant;
use crate::annotate::seqvars::consequence::terms::{ANN_NMD, AnnField, Consequence};
use crate::errors::SeqvarsError;
use crate::pbs::txs::Strand;

/// Premature stop codons starting within this many coding bases escape NMD by re-initiation.
const START_PROXIMAL_BASES: i32 = 100;
/// Premature stop codons starting at most this many bases before the last exon-exon junction
/// escape NMD.
const LAST_JUNCTION_BASES: i32 = 50;

/// The only value of VEP's NMD plugin.
const VEP_ESCAPING: &str = "NMD_escaping_variant";

/// Predicts whether the premature stop codon (PTC) of stop-gained and frameshift variants triggers
/// nonsense-mediated mRNA decay (NMD), after VEP's NMD plugin.
///
/// The rules apply to the PTC: for a frameshift, the first stop of the new reading frame. Every
/// stop-gained or frameshift annotation gets an [`NmdPrediction`]; splice variants get none as
/// their PTC is unknown.
///
/// With `vep_compatible`, the plugin mimics VEP instead: it applies VEP's rules and thresholds
/// to the variant position and reports `NMD_escaping_variant` or nothing.
#[derive(Debug, Default)]
pub struct NmdPlugin {
    pub vep_compatible: bool,
}

/// The `nmd` column value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, strum::Display)]
#[strum(serialize_all = "snake_case")]
pub enum NmdPrediction {
    /// The PTC lies in the last exon (which includes single-exon transcripts), or starts within
    /// the last 50 bases of the penultimate exon or the first 100 coding bases.
    Escaping,
    /// None of the escape rules applies.
    Triggering,
    /// The mutant transcript has no stop codon, so non-stop decay applies rather than NMD.
    NoStopCodon,
    /// The transcript model does not allow locating the PTC.
    Unknown,
}

impl TranscriptPlugin for NmdPlugin {
    fn columns(&self) -> Vec<String> {
        vec![ANN_NMD.into()]
    }

    fn annotate(
        &self,
        variants: &[VcfVariant],
        txs: &mut [TxAnnotation<'_>],
    ) -> Result<(), SeqvarsError> {
        for TxAnnotation { context, ann } in txs {
            let value = if self.vep_compatible {
                (vep_escapes(context, ann, variants) == Some(true)).then(|| VEP_ESCAPING.into())
            } else {
                let truncating = ann.consequences.iter().any(|csq| {
                    matches!(
                        csq,
                        Consequence::StopGained | Consequence::FrameshiftVariant
                    )
                });
                truncating.then(|| {
                    predict(context, ann)
                        .unwrap_or(NmdPrediction::Unknown)
                        .to_string()
                })
            };
            if let Some(value) = value {
                ann.custom_fields.insert(ANN_NMD.into(), Some(value));
            }
        }
        Ok(())
    }
}

/// `None` if the PTC or the exon structure is unknown.
fn predict(context: &TxContext<'_>, ann: &AnnField) -> Option<NmdPrediction> {
    let exons = &context.tx?.genome_alignments.first()?.exons;
    let ref_data = context.ref_data()?;
    let alt_data = context.alt_data()?;
    let Some(stop) = alt_data.aa_sequence.find('*') else {
        return Some(NmdPrediction::NoStopCodon);
    };

    // 1-based positions on the mutant transcript.
    let ptc = ref_data.cds_start + 3 * stop as i32;
    // Despite its name, `alt_cds_start_i` is the exon's 1-based transcript start.
    let last_exon_start = exons.iter().filter_map(|exon| exon.alt_cds_start_i).max()?;
    // Insertions and deletions before the last junction move it.
    let last_junction = if ann.cdna_pos.as_ref()?.ord < last_exon_start {
        last_exon_start + alt_data.transcript_sequence.len() as i32
            - ref_data.transcript_sequence.len() as i32
    } else {
        last_exon_start
    };

    Some(classify(ptc - ref_data.cds_start + 1, last_junction - ptc))
}

/// The rules for a PTC starting at `ptc_c` (c.) and `before_junction` bases before the last
/// exon-exon junction (zero or less: in the last exon).
fn classify(ptc_c: i32, before_junction: i32) -> NmdPrediction {
    if before_junction <= LAST_JUNCTION_BASES || ptc_c <= START_PROXIMAL_BASES {
        NmdPrediction::Escaping
    } else {
        NmdPrediction::Triggering
    }
}

/// VEP's NMD plugin, which applies its rules to the variant: `Some(true)` if it escapes NMD.
fn vep_escapes(context: &TxContext<'_>, ann: &AnnField, variants: &[VcfVariant]) -> Option<bool> {
    // VEP neither knows compound variants nor annotates variants without CDS coordinates.
    let [variant] = variants else {
        return None;
    };
    let triggers = [
        Consequence::StopGained,
        Consequence::FrameshiftVariant,
        Consequence::SpliceDonorVariant,
        Consequence::SpliceAcceptorVariant,
    ];
    if !ann.consequences.iter().any(|csq| triggers.contains(csq)) || ann.cds_pos.is_none() {
        return None;
    }
    let HgvsVariant::CdsVariant { loc_edit, .. } = context.var_c.as_ref()? else {
        return None;
    };
    let loc = loc_edit.loc.inner();
    // VEP ends an insertion at the base before it.
    let cds_end = if loc_edit.edit.inner().is_ins() {
        loc.start.base
    } else {
        loc.end.base
    };

    let alignment = context.tx?.genome_alignments.first()?;
    let exons = alignment
        .exons
        .iter()
        .sorted_by_key(|exon| exon.alt_cds_start_i)
        .collect_vec();
    // 1-based genomic end; trimming the VCF anchor base only moves the start.
    let end = variant.position + variant.reference.len() as i32 - 1;
    let last = exons.last()?;
    let in_last_exon = last.alt_start_i < end && end <= last.alt_end_i;
    let near_last_junction = exons.len() > 1 && {
        let penultimate = exons[exons.len() - 2];
        if alignment.strand == Strand::Minus as i32 {
            (penultimate.alt_start_i + 1..=penultimate.alt_start_i + 52).contains(&end)
        } else {
            (penultimate.alt_end_i - 51..=penultimate.alt_end_i).contains(&end)
        }
    };

    Some(cds_end <= 101 || in_last_exon || near_last_junction || exons.len() == 1)
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;
    use std::sync::Arc;

    use super::*;
    use crate::annotate::seqvars::consequence::logic::ConsequencePredictor;
    use crate::annotate::seqvars::consequence::{ConfigBuilder, load_tx_db};
    use crate::annotate::seqvars::provider::Provider as MehariProvider;

    #[rstest::rstest]
    #[case(100, 1000, NmdPrediction::Escaping)] // last start-proximal base
    #[case(101, 1000, NmdPrediction::Triggering)]
    #[case(500, 50, NmdPrediction::Escaping)]
    #[case(500, 51, NmdPrediction::Triggering)]
    #[case(500, 0, NmdPrediction::Escaping)] // first base of the last exon
    fn classify_rules(
        #[case] ptc_c: i32,
        #[case] before_junction: i32,
        #[case] expected: NmdPrediction,
    ) {
        assert_eq!(classify(ptc_c, before_junction), expected);
    }

    /// BRCA1 NM_007294.4: the CDS starts at n.114, the last exon at n.5581 (c.5468).
    #[rstest::rstest]
    #[case(
        "17:41276049:A:T",
        Consequence::StopGained,
        Some("escaping"),
        Some(VEP_ESCAPING)
    )] // c.65
    #[case("17:41267765:T:A", Consequence::StopGained, Some("triggering"), None)] // c.112
    #[case("17:41244528:G:T", Consequence::StopGained, Some("triggering"), None)] // c.3020
    #[case(
        "17:41199696:G:A",
        Consequence::StopGained,
        Some("escaping"),
        Some(VEP_ESCAPING)
    )] // c.5431
    // c.5411del lies 57 bases before the last junction, its stop in the last exon.
    #[case(
        "17:41199715:GA:G",
        Consequence::FrameshiftVariant,
        Some("escaping"),
        None
    )]
    #[case(
        "17:41197811:G:A",
        Consequence::StopGained,
        Some("escaping"),
        Some(VEP_ESCAPING)
    )] // c.5476
    #[case("17:41215920:G:T", Consequence::MissenseVariant, None, None)]
    fn brca1(
        #[case] spdi: &str,
        #[case] consequence: Consequence,
        #[case] expected: Option<&str>,
        #[case] expected_vep: Option<&str>,
    ) -> Result<(), anyhow::Error> {
        let [chromosome, position, reference, alternative] =
            spdi.split(':').collect::<Vec<_>>()[..]
        else {
            unreachable!()
        };
        let variant = VcfVariant {
            chromosome: chromosome.into(),
            position: position.parse()?,
            reference: reference.into(),
            alternative: alternative.into(),
        };

        for (vep_compatible, expected) in [(false, expected), (true, expected_vep)] {
            let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
            let provider = Arc::new(MehariProvider::new(
                tx_db,
                None::<PathBuf>,
                true,
                Default::default(),
            ));
            let predictor = ConsequencePredictor::new(
                provider,
                ConfigBuilder::default()
                    .plugins(vec![Arc::new(NmdPlugin { vep_compatible })])
                    .build()?,
            );

            let anns = predictor.predict(&variant)?.unwrap();
            let ann = anns
                .iter()
                .find(|ann| ann.feature_id == "NM_007294.4")
                .unwrap();

            assert!(ann.consequences.contains(&consequence), "{ann:?}");
            assert_eq!(
                ann.custom_fields.get(ANN_NMD).cloned().flatten().as_deref(),
                expected,
                "vep_compatible = {vep_compatible}"
            );
        }

        Ok(())
    }
}
