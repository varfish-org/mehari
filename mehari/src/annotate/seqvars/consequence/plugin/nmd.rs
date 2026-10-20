use super::{TranscriptPlugin, TxAnnotation, TxContext};
use crate::annotate::seqvars::consequence::VcfVariant;
use crate::annotate::seqvars::consequence::terms::{ANN_NMD, AnnField, Consequence};
use crate::errors::SeqvarsError;

/// Premature stop codons starting within this many coding bases escape NMD by re-initiation.
const START_PROXIMAL_BASES: i32 = 100;
/// Premature stop codons starting at most this many bases before the last exon-exon junction
/// escape NMD.
const LAST_JUNCTION_BASES: i32 = 50;

/// Flags stop-gained and frameshift variants whose premature stop codon (PTC) escapes
/// nonsense-mediated mRNA decay (NMD), with the rules of VEP's NMD plugin.
///
/// The PTC escapes if it starts within the first 100 coding bases, within the last 50 bases of
/// the penultimate exon, or in the last exon (which covers single-exon transcripts). Unlike VEP,
/// the rules apply to the PTC rather than the variant: for a frameshift, that is the first stop
/// of the new reading frame. Splice variants are not flagged as their PTC is unknown.
#[derive(Debug, Default)]
pub struct NmdPlugin;

impl TranscriptPlugin for NmdPlugin {
    fn columns(&self) -> Vec<String> {
        vec![ANN_NMD.into()]
    }

    fn annotate(
        &self,
        _variants: &[VcfVariant],
        txs: &mut [TxAnnotation<'_>],
    ) -> Result<(), SeqvarsError> {
        for TxAnnotation { context, ann } in txs {
            let truncating = ann.consequences.iter().any(|csq| {
                matches!(
                    csq,
                    Consequence::StopGained | Consequence::FrameshiftVariant
                )
            });
            if truncating && escapes_nmd(context, ann) == Some(true) {
                ann.custom_fields
                    .insert(ANN_NMD.into(), Some("NMD_escaping_variant".into()));
            }
        }
        Ok(())
    }
}

/// `None` if the PTC or the exon structure is unknown.
fn escapes_nmd(context: &TxContext<'_>, ann: &AnnField) -> Option<bool> {
    let exons = &context.tx?.genome_alignments.first()?.exons;
    let ref_data = context.ref_data()?;
    let alt_data = context.alt_data()?;

    // 1-based positions on the mutant transcript.
    let ptc = ref_data.cds_start + 3 * alt_data.aa_sequence.find('*')? as i32;
    // Despite its name, `alt_cds_start_i` is the exon's 1-based transcript start.
    let last_exon_start = exons.iter().filter_map(|exon| exon.alt_cds_start_i).max()?;
    // Insertions and deletions before the last junction move it.
    let last_junction = if ann.cdna_pos.as_ref()?.ord < last_exon_start {
        last_exon_start + alt_data.transcript_sequence.len() as i32
            - ref_data.transcript_sequence.len() as i32
    } else {
        last_exon_start
    };

    Some(escapes(ptc - ref_data.cds_start + 1, last_junction - ptc))
}

/// The rules for a PTC starting at `ptc_c` (c.) and `before_junction` bases before the last
/// exon-exon junction (zero or less: in the last exon).
fn escapes(ptc_c: i32, before_junction: i32) -> bool {
    ptc_c <= START_PROXIMAL_BASES || before_junction <= LAST_JUNCTION_BASES
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
    #[case(100, 1000, true)] // last start-proximal base
    #[case(101, 1000, false)]
    #[case(500, 50, true)] // first of the last 50 bases of the penultimate exon
    #[case(500, 51, false)]
    #[case(500, 0, true)] // first base of the last exon
    fn escape_rules(#[case] ptc_c: i32, #[case] before_junction: i32, #[case] expected: bool) {
        assert_eq!(escapes(ptc_c, before_junction), expected);
    }

    /// BRCA1 NM_007294.4: the CDS starts at n.114, the last exon at n.5581 (c.5468).
    #[rstest::rstest]
    #[case(
        "17:41276049:A:T",
        Consequence::StopGained,
        Some("NMD_escaping_variant")
    )] // c.65
    #[case("17:41267765:T:A", Consequence::StopGained, None)] // c.112
    #[case("17:41244528:G:T", Consequence::StopGained, None)] // c.3020
    #[case(
        "17:41199696:G:A",
        Consequence::StopGained,
        Some("NMD_escaping_variant")
    )] // c.5431
    // c.5411del lies 57 bases before the last junction, its stop in the last exon.
    #[case(
        "17:41199715:GA:G",
        Consequence::FrameshiftVariant,
        Some("NMD_escaping_variant")
    )]
    #[case(
        "17:41197811:G:A",
        Consequence::StopGained,
        Some("NMD_escaping_variant")
    )] // c.5476
    #[case("17:41215920:G:T", Consequence::MissenseVariant, None)]
    fn brca1(
        #[case] spdi: &str,
        #[case] consequence: Consequence,
        #[case] expected: Option<&str>,
    ) -> Result<(), anyhow::Error> {
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
                .plugins(vec![Arc::new(NmdPlugin)])
                .build()?,
        );

        let [chromosome, position, reference, alternative] =
            spdi.split(':').collect::<Vec<_>>()[..]
        else {
            unreachable!()
        };
        let anns = predictor
            .predict(&VcfVariant {
                chromosome: chromosome.into(),
                position: position.parse()?,
                reference: reference.into(),
                alternative: alternative.into(),
            })?
            .unwrap();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_007294.4")
            .unwrap();

        assert!(ann.consequences.contains(&consequence), "{ann:?}");
        assert_eq!(
            ann.custom_fields.get(ANN_NMD).cloned().flatten().as_deref(),
            expected
        );

        Ok(())
    }
}
