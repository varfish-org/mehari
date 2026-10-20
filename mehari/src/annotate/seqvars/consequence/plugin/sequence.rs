use super::{TranscriptPlugin, TxAnnotation};
use crate::annotate::seqvars::consequence::terms::{
    ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_TX_SEQ_ALT, ANN_TX_SEQ_REF,
};
use crate::annotate::seqvars::consequence::{SequenceReporting, VcfVariant};
use crate::errors::SeqvarsError;

/// Reports reference and/or alternative transcript and protein sequences.
#[derive(Debug)]
pub struct SequencePlugin {
    cdna: SequenceReporting,
    protein: SequenceReporting,
}

impl SequencePlugin {
    /// `None` if no sequence is to be reported.
    pub fn from_settings(cdna: SequenceReporting, protein: SequenceReporting) -> Option<Self> {
        (cdna != SequenceReporting::None || protein != SequenceReporting::None)
            .then_some(Self { cdna, protein })
    }
}

impl TranscriptPlugin for SequencePlugin {
    fn columns(&self) -> Vec<String> {
        [
            (self.cdna.includes_ref(), ANN_TX_SEQ_REF),
            (self.cdna.includes_alt(), ANN_TX_SEQ_ALT),
            (self.protein.includes_ref(), ANN_AA_SEQ_REF),
            (self.protein.includes_alt(), ANN_AA_SEQ_ALT),
        ]
        .into_iter()
        .filter(|(enabled, _)| *enabled)
        .map(|(_, column)| column.to_string())
        .collect()
    }

    fn annotate(
        &self,
        _variants: &[VcfVariant],
        txs: &mut [TxAnnotation<'_>],
    ) -> Result<(), SeqvarsError> {
        for TxAnnotation { context, ann } in txs {
            if context.var_c.is_none() {
                continue;
            }
            let Some(ref_data) = context.ref_data() else {
                continue;
            };
            let fields = &mut ann.custom_fields;
            if self.cdna.includes_ref() {
                fields.insert(
                    ANN_TX_SEQ_REF.into(),
                    Some(ref_data.transcript_sequence.to_string()),
                );
            }
            if self.protein.includes_ref() {
                fields.insert(
                    ANN_AA_SEQ_REF.into(),
                    Some(ref_data.aa_sequence.to_string()),
                );
            }
            if (self.cdna.includes_alt() || self.protein.includes_alt())
                && let Some(alt_data) = context.alt_data()
            {
                if self.cdna.includes_alt() {
                    fields.insert(
                        ANN_TX_SEQ_ALT.into(),
                        Some(alt_data.transcript_sequence.to_string()),
                    );
                }
                if self.protein.includes_alt() {
                    fields.insert(
                        ANN_AA_SEQ_ALT.into(),
                        Some(alt_data.aa_sequence.to_string()),
                    );
                }
            }
        }
        Ok(())
    }
}
