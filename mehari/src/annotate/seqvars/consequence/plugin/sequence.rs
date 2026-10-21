use super::{TranscriptPlugin, TxAnnotation};
use crate::annotate::seqvars::consequence::VcfVariant;
use crate::annotate::seqvars::consequence::terms::{
    ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_TX_SEQ_ALT, ANN_TX_SEQ_REF,
};
use crate::errors::SeqvarsError;

/// Reports one transcript or protein sequence, without (`Ref`) or with (`Alt`) the variant.
///
/// Alternative sequences are only reported for coding variants.
#[derive(Debug, Clone, Copy)]
pub enum SequencePlugin {
    TxRef,
    TxAlt,
    AaRef,
    AaAlt,
}

impl SequencePlugin {
    fn column(self) -> &'static str {
        match self {
            Self::TxRef => ANN_TX_SEQ_REF,
            Self::TxAlt => ANN_TX_SEQ_ALT,
            Self::AaRef => ANN_AA_SEQ_REF,
            Self::AaAlt => ANN_AA_SEQ_ALT,
        }
    }
}

impl TranscriptPlugin for SequencePlugin {
    fn columns(&self) -> Vec<String> {
        vec![self.column().into()]
    }

    fn annotate(
        &self,
        _variants: &[VcfVariant],
        txs: &mut Vec<TxAnnotation>,
    ) -> Result<(), SeqvarsError> {
        for TxAnnotation { context, ann } in txs {
            if context.var_c.is_none() {
                continue;
            }
            let sequence = match self {
                Self::TxRef => context
                    .ref_data()
                    .map(|data| data.transcript_sequence.to_string()),
                Self::AaRef => context.ref_data().map(|data| data.aa_sequence.to_string()),
                Self::TxAlt => context
                    .alt_data()
                    .map(|data| data.transcript_sequence.clone()),
                Self::AaAlt => context.alt_data().map(|data| data.aa_sequence.clone()),
            };
            if let Some(sequence) = sequence {
                ann.custom_fields
                    .insert(self.column().into(), Some(sequence));
            }
        }
        Ok(())
    }
}
