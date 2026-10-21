use super::{TranscriptPlugin, TxAnnotation, TxContext};
use crate::annotate::seqvars::consequence::VcfVariant;
use crate::annotate::seqvars::consequence::terms::{
    ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_TX_SEQ_ALT, ANN_TX_SEQ_REF,
};
use crate::annotate::seqvars::provider::PbsTranscriptExt;
use crate::errors::SeqvarsError;

/// Reports one transcript or protein sequence, without (`Ref`) or with (`Alt`) the variant.
///
/// Alternative sequences are reported for variants in the CDS or a UTR, unless the variant changes
/// splicing. Transcript sequences leave out the `A` bases that `db create` appends to complete a
/// stop codon.
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
                    .map(|data| without_padding(context, &data.transcript_sequence)),
                Self::AaRef => context.ref_data().map(|data| data.aa_sequence.to_string()),
                Self::TxAlt => context
                    .alt_data()
                    .map(|data| without_padding(context, &data.transcript_sequence)),
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

/// `seq` without the `A` bases that `db create` appended to the stored sequence.
fn without_padding(context: &TxContext, seq: &str) -> String {
    match (context.tx(), context.ref_data()) {
        (Some(tx), Some(ref_data)) => {
            tx.without_stop_codon_padding(seq, ref_data.transcript_sequence.len(), tx.tx_len())
        }
        _ => seq.to_string(),
    }
}
