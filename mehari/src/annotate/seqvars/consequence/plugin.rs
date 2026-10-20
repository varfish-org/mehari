//! Plugins that fill custom ANN columns per transcript.

use std::cell::OnceCell;
use std::fmt::Debug;
use std::sync::Arc;

use hgvs::mapper::altseq::{
    AltSeqBuilder, AltTranscriptData, RefTranscriptData, ref_transcript_data_cached,
};
use hgvs::parser::HgvsVariant;

use super::SequenceReporting;
use super::terms::{ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_TX_SEQ_ALT, ANN_TX_SEQ_REF, AnnField};
use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::errors::SeqvarsError;
use crate::pbs::txs::Transcript;

/// Fills custom ANN columns for the transcripts of a variant.
pub trait TranscriptPlugin: Debug + Send + Sync {
    /// The columns this plugin fills, in output order.
    fn columns(&self) -> Vec<String>;

    /// Called once per variant with all of its reported transcripts.
    ///
    /// Writes into `ann.custom_fields`; keys not listed in `columns` are not output.
    fn annotate(&self, txs: &mut [TxAnnotation<'_>]) -> Result<(), SeqvarsError>;
}

/// The annotation of one transcript and the inputs it was derived from.
pub struct TxAnnotation<'a> {
    pub context: TxContext<'a>,
    pub ann: AnnField,
}

impl AsRef<AnnField> for TxAnnotation<'_> {
    fn as_ref(&self) -> &AnnField {
        &self.ann
    }
}

/// Inputs for annotating one transcript.
///
/// Sequence data is computed on first use and then shared between plugins.
pub struct TxContext<'a> {
    provider: &'a Arc<MehariProvider>,
    pub tx: &'a Transcript,
    /// The variant on the transcript: `c.` for coding, `n.` for non-coding transcripts,
    /// `None` if the variant lies up- or downstream.
    pub var_c: Option<HgvsVariant>,
    ref_data: OnceCell<Option<RefTranscriptData>>,
    alt_data: OnceCell<Option<AltTranscriptData>>,
}

impl<'a> TxContext<'a> {
    pub(crate) fn new(
        provider: &'a Arc<MehariProvider>,
        tx: &'a Transcript,
        var_c: Option<HgvsVariant>,
    ) -> Self {
        Self {
            provider,
            tx,
            var_c,
            ref_data: OnceCell::new(),
            alt_data: OnceCell::new(),
        }
    }

    /// Reuse reference data the caller has computed already.
    pub(crate) fn with_ref_data(mut self, ref_data: RefTranscriptData) -> Self {
        self.ref_data = OnceCell::from(Some(ref_data));
        self
    }

    /// Reference transcript and protein sequence.
    pub fn ref_data(&self) -> Option<&RefTranscriptData> {
        self.ref_data
            .get_or_init(|| {
                ref_transcript_data_cached(self.provider.clone(), &self.tx.id, None).ok()
            })
            .as_ref()
    }

    /// Transcript and protein sequence with the variant applied; coding variants only.
    pub fn alt_data(&self) -> Option<&AltTranscriptData> {
        self.alt_data
            .get_or_init(|| {
                let var_c = self
                    .var_c
                    .as_ref()
                    .filter(|var| matches!(var, HgvsVariant::CdsVariant { .. }))?;
                AltSeqBuilder::new(var_c.clone(), self.ref_data()?)
                    .build_altseq()
                    .ok()?
                    .into_iter()
                    .next()
            })
            .as_ref()
    }
}

/// Reports reference and/or alternative transcript and protein sequences.
#[derive(Debug)]
pub(crate) struct SequencePlugin {
    pub(crate) cdna: SequenceReporting,
    pub(crate) protein: SequenceReporting,
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

    fn annotate(&self, txs: &mut [TxAnnotation<'_>]) -> Result<(), SeqvarsError> {
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
