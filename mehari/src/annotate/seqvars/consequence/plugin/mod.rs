//! Plugins that fill custom ANN columns per transcript.

use std::cell::OnceCell;
use std::fmt::Debug;
use std::sync::Arc;

use hgvs::mapper::altseq::{
    AltSeqBuilder, AltTranscriptData, RefTranscriptData, ref_transcript_data_cached,
};
use hgvs::parser::HgvsVariant;

use super::VcfVariant;
use super::terms::AnnField;
use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::errors::SeqvarsError;
use crate::pbs::txs::Transcript;

mod compound;
mod nmd;
mod sequence;

pub use compound::CompoundPlugin;
pub use nmd::NmdPlugin;
pub use sequence::SequencePlugin;

/// Fills custom ANN columns for the transcripts of a variant.
pub trait TranscriptPlugin: Debug + Send + Sync {
    /// The columns this plugin fills, in output order.
    fn columns(&self) -> Vec<String>;

    /// Called once per prediction with all of its reported annotations.
    ///
    /// `variants` holds the input variant, or the variants of a compound group.
    /// Writes into `ann.custom_fields`; keys not listed in `columns` are not output.
    fn annotate(
        &self,
        variants: &[VcfVariant],
        txs: &mut [TxAnnotation<'_>],
    ) -> Result<(), SeqvarsError>;
}

/// The annotation of one transcript and the inputs it was derived from.
pub struct TxAnnotation<'a> {
    pub context: TxContext<'a>,
    pub ann: AnnField,
}

/// Inputs for annotating one transcript.
///
/// Sequence data is computed on first use and then shared between plugins.
pub struct TxContext<'a> {
    provider: &'a Arc<MehariProvider>,
    /// `None` for intergenic variants.
    pub tx: Option<&'a Transcript>,
    /// The variant on the transcript: `c.` for coding, `n.` for non-coding transcripts,
    /// `None` if the variant lies up- or downstream.
    pub var_c: Option<HgvsVariant>,
    ref_data: OnceCell<Option<RefTranscriptData>>,
    alt_data: OnceCell<Option<AltTranscriptData>>,
}

impl<'a> TxContext<'a> {
    pub(crate) fn new(
        provider: &'a Arc<MehariProvider>,
        tx: Option<&'a Transcript>,
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
                let tx = self.tx?;
                ref_transcript_data_cached(self.provider.clone(), &tx.id, None).ok()
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
