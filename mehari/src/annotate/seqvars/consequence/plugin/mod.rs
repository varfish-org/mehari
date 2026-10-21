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
pub use nmd::{NmdPlugin, NmdPrediction};
pub use sequence::SequencePlugin;

/// Plugins to select by name, e.g. with `--plugin`.
#[derive(
    Debug,
    Copy,
    Clone,
    PartialEq,
    Eq,
    clap::ValueEnum,
    parse_display::FromStr,
    parse_display::Display,
)]
#[display(style = "kebab-case")]
pub enum BuiltinPlugin {
    /// Reference transcript sequence (`tx_sequence_ref` column).
    TxSequenceRef,
    /// Transcript sequence with the variant applied (`tx_sequence_alt` column).
    TxSequenceAlt,
    /// Reference protein sequence (`aa_sequence_ref` column).
    AaSequenceRef,
    /// Protein sequence with the variant applied (`aa_sequence_alt` column).
    AaSequenceAlt,
    /// Predict whether premature stop codons trigger nonsense-mediated decay (`nmd` column).
    Nmd,
}

impl BuiltinPlugin {
    /// `vep_compatible` makes the plugin mimic its VEP counterpart.
    pub fn build(self, vep_compatible: bool) -> Arc<dyn TranscriptPlugin> {
        match self {
            Self::TxSequenceRef => Arc::new(SequencePlugin::TxRef),
            Self::TxSequenceAlt => Arc::new(SequencePlugin::TxAlt),
            Self::AaSequenceRef => Arc::new(SequencePlugin::AaRef),
            Self::AaSequenceAlt => Arc::new(SequencePlugin::AaAlt),
            Self::Nmd => Arc::new(NmdPlugin { vep_compatible }),
        }
    }
}

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
        txs: &mut [TxAnnotation],
    ) -> Result<(), SeqvarsError>;
}

/// The annotation of one transcript and the inputs it was derived from.
pub struct TxAnnotation {
    pub context: TxContext,
    pub ann: AnnField,
}

/// Inputs for annotating one transcript.
///
/// Sequence data is computed on first use and then shared between plugins. The context owns
/// its data, so a wasm host can hand it to plugins as a resource.
pub struct TxContext {
    provider: Arc<MehariProvider>,
    tx_id: Option<String>,
    /// The variant on the transcript: `c.` for coding, `n.` for non-coding transcripts,
    /// `None` if the variant lies up- or downstream.
    pub var_c: Option<HgvsVariant>,
    ref_data: OnceCell<Option<RefTranscriptData>>,
    alt_data: OnceCell<Option<AltTranscriptData>>,
}

impl TxContext {
    pub(crate) fn new(
        provider: Arc<MehariProvider>,
        tx: Option<&Transcript>,
        var_c: Option<HgvsVariant>,
    ) -> Self {
        Self {
            provider,
            tx_id: tx.map(|tx| tx.id.clone()),
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

    /// The transcript; `None` for intergenic variants.
    pub fn tx(&self) -> Option<&Transcript> {
        self.provider.get_tx(self.tx_id.as_deref()?)
    }

    /// Reference transcript and protein sequence.
    pub fn ref_data(&self) -> Option<&RefTranscriptData> {
        self.ref_data
            .get_or_init(|| {
                let tx_id = self.tx_id.as_deref()?;
                ref_transcript_data_cached(self.provider.clone(), tx_id, None).ok()
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
