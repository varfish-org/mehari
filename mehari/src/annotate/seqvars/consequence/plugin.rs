//! Plugins that fill custom ANN columns per transcript.

use std::cell::OnceCell;
use std::fmt::Debug;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use hgvs::mapper::altseq::{
    AltSeqBuilder, AltTranscriptData, RefTranscriptData, ref_transcript_data_cached,
};
use hgvs::parser::HgvsVariant;
use itertools::Itertools;

use super::terms::{
    ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_COMPOUND_IDS, ANN_COMPOUND_VARIANTS, ANN_TX_SEQ_ALT,
    ANN_TX_SEQ_REF, AnnField,
};
use super::{SequenceReporting, VcfVariant};
use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::errors::SeqvarsError;
use crate::pbs::txs::Transcript;

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

/// Labels single variants `Single` and compound groups `comp_<n>`, plus the input variants.
///
/// Groups are numbered in prediction order.
#[derive(Debug, Default)]
pub struct CompoundPlugin {
    groups: AtomicUsize,
}

impl TranscriptPlugin for CompoundPlugin {
    fn columns(&self) -> Vec<String> {
        vec![ANN_COMPOUND_IDS.into(), ANN_COMPOUND_VARIANTS.into()]
    }

    fn annotate(
        &self,
        variants: &[VcfVariant],
        txs: &mut [TxAnnotation<'_>],
    ) -> Result<(), SeqvarsError> {
        // Only groups with annotations get a number.
        if txs.is_empty() {
            return Ok(());
        }
        let id = if variants.len() > 1 {
            format!("comp_{}", self.groups.fetch_add(1, Ordering::Relaxed) + 1)
        } else {
            "Single".to_string()
        };
        let variants = variants
            .iter()
            .map(|var| {
                format!(
                    "{}:{}:{}:{}",
                    var.chromosome, var.position, var.reference, var.alternative
                )
            })
            .join(",");
        for TxAnnotation { ann, .. } in txs {
            ann.custom_fields
                .insert(ANN_COMPOUND_IDS.into(), Some(id.clone()));
            ann.custom_fields
                .insert(ANN_COMPOUND_VARIANTS.into(), Some(variants.clone()));
        }
        Ok(())
    }
}
