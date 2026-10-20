use std::sync::atomic::{AtomicUsize, Ordering};

use itertools::Itertools;

use super::{TranscriptPlugin, TxAnnotation};
use crate::annotate::seqvars::consequence::VcfVariant;
use crate::annotate::seqvars::consequence::terms::{ANN_COMPOUND_IDS, ANN_COMPOUND_VARIANTS};
use crate::errors::SeqvarsError;

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
        txs: &mut [TxAnnotation],
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
