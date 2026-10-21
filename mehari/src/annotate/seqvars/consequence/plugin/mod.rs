//! Plugins that fill custom ANN columns per transcript.

use std::cell::OnceCell;
use std::fmt::Debug;
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::Arc;

use hgvs::mapper::altseq::{
    AltSeqBuilder, AltTranscriptData, RefTranscriptData, ref_transcript_data_cached,
};
use hgvs::parser::HgvsVariant;
use itertools::Itertools;
use strum::VariantArray;

use super::VcfVariant;
use super::terms::AnnField;
use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::db::TranscriptDatabase;
use crate::errors::SeqvarsError;
use crate::pbs::txs::{Transcript, TxSeqDatabase};

mod compound;
mod nmd;
mod sequence;
#[cfg(feature = "wasm")]
mod wasm;

pub use compound::CompoundPlugin;
pub use nmd::{NmdPlugin, NmdPrediction};
pub use sequence::SequencePlugin;
#[cfg(feature = "wasm")]
pub use wasm::WasmPlugin;

/// A plugin to run: a built-in by name, or a WebAssembly component by its path, followed by its
/// options as `,key=value`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PluginSpec {
    Builtin(BuiltinPlugin),
    Wasm {
        path: PathBuf,
        options: Vec<(String, String)>,
    },
}

impl FromStr for PluginSpec {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut parts = s.split(',');
        let name = parts.next().unwrap_or_default();
        let options = parts
            .map(|option| {
                option
                    .split_once('=')
                    .map(|(key, value)| (key.to_string(), value.to_string()))
                    .ok_or_else(|| {
                        format!(
                            "plugin option {option:?} is not key=value; \
                             repeat the key for several values, e.g. genes=A,genes=B"
                        )
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        if let Ok(plugin) = name.parse() {
            if !options.is_empty() {
                return Err(format!("built-in plugin {name} takes no options"));
            }
            Ok(Self::Builtin(plugin))
        } else if name.ends_with(".wasm") {
            Ok(Self::Wasm {
                path: name.into(),
                options,
            })
        } else {
            Err(format!(
                "unknown plugin {name:?}, expected one of {} or a path ending in .wasm",
                BuiltinPlugin::VARIANTS.iter().join(", ")
            ))
        }
    }
}

impl PluginSpec {
    pub fn build(&self, run: &RunInfo) -> Result<Arc<dyn TranscriptPlugin>, SeqvarsError> {
        match self {
            Self::Builtin(plugin) => Ok(plugin.build(run.vep_compatible)),
            #[cfg(feature = "wasm")]
            Self::Wasm { path, options } => {
                Ok(Arc::new(WasmPlugin::load(path, options, run.clone())?))
            }
            #[cfg(not(feature = "wasm"))]
            Self::Wasm { path, .. } => Err(SeqvarsError::Plugin(format!(
                "{}: mehari was built without the `wasm` feature",
                path.display()
            ))),
        }
    }
}

/// What plugins learn about the run.
#[derive(Debug, Clone, Default)]
pub struct RunInfo {
    /// Genome assembly of the transcript database, lowercase, e.g. `grch38`.
    pub assembly: String,
    /// Mimic VEP, see `--vep-compatibility-mode`.
    pub vep_compatible: bool,
}

impl RunInfo {
    pub fn new(tx_db: &TxSeqDatabase, vep_compatible: bool) -> Self {
        Self {
            // Older databases spell it `GRCh38`.
            assembly: tx_db.assembly().to_ascii_lowercase(),
            vep_compatible,
        }
    }
}

/// Plugins to select by name.
#[derive(
    Debug, Copy, Clone, PartialEq, Eq, VariantArray, parse_display::FromStr, parse_display::Display,
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
    /// `variants` holds the input variant as given in the VCF, or the variants of a compound
    /// group. The HGVS descriptions in `ann` and `context.var_c` follow HGVS normalization, which
    /// shifts towards the transcript's 3' end and so depends on its strand.
    /// Writes into `ann.custom_fields`; keys not listed in `columns` are not output. `txs` is a
    /// `Vec` so that a plugin can move the annotations, e.g. to lend them to a WebAssembly
    /// component.
    fn annotate(
        &self,
        variants: &[VcfVariant],
        txs: &mut Vec<TxAnnotation>,
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
    /// The variant on the transcript in `n.` coordinates.
    pub var_n: Option<HgvsVariant>,
    ref_data: OnceCell<Option<RefTranscriptData>>,
    alt_data: OnceCell<Option<AltTranscriptData>>,
}

impl TxContext {
    pub(crate) fn new(
        provider: Arc<MehariProvider>,
        tx: Option<&Transcript>,
        var_c: Option<HgvsVariant>,
        var_n: Option<HgvsVariant>,
    ) -> Self {
        Self {
            provider,
            tx_id: tx.map(|tx| tx.id.clone()),
            var_c,
            var_n,
            ref_data: OnceCell::new(),
            alt_data: OnceCell::new(),
        }
    }

    /// Reuse reference data the caller has computed already.
    pub(crate) fn with_ref_data(mut self, ref_data: RefTranscriptData) -> Self {
        self.ref_data = OnceCell::from(Some(ref_data));
        self
    }

    /// 1-based transcript positions of the first and the last variant base, if both lie in
    /// exons. Insertions span the two bases they lie between.
    pub fn tx_span(&self) -> Option<(i32, i32)> {
        let HgvsVariant::TxVariant { loc_edit, .. } = self.var_n.as_ref()? else {
            return None;
        };
        let loc = loc_edit.loc.inner();
        let exonic = |offset: Option<i32>| offset.unwrap_or(0) == 0;
        (exonic(loc.start.offset) && exonic(loc.end.offset))
            .then_some((loc.start.base, loc.end.base))
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn parse_plugin_spec() {
        assert_eq!("nmd".parse(), Ok(PluginSpec::Builtin(BuiltinPlugin::Nmd)));
        assert_eq!(
            "plugins/nmd.wasm".parse(),
            Ok(PluginSpec::Wasm {
                path: "plugins/nmd.wasm".into(),
                options: vec![],
            })
        );
        assert_eq!(
            "x.wasm,data=/d,panel=a=b".parse(),
            Ok(PluginSpec::Wasm {
                path: "x.wasm".into(),
                options: vec![("data".into(), "/d".into()), ("panel".into(), "a=b".into())],
            })
        );
        assert!("nmdd".parse::<PluginSpec>().is_err());
        assert!("nmd,a=b".parse::<PluginSpec>().is_err());
        assert!("x.wasm,flag".parse::<PluginSpec>().is_err());
        assert_eq!(
            "x.wasm,genes=A,genes=B".parse(),
            Ok(PluginSpec::Wasm {
                path: "x.wasm".into(),
                options: vec![("genes".into(), "A".into()), ("genes".into(), "B".into())],
            })
        );
        let error = "x.wasm,genes=A,B".parse::<PluginSpec>().unwrap_err();
        assert!(error.contains("repeat the key"), "{error}");
    }

    #[test]
    fn run_info_assembly_is_lowercase() {
        let tx_db = TxSeqDatabase {
            source_version: vec![crate::pbs::txs::SourceVersion {
                assembly: "GRCh38".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        assert_eq!(RunInfo::new(&tx_db, false).assembly, "grch38");
    }
}
