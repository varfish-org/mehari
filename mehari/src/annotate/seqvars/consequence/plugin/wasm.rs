//! Plugins compiled to WebAssembly components, implementing the `plugin` world of
//! `wit/plugin.wit`.

use std::path::{Path, PathBuf};
use std::sync::Mutex;

use itertools::Itertools;
use wasmtime::component::{Component, HasSelf, Linker, Resource, ResourceTable};
use wasmtime::{Engine, Store};
use wasmtime_wasi::{WasiCtx, WasiCtxBuilder, WasiCtxView, WasiView};

use super::{TranscriptPlugin, TxAnnotation, TxContext};
use crate::annotate::seqvars::consequence::VcfVariant;
use crate::annotate::seqvars::consequence::terms::AnnField;
use crate::errors::SeqvarsError;
use crate::pbs::txs::Strand;

mod bindings {
    wasmtime::component::bindgen!({
        path: "wit",
        world: "plugin",
        with: {
            "mehari:plugin/types.contexts": super::TxContexts,
        },
    });
}

use bindings::Plugin;
use bindings::mehari::plugin::types;

/// A plugin loaded from a WebAssembly component.
///
/// Plugins get no file, network, or environment access; stderr is passed through so their
/// panic messages show.
#[derive(derivative::Derivative)]
#[derivative(Debug)]
pub struct WasmPlugin {
    path: PathBuf,
    #[derivative(Debug = "ignore")]
    engine: Engine,
    #[derivative(Debug = "ignore")]
    pre: bindings::PluginPre<HostState>,
    columns: Vec<String>,
    /// Instances not in use; a call takes one or creates it, so there is one per thread.
    #[derivative(Debug = "ignore")]
    instances: Mutex<Vec<Instance>>,
}

type Instance = (Store<HostState>, Plugin);

/// The contexts of one call, lent to the plugin.
pub struct TxContexts(Vec<TxContext>);

struct HostState {
    table: ResourceTable,
    wasi: WasiCtx,
}

impl WasiView for HostState {
    fn ctx(&mut self) -> WasiCtxView<'_> {
        WasiCtxView {
            ctx: &mut self.wasi,
            table: &mut self.table,
        }
    }
}

impl WasmPlugin {
    /// Compile the component at `path` and query its columns.
    pub fn load(path: &Path) -> Result<Self, SeqvarsError> {
        let error = |e: wasmtime::Error| plugin_error(path, e);
        let engine = Engine::default();
        let component = Component::from_file(&engine, path).map_err(error)?;
        let mut linker = Linker::new(&engine);
        wasmtime_wasi::p2::add_to_linker_sync(&mut linker).map_err(error)?;
        Plugin::add_to_linker::<_, HasSelf<_>>(&mut linker, |state| state).map_err(error)?;
        let pre = bindings::PluginPre::new(linker.instantiate_pre(&component).map_err(error)?)
            .map_err(error)?;

        let (mut store, plugin) = instantiate(&engine, &pre).map_err(error)?;
        let columns = plugin.call_columns(&mut store).map_err(error)?;

        Ok(Self {
            path: path.into(),
            engine,
            pre,
            columns,
            instances: Mutex::new(vec![(store, plugin)]),
        })
    }

    fn checkout(&self) -> Result<Instance, SeqvarsError> {
        let instance = self
            .instances
            .lock()
            .map_err(|_| SeqvarsError::Plugin(format!("{}: poisoned lock", self.path.display())))?
            .pop();
        match instance {
            Some(instance) => Ok(instance),
            None => instantiate(&self.engine, &self.pre).map_err(|e| plugin_error(&self.path, e)),
        }
    }

    /// Hand the transcript contexts to the plugin for the call and take them back afterwards.
    fn call(
        &self,
        (store, plugin): &mut Instance,
        variants: &[VcfVariant],
        txs: &mut [TxAnnotation],
    ) -> Result<Vec<Vec<Option<String>>>, SeqvarsError> {
        let variants = variants
            .iter()
            .map(|var| types::VcfVariant {
                chromosome: var.chromosome.clone(),
                position: var.position,
                reference: var.reference.clone(),
                alternative: var.alternative.clone(),
            })
            .collect_vec();

        let annotations = txs.iter().map(|tx| annotation(&tx.ann)).collect_vec();
        let contexts = txs
            .iter_mut()
            .map(|tx| {
                let placeholder = TxContext::new(tx.context.provider.clone(), None, None);
                std::mem::replace(&mut tx.context, placeholder)
            })
            .collect();
        let handle = store
            .data_mut()
            .table
            .push(TxContexts(contexts))
            .map_err(|e| plugin_error(&self.path, e.into()))?;

        let result = plugin.call_annotate(
            &mut *store,
            &variants,
            &annotations,
            Resource::new_borrow(handle.rep()),
        );

        if let Ok(TxContexts(contexts)) = store.data_mut().table.delete(handle) {
            for (tx, context) in txs.iter_mut().zip(contexts) {
                tx.context = context;
            }
        }
        result.map_err(|e| plugin_error(&self.path, e))
    }
}

impl TranscriptPlugin for WasmPlugin {
    fn columns(&self) -> Vec<String> {
        self.columns.clone()
    }

    fn annotate(
        &self,
        variants: &[VcfVariant],
        txs: &mut [TxAnnotation],
    ) -> Result<(), SeqvarsError> {
        let mut instance = self.checkout()?;
        // An instance that trapped is not reused.
        let rows = self.call(&mut instance, variants, txs)?;
        if let Ok(mut instances) = self.instances.lock() {
            instances.push(instance);
        }

        if rows.len() != txs.len() || rows.iter().any(|row| row.len() != self.columns.len()) {
            return Err(SeqvarsError::Plugin(format!(
                "{}: expected {} rows of {} values",
                self.path.display(),
                txs.len(),
                self.columns.len()
            )));
        }
        for (tx, row) in txs.iter_mut().zip(rows) {
            for (column, value) in self.columns.iter().zip(row) {
                if value.is_some() {
                    tx.ann.custom_fields.insert(column.clone(), value);
                }
            }
        }
        Ok(())
    }
}

fn instantiate(
    engine: &Engine,
    pre: &bindings::PluginPre<HostState>,
) -> wasmtime::Result<Instance> {
    let state = HostState {
        table: ResourceTable::new(),
        wasi: WasiCtxBuilder::new().inherit_stderr().build(),
    };
    let mut store = Store::new(engine, state);
    let plugin = pre.instantiate(&mut store)?;
    Ok((store, plugin))
}

fn plugin_error(path: &Path, error: wasmtime::Error) -> SeqvarsError {
    SeqvarsError::Plugin(format!("{}: {error:#}", path.display()))
}

fn annotation(ann: &AnnField) -> types::Annotation {
    types::Annotation {
        feature_id: ann.feature_id.clone(),
        consequences: ann.consequences.iter().map(|csq| csq.to_string()).collect(),
        hgvs_g: ann.hgvs_g.clone(),
        hgvs_n: ann.hgvs_n.clone(),
        hgvs_c: ann.hgvs_c.clone(),
        hgvs_p: ann.hgvs_p.clone(),
        cdna_pos: ann.cdna_pos.as_ref().map(|pos| pos.ord),
    }
}

impl types::Host for HostState {}

impl HostState {
    fn context(&self, contexts: &Resource<TxContexts>, index: u32) -> Option<&TxContext> {
        self.table.get(contexts).ok()?.0.get(index as usize)
    }
}

impl types::HostContexts for HostState {
    fn transcript(
        &mut self,
        contexts: Resource<TxContexts>,
        index: u32,
    ) -> Option<types::Transcript> {
        let tx = self.context(&contexts, index)?.tx()?;
        let alignment = tx.genome_alignments.first()?;
        let exons = alignment
            .exons
            .iter()
            .sorted_by_key(|exon| exon.alt_cds_start_i)
            .map(|exon| {
                // Despite their names, `alt_cds_*_i` are the exon's 1-based transcript positions.
                Some(types::Exon {
                    tx_start: exon.alt_cds_start_i?,
                    tx_end: exon.alt_cds_end_i?,
                    genome_start: exon.alt_start_i + 1,
                    genome_end: exon.alt_end_i,
                })
            })
            .collect::<Option<Vec<_>>>()?;
        Some(types::Transcript {
            id: tx.id.clone(),
            strand: if alignment.strand == Strand::Minus as i32 {
                types::Strand::Minus
            } else {
                types::Strand::Plus
            },
            exons,
        })
    }

    fn reference(
        &mut self,
        contexts: Resource<TxContexts>,
        index: u32,
    ) -> Option<types::Sequences> {
        let data = self.context(&contexts, index)?.ref_data()?;
        Some(types::Sequences {
            transcript: data.transcript_sequence.to_string(),
            protein: data.aa_sequence.to_string(),
            cds_start: data.cds_start,
        })
    }

    fn alternative(
        &mut self,
        contexts: Resource<TxContexts>,
        index: u32,
    ) -> Option<types::Sequences> {
        let data = self.context(&contexts, index)?.alt_data()?;
        Some(types::Sequences {
            transcript: data.transcript_sequence.clone(),
            protein: data.aa_sequence.clone(),
            cds_start: data.cds_start,
        })
    }

    fn drop(&mut self, contexts: Resource<TxContexts>) -> wasmtime::Result<()> {
        self.table.delete(contexts)?;
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use std::sync::Arc;

    use rayon::prelude::*;

    use super::*;
    use crate::annotate::seqvars::consequence::logic::ConsequencePredictor;
    use crate::annotate::seqvars::consequence::plugin::NmdPlugin;
    use crate::annotate::seqvars::consequence::terms::ANN_NMD;
    use crate::annotate::seqvars::consequence::{ConfigBuilder, load_tx_db};
    use crate::annotate::seqvars::provider::Provider as MehariProvider;

    /// Built from `plugins/nmd`, see its README.
    const NMD_WASM: &str = "tests/data/plugins/nmd.wasm";

    #[test]
    fn load_reports_columns() -> Result<(), anyhow::Error> {
        let plugin = WasmPlugin::load(Path::new(NMD_WASM))?;
        assert_eq!(plugin.columns(), [ANN_NMD]);
        Ok(())
    }

    #[test]
    fn load_missing_file_fails() {
        assert!(WasmPlugin::load(Path::new("tests/data/plugins/missing.wasm")).is_err());
    }

    /// The reference plugin predicts the same as the built-in one, also when called from
    /// several threads.
    #[test]
    fn nmd_plugin_matches_builtin() -> Result<(), anyhow::Error> {
        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            None::<PathBuf>,
            true,
            Default::default(),
        ));
        let predictor = |plugin: Arc<dyn TranscriptPlugin>| -> Result<_, anyhow::Error> {
            Ok(ConsequencePredictor::new(
                provider.clone(),
                ConfigBuilder::default().plugins(vec![plugin]).build()?,
            ))
        };
        let builtin = predictor(Arc::new(NmdPlugin::default()))?;
        let wasm = predictor(Arc::new(WasmPlugin::load(Path::new(NMD_WASM))?))?;

        // BRCA1 stops and frameshifts at the rule boundaries, and a missense variant.
        let variants = [
            "17:41276049:A:T",
            "17:41276042:A:T",
            "17:41267765:T:A",
            "17:41267760:A:T",
            "17:41244528:G:T",
            "17:41244522:G:T",
            "17:41199696:G:A",
            "17:41199715:GA:G",
            "17:41199714:GG:G",
            "17:41199711:GG:G",
            "17:41199709:TG:T",
            "17:41197818:GG:G",
            "17:41215920:G:GA",
            "17:41197811:G:A",
            "17:41197803:A:T",
            "17:41215920:G:T",
            "17:41197701:G:C",
        ]
        .map(|spdi| {
            let [chromosome, position, reference, alternative] =
                spdi.split(':').collect::<Vec<_>>()[..]
            else {
                unreachable!()
            };
            VcfVariant {
                chromosome: chromosome.into(),
                position: position.parse().unwrap(),
                reference: reference.into(),
                alternative: alternative.into(),
            }
        });

        let nmd = |predictor: &ConsequencePredictor, variant: &VcfVariant| {
            predictor
                .predict(variant)
                .unwrap()
                .unwrap()
                .into_iter()
                .map(|ann| (ann.feature_id, ann.custom_fields.get(ANN_NMD).cloned()))
                .collect_vec()
        };
        let expected = variants.iter().map(|var| nmd(&builtin, var)).collect_vec();
        let actual = variants
            .par_iter()
            .map(|var| nmd(&wasm, var))
            .collect::<Vec<_>>();

        assert_eq!(actual, expected);
        assert!(expected.iter().flatten().any(|(_, value)| value.is_some()));
        Ok(())
    }
}
