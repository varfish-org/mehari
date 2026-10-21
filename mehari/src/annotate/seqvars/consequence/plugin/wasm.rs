//! Plugins compiled to WebAssembly components, implementing the `plugin` world of
//! `wit/plugin.wit`.

use std::path::{Path, PathBuf};
use std::sync::Mutex;

use hgvs::sequences::TranslationTable;
use itertools::Itertools;
use wasmtime::component::{Component, HasSelf, Linker, Resource, ResourceTable};
use wasmtime::{Config, Engine, Store};
use wasmtime_wasi::{FsPerms, WasiCtx, WasiCtxBuilder, WasiCtxView, WasiView};

use super::{RunInfo, TranscriptPlugin, TxAnnotation, TxContext};
use crate::annotate::seqvars::consequence::VcfVariant;
use crate::annotate::seqvars::consequence::terms::{AnnField, Pos};
use crate::errors::SeqvarsError;
use crate::pbs::txs::{Strand, Transcript};

mod bindings {
    wasmtime::component::bindgen!({
        path: "wit",
        world: "plugin",
        imports: { default: trappable },
        with: {
            "mehari:plugin/types.annotations": super::Annotations,
        },
    });
}

use bindings::Plugin;
use bindings::mehari::plugin::{host, types};

/// The option naming a directory that the plugin may read, at [`DATA_DIR`].
pub const DATA_OPTION: &str = "data";
/// Where plugins find the directory given with [`DATA_OPTION`].
pub const DATA_DIR: &str = "/data";

/// A plugin loaded from a WebAssembly component.
///
/// Plugins get no network or environment access, and read access only to the directory given
/// with [`DATA_OPTION`]; stderr is passed through so their errors show.
#[derive(derivative::Derivative)]
#[derivative(Debug)]
pub struct WasmPlugin {
    path: PathBuf,
    /// Passed to `init`.
    options: Vec<(String, String)>,
    data: Option<PathBuf>,
    run: RunInfo,
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

/// The annotations of one call, lent to the plugin. One resource for all annotations rather
/// than one each, as wit-bindgen cannot lift lists of borrowed handles.
pub struct Annotations(Vec<TxAnnotation>);

struct HostState {
    table: ResourceTable,
    wasi: WasiCtx,
    run: RunInfo,
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
    /// Compile the component at `path` and initialize it with `options`.
    pub fn load(
        path: &Path,
        options: &[(String, String)],
        run: RunInfo,
    ) -> Result<Self, SeqvarsError> {
        let error = |e: wasmtime::Error| plugin_error(path, e);
        let (data, options): (Vec<_>, Vec<_>) = options
            .iter()
            .cloned()
            .partition(|(key, _)| key == DATA_OPTION);
        let data = match &data[..] {
            [] => None,
            [(_, dir)] => Some(dir.into()),
            _ => return Err(plugin_message(path, "option `data` given more than once")),
        };

        // Plugins print their errors to stderr; wasm backtraces would bury them.
        let engine = Engine::new(Config::new().wasm_backtrace_max_frames(None)).map_err(error)?;
        let component = Component::from_file(&engine, path).map_err(error)?;
        let mut linker = Linker::new(&engine);
        wasmtime_wasi::p2::add_to_linker_sync(&mut linker).map_err(error)?;
        Plugin::add_to_linker::<_, HasSelf<_>>(&mut linker, |state| state).map_err(error)?;
        let pre = bindings::PluginPre::new(linker.instantiate_pre(&component).map_err(error)?)
            .map_err(error)?;

        let mut plugin = Self {
            path: path.into(),
            options,
            data,
            run,
            engine,
            pre,
            columns: Vec::new(),
            instances: Mutex::default(),
        };
        let (instance, columns) = outside_runtime(path, || plugin.instantiate())?;
        plugin.columns = columns;
        plugin.instances = Mutex::new(vec![instance]);
        Ok(plugin)
    }

    /// A new instance and the columns its `init` returned.
    fn instantiate(&self) -> Result<(Instance, Vec<String>), SeqvarsError> {
        let error = |e: wasmtime::Error| plugin_error(&self.path, e);
        let mut wasi = WasiCtxBuilder::new();
        wasi.inherit_stderr();
        if let Some(dir) = &self.data {
            wasi.preopened_dir(dir, DATA_DIR, FsPerms::ReadOnly)
                .map_err(|e| {
                    plugin_message(
                        &self.path,
                        &format!("data directory {}: {e:#}", dir.display()),
                    )
                })?;
        }
        let state = HostState {
            table: ResourceTable::new(),
            wasi: wasi.build(),
            run: self.run.clone(),
        };
        let mut store = Store::new(&self.engine, state);
        let plugin = self.pre.instantiate(&mut store).map_err(error)?;
        let columns = plugin
            .call_init(&mut store, &self.options)
            .map_err(error)?
            .map_err(|message| plugin_message(&self.path, &message))?;
        Ok(((store, plugin), columns))
    }

    fn checkout(&self) -> Result<Instance, SeqvarsError> {
        let instance = self
            .instances
            .lock()
            .map_err(|_| plugin_message(&self.path, "poisoned lock"))?
            .pop();
        if let Some(instance) = instance {
            return Ok(instance);
        }
        let (instance, columns) = outside_runtime(&self.path, || self.instantiate())?;
        if columns != self.columns {
            return Err(plugin_message(&self.path, "init returned other columns"));
        }
        Ok(instance)
    }

    /// Lend the annotations to the plugin for the call and take them back afterwards.
    fn call(
        &self,
        (store, plugin): &mut Instance,
        variants: &[VcfVariant],
        txs: &mut Vec<TxAnnotation>,
    ) -> Result<Vec<Vec<Option<String>>>, SeqvarsError> {
        let error = |e: wasmtime::Error| plugin_error(&self.path, e);
        let variants = variants
            .iter()
            .map(|var| types::VcfVariant {
                chromosome: var.chromosome.clone(),
                position: var.position,
                reference: var.reference.clone(),
                alternative: var.alternative.clone(),
            })
            .collect_vec();

        let handle = store
            .data_mut()
            .table
            .push(Annotations(std::mem::take(txs)))
            .map_err(|e| error(e.into()))?;
        let result =
            plugin.call_annotate(&mut *store, &variants, Resource::new_borrow(handle.rep()));
        let Annotations(lent) = store
            .data_mut()
            .table
            .delete(handle)
            .map_err(|e| error(e.into()))?;
        *txs = lent;

        result
            .map_err(error)?
            .map_err(|message| plugin_message(&self.path, &message))
    }
}

impl TranscriptPlugin for WasmPlugin {
    fn columns(&self) -> Vec<String> {
        self.columns.clone()
    }

    fn annotate(
        &self,
        variants: &[VcfVariant],
        txs: &mut Vec<TxAnnotation>,
    ) -> Result<(), SeqvarsError> {
        let mut instance = self.checkout()?;
        // An instance that trapped is not reused.
        let rows = outside_runtime(&self.path, || self.call(&mut instance, variants, txs))?;
        if let Ok(mut instances) = self.instances.lock() {
            instances.push(instance);
        }

        if rows.len() != txs.len() || rows.iter().any(|row| row.len() != self.columns.len()) {
            return Err(plugin_message(
                &self.path,
                &format!(
                    "expected {} rows of {} values",
                    txs.len(),
                    self.columns.len()
                ),
            ));
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

/// Run `f`, which calls into a plugin, outside any tokio runtime: wasmtime-wasi blocks on tokio
/// for I/O, which panics on a thread that drives a runtime.
fn outside_runtime<R: Send>(
    path: &Path,
    f: impl FnOnce() -> Result<R, SeqvarsError> + Send,
) -> Result<R, SeqvarsError> {
    match tokio::runtime::Handle::try_current() {
        Err(_) => f(),
        Ok(handle) if handle.runtime_flavor() == tokio::runtime::RuntimeFlavor::CurrentThread => {
            std::thread::scope(|scope| scope.spawn(f).join())
                .unwrap_or_else(|_| Err(plugin_message(path, "plugin thread panicked")))
        }
        Ok(_) => tokio::task::block_in_place(f),
    }
}

fn plugin_error(path: &Path, error: wasmtime::Error) -> SeqvarsError {
    plugin_message(path, &format!("{error:#}"))
}

fn plugin_message(path: &Path, message: &str) -> SeqvarsError {
    SeqvarsError::Plugin(format!("{}: {message}", path.display()))
}

impl host::Host for HostState {
    fn assembly(&mut self) -> wasmtime::Result<String> {
        Ok(self.run.assembly.clone())
    }

    fn version(&mut self) -> wasmtime::Result<String> {
        Ok(crate::built_info::PKG_VERSION.into())
    }

    fn vep_compatible(&mut self) -> wasmtime::Result<bool> {
        Ok(self.run.vep_compatible)
    }
}

impl types::Host for HostState {}

impl HostState {
    fn tx(
        &self,
        annotations: &Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<&TxAnnotation> {
        let Annotations(txs) = self.table.get(annotations)?;
        txs.get(index as usize).ok_or_else(|| {
            wasmtime::format_err!(
                "annotation index {index} out of range, there are {}",
                txs.len()
            )
        })
    }

    fn ann(&self, annotations: &Resource<Annotations>, index: u32) -> wasmtime::Result<&AnnField> {
        Ok(&self.tx(annotations, index)?.ann)
    }

    fn context(
        &self,
        annotations: &Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<&TxContext> {
        Ok(&self.tx(annotations, index)?.context)
    }
}

impl types::HostAnnotations for HostState {
    fn len(&mut self, annotations: Resource<Annotations>) -> wasmtime::Result<u32> {
        Ok(u32::try_from(self.table.get(&annotations)?.0.len())?)
    }

    fn feature_id(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<String> {
        Ok(self.ann(&a, index)?.feature_id.clone())
    }

    fn feature_type(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<String> {
        Ok(self.ann(&a, index)?.feature_type.to_string())
    }

    fn consequences(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Vec<String>> {
        Ok(strings(&self.ann(&a, index)?.consequences))
    }

    fn impact(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<String> {
        Ok(self.ann(&a, index)?.putative_impact.to_string())
    }

    fn gene_symbol(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<String> {
        Ok(self.ann(&a, index)?.gene_symbol.clone())
    }

    fn gene_id(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<String> {
        Ok(self.ann(&a, index)?.gene_id.clone())
    }

    fn biotypes(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<Vec<String>> {
        Ok(strings(&self.ann(&a, index)?.feature_biotype))
    }

    fn tags(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<Vec<String>> {
        Ok(strings(&self.ann(&a, index)?.feature_tags))
    }

    fn rank(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Option<types::Rank>> {
        Ok(self.ann(&a, index)?.rank.as_ref().map(|rank| types::Rank {
            ord: rank.ord,
            total: rank.total,
        }))
    }

    fn hgvs_g(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<Option<String>> {
        Ok(self.ann(&a, index)?.hgvs_g.clone())
    }

    fn hgvs_n(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<Option<String>> {
        Ok(self.ann(&a, index)?.hgvs_n.clone())
    }

    fn hgvs_c(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<Option<String>> {
        Ok(self.ann(&a, index)?.hgvs_c.clone())
    }

    fn hgvs_p(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<Option<String>> {
        Ok(self.ann(&a, index)?.hgvs_p.clone())
    }

    fn cdna_pos(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Option<types::Position>> {
        Ok(self.ann(&a, index)?.cdna_pos.as_ref().map(position))
    }

    fn cds_pos(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Option<types::Position>> {
        Ok(self.ann(&a, index)?.cds_pos.as_ref().map(position))
    }

    fn protein_pos(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Option<types::Position>> {
        Ok(self.ann(&a, index)?.protein_pos.as_ref().map(position))
    }

    fn distance(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<Option<i32>> {
        Ok(self.ann(&a, index)?.distance)
    }

    fn warnings(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<Vec<String>> {
        let messages = &self.ann(&a, index)?.messages;
        Ok(messages.iter().flatten().map(ToString::to_string).collect())
    }

    fn field(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
        column: String,
    ) -> wasmtime::Result<Option<String>> {
        let fields = &self.ann(&a, index)?.custom_fields;
        Ok(fields.get(&column).cloned().flatten())
    }

    fn transcript(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Option<types::Transcript>> {
        Ok(self.context(&a, index)?.tx().and_then(transcript))
    }

    fn reference(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Option<types::Sequences>> {
        let data = self.context(&a, index)?.ref_data();
        Ok(data.map(|data| types::Sequences {
            transcript: data.transcript_sequence.to_string(),
            protein: data.aa_sequence.to_string(),
            cds_start: data.cds_start,
        }))
    }

    fn alternative(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Option<types::Sequences>> {
        let data = self.context(&a, index)?.alt_data();
        Ok(data.map(|data| types::Sequences {
            transcript: data.transcript_sequence.clone(),
            protein: data.aa_sequence.clone(),
            cds_start: data.cds_start,
        }))
    }

    fn cds_end(&mut self, a: Resource<Annotations>, index: u32) -> wasmtime::Result<Option<i32>> {
        Ok(self
            .context(&a, index)?
            .ref_data()
            .map(|data| data.cds_stop))
    }

    fn translation(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Option<types::Translation>> {
        let data = self.context(&a, index)?.ref_data();
        Ok(data.map(|data| types::Translation {
            table: match data.translation_table {
                TranslationTable::Standard => "standard",
                TranslationTable::Selenocysteine => "selenocysteine",
                TranslationTable::VertebrateMitochondrial => "vertebrate-mitochondrial",
            }
            .into(),
            exceptions: data
                .translation_exceptions
                .iter()
                .map(|exception| types::TranslationException {
                    position: exception.position,
                    amino_acid: exception.amino_acid.to_string(),
                })
                .collect(),
        }))
    }

    fn frameshift(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Option<bool>> {
        let data = self.context(&a, index)?.alt_data();
        Ok(data.map(|data| data.is_frameshift))
    }

    fn first_changed_codon(
        &mut self,
        a: Resource<Annotations>,
        index: u32,
    ) -> wasmtime::Result<Option<i32>> {
        let data = self.context(&a, index)?.alt_data();
        Ok(data.and_then(|data| data.variant_start_aa))
    }

    fn drop(&mut self, annotations: Resource<Annotations>) -> wasmtime::Result<()> {
        self.table.delete(annotations)?;
        Ok(())
    }
}

fn strings(items: &[impl ToString]) -> Vec<String> {
    items.iter().map(ToString::to_string).collect()
}

fn position(pos: &Pos) -> types::Position {
    types::Position {
        ord: pos.ord,
        total: pos.total,
    }
}

fn transcript(tx: &Transcript) -> Option<types::Transcript> {
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
                cigar: exon.cigar.clone(),
            })
        })
        .collect::<Option<Vec<_>>>()?;
    Some(types::Transcript {
        id: tx.id.clone(),
        protein: tx.protein.clone(),
        contig: alignment.contig.clone(),
        genome_build: alignment.genome_build.clone(),
        strand: if alignment.strand == Strand::Minus as i32 {
            types::Strand::Minus
        } else {
            types::Strand::Plus
        },
        genome_cds_start: alignment.cds_start.map(|start| start + 1),
        genome_cds_end: alignment.cds_end,
        exons,
    })
}

#[cfg(test)]
mod test {
    use std::collections::BTreeMap;
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
    /// Built from `plugins/nmd-python`, see its README; too large to commit.
    const NMD_PYTHON_WASM: &str = "../plugins/nmd-python/nmd.wasm";
    /// Built from `plugins/probe`, see its README.
    const PROBE_WASM: &str = "tests/data/plugins/probe.wasm";

    fn load(path: &str, options: &[(&str, &str)]) -> Result<WasmPlugin, SeqvarsError> {
        let options = options
            .iter()
            .map(|(key, value)| (key.to_string(), value.to_string()))
            .collect_vec();
        let run = RunInfo {
            assembly: "grch37".into(),
            vep_compatible: true,
        };
        WasmPlugin::load(Path::new(path), &options, run)
    }

    fn predictor(
        plugins: Vec<Arc<dyn TranscriptPlugin>>,
    ) -> Result<ConsequencePredictor, anyhow::Error> {
        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            None::<PathBuf>,
            true,
            Default::default(),
        ));
        Ok(ConsequencePredictor::new(
            provider,
            ConfigBuilder::default().plugins(plugins).build()?,
        ))
    }

    fn variant(spdi: &str) -> VcfVariant {
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
    }

    /// The probe's report for BRCA1 NM_007294.4 and the annotation it was made for; the built-in
    /// `nmd` plugin runs first.
    fn probe(
        spdi: &str,
        options: &[(&str, &str)],
    ) -> Result<(BTreeMap<String, String>, AnnField), anyhow::Error> {
        let plugins: Vec<Arc<dyn TranscriptPlugin>> = vec![
            Arc::new(NmdPlugin::default()),
            Arc::new(load(PROBE_WASM, options)?),
        ];
        let ann = predictor(plugins)?
            .predict(&variant(spdi))?
            .unwrap()
            .into_iter()
            .find(|ann| ann.feature_id == "NM_007294.4")
            .unwrap();
        let report = ann.custom_fields["probe"]
            .as_deref()
            .unwrap()
            .split(';')
            .map(|pair| {
                let (key, value) = pair.split_once('=').unwrap();
                (key.to_string(), value.to_string())
            })
            .collect();
        Ok((report, ann))
    }

    #[test]
    fn load_reports_columns() -> Result<(), anyhow::Error> {
        assert_eq!(load(NMD_WASM, &[])?.columns(), [ANN_NMD]);
        Ok(())
    }

    #[test]
    fn load_missing_file_fails() {
        assert!(load("tests/data/plugins/missing.wasm", &[]).is_err());
    }

    #[test]
    fn nmd_plugin_matches_builtin() -> Result<(), anyhow::Error> {
        assert_matches_builtin(NMD_WASM)
    }

    #[test]
    #[ignore = "needs plugins/nmd-python/nmd.wasm, see its README"]
    fn nmd_python_plugin_matches_builtin() -> Result<(), anyhow::Error> {
        assert_matches_builtin(NMD_PYTHON_WASM)
    }

    /// The NMD plugin at `path` predicts the same as the built-in one, also when called from
    /// several threads.
    fn assert_matches_builtin(path: &str) -> Result<(), anyhow::Error> {
        let builtin = predictor(vec![Arc::new(NmdPlugin::default())])?;
        let wasm = predictor(vec![Arc::new(load(path, &[])?)])?;

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
        .map(variant);

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

    /// Every function of the interface, for the missense variant c.5123C>A (p.Ala1708Glu).
    #[test]
    fn probe_reports_everything() -> Result<(), anyhow::Error> {
        let (report, ann) = probe("17:41215920:G:T", &[])?;
        let join = |items: Vec<String>| items.join("&");
        let show = |value: Option<String>| value.unwrap_or_default();
        let expected = [
            ("assembly", "grch37".to_string()),
            ("version", crate::built_info::PKG_VERSION.to_string()),
            ("vep", "true".into()),
            ("variants", "1".into()),
            // As mehari reports them.
            ("feature-id", ann.feature_id.clone()),
            ("feature-type", ann.feature_type.to_string()),
            ("consequences", join(strings(&ann.consequences))),
            ("impact", ann.putative_impact.to_string()),
            ("gene-symbol", ann.gene_symbol.clone()),
            ("gene-id", ann.gene_id.clone()),
            ("biotypes", join(strings(&ann.feature_biotype))),
            ("tags", join(strings(&ann.feature_tags))),
            ("rank", show(ann.rank.as_ref().map(ToString::to_string))),
            ("hgvs-g", show(ann.hgvs_g.clone())),
            ("hgvs-n", show(ann.hgvs_n.clone())),
            ("hgvs-c", show(ann.hgvs_c.clone())),
            ("hgvs-p", show(ann.hgvs_p.clone())),
            (
                "cdna-pos",
                show(ann.cdna_pos.as_ref().map(ToString::to_string)),
            ),
            (
                "cds-pos",
                show(ann.cds_pos.as_ref().map(ToString::to_string)),
            ),
            (
                "protein-pos",
                show(ann.protein_pos.as_ref().map(ToString::to_string)),
            ),
            ("distance", show(ann.distance.map(|d| d.to_string()))),
            ("warnings", String::new()),
            // c.1 is n.114 and chr17:41276113, the stop codon ends at c.5592.
            ("cds-end", "5705".into()),
            ("frameshift", "false".into()),
            ("first-changed-codon", "1708".into()),
            ("tx-id", "NM_007294.4".into()),
            ("tx-protein", "NP_009225.1".into()),
            ("tx-contig", "NC_000017.10".into()),
            ("tx-build", "grch37".into()),
            ("tx-strand", "Strand::Minus".into()),
            ("tx-genome-cds", "41197695-41276113".into()),
            ("tx-exons", "23".into()),
            ("tx-first-exon", "1-94:41277288-41277381:94M".into()),
            ("ref", "7091/1864/114".into()),
            ("alt", "7091/1864/114".into()),
            ("translation", "standard/".into()),
        ]
        .map(|(key, value)| (key.to_string(), value))
        .into_iter()
        .collect::<BTreeMap<_, _>>();

        assert_eq!(report, expected);
        assert_eq!(report["hgvs-p"], "p.Ala1708Glu");
        Ok(())
    }

    /// Plugins see the columns of earlier plugins.
    #[test]
    fn probe_reads_earlier_columns() -> Result<(), anyhow::Error> {
        // c.65T>A, p.Leu22Ter.
        let (report, _) = probe("17:41276049:A:T", &[("field", ANN_NMD)])?;
        assert_eq!(report["field"], "escaping");
        Ok(())
    }

    /// The data directory can be read but not written or left.
    #[rstest::rstest]
    #[case("hello.txt", "hello")]
    #[case("../secret.txt", "error")]
    fn probe_reads_data_dir(
        #[case] file: &str,
        #[case] expected: &str,
    ) -> Result<(), anyhow::Error> {
        let dir = tempfile::tempdir()?;
        let data = dir.path().join("data");
        std::fs::create_dir(&data)?;
        std::fs::write(data.join("hello.txt"), "hello")?;
        std::fs::write(dir.path().join("secret.txt"), "secret")?;
        let data_option = data.to_string_lossy();
        let options = [
            ("data", data_option.as_ref()),
            ("read", file),
            ("write", "new.txt"),
        ];
        let (report, _) = probe("17:41215920:G:T", &options)?;
        assert!(report["read"].starts_with(expected), "{}", report["read"]);
        assert_eq!(report["write"], "false");
        assert!(!data.join("new.txt").exists());
        Ok(())
    }

    #[test]
    fn probe_without_data_dir_reads_nothing() -> Result<(), anyhow::Error> {
        let (report, _) = probe("17:41215920:G:T", &[("read", "etc/passwd")])?;
        assert!(report["read"].starts_with("error"), "{}", report["read"]);
        Ok(())
    }

    /// File access works from async code, where wasmtime-wasi cannot block on its own.
    async fn probe_reads_from_async_code() -> Result<(), anyhow::Error> {
        let dir = tempfile::tempdir()?;
        std::fs::write(dir.path().join("hello.txt"), "hello")?;
        let data = dir.path().to_string_lossy();
        let options = [
            ("data", data.as_ref()),
            ("init-read", "hello.txt"),
            ("read", "hello.txt"),
        ];
        let (report, _) = probe("17:41215920:G:T", &options)?;
        assert_eq!(report["read"], "hello");
        Ok(())
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn probe_reads_from_multi_thread_runtime() -> Result<(), anyhow::Error> {
        probe_reads_from_async_code().await
    }

    #[tokio::test(flavor = "current_thread")]
    async fn probe_reads_from_current_thread_runtime() -> Result<(), anyhow::Error> {
        probe_reads_from_async_code().await
    }

    #[rstest::rstest]
    #[case(&[("fail-init", "init went wrong")], "init went wrong")]
    #[case(&[("data", "tests/data/plugins/missing")], "data directory")]
    #[case(&[("data", "a"), ("data", "b")], "more than once")]
    fn load_fails(#[case] options: &[(&str, &str)], #[case] message: &str) {
        let error = load(PROBE_WASM, options).unwrap_err().to_string();
        assert!(error.contains(message), "{error}");
    }

    #[rstest::rstest]
    #[case(&[("fail", "annotate went wrong")], "annotate went wrong")]
    #[case(&[("bad-index", "")], "out of range")]
    fn annotate_fails(#[case] options: &[(&str, &str)], #[case] message: &str) {
        let error = probe("17:41215920:G:T", options).unwrap_err().to_string();
        assert!(error.contains(message), "{error}");
    }
}
