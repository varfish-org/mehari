//! Plugins compiled to WebAssembly components, implementing the `plugin` world of
//! `wit/plugin.wit`.

use std::cell::RefCell;
use std::path::{Path, PathBuf};

use hgvs::sequences::TranslationTable;
use itertools::Itertools;
use thread_local::ThreadLocal;
use wasmtime::component::{Component, HasSelf, Linker, Resource, ResourceTable};
use wasmtime::{Config, Engine, Store, bail, ensure, format_err};
use wasmtime_wasi::{FsPerms, WasiCtx, WasiCtxBuilder, WasiCtxView, WasiView};

use super::{RunInfo, TranscriptPlugin, TxAnnotation, TxContext};
use crate::annotate::seqvars::consequence::VcfVariant;
use crate::annotate::seqvars::consequence::terms::{AnnField, Pos};
use crate::errors::SeqvarsError;
use crate::pbs::txs::{GenomeAlignment, Strand};

mod bindings {
    wasmtime::component::bindgen!({
        path: "wit",
        world: "plugin",
        imports: { default: trappable },
        with: {
            "mehari:plugin/types.setup": super::Setup,
            "mehari:plugin/types.prediction": super::Prediction,
            "mehari:plugin/types.annotation": super::Annotation,
        },
    });
}

use bindings::Plugin;
use bindings::mehari::plugin::{host, types};

/// The option naming a directory that the plugin may read, at [`DATA_DIR`].
pub const DATA_OPTION: &str = "data";
/// Where plugins find the directory given with [`DATA_OPTION`].
pub const DATA_DIR: &str = "/data";
/// Option names starting with this are reserved for mehari.
pub const RESERVED_OPTION_PREFIX: &str = "mehari-";

/// A plugin loaded from a WebAssembly component.
///
/// Plugins get no network or environment access, and read access only to the directory given
/// with the `data` option; stderr is passed through so their errors show.
///
/// Experimental, see the [module](super) documentation.
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
    declaration: Declaration,
    /// One instance per thread, created by its first call; `None` after an instance trapped.
    #[derivative(Debug = "ignore")]
    instances: ThreadLocal<RefCell<Option<Instance>>>,
}

type Instance = (Store<HostState>, Plugin);

/// A column a plugin declared, and the type of its values.
type Column = (String, ColumnType);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ColumnType {
    value: ValueType,
    /// Holds a list of values, written joined by `&`.
    list: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, strum::Display)]
#[strum(serialize_all = "lowercase")]
enum ValueType {
    String,
    Integer,
    Float,
    Boolean,
}

impl ColumnType {
    const fn single(value: ValueType) -> Self {
        Self { value, list: false }
    }

    const fn list(value: ValueType) -> Self {
        Self { value, list: true }
    }
}

impl std::fmt::Display for ColumnType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.list {
            write!(f, "{} list", self.value)
        } else {
            write!(f, "{}", self.value)
        }
    }
}

/// What `init` declared; every instance must declare the same.
#[derive(Debug, Clone, Default, PartialEq)]
struct Declaration {
    columns: Vec<Column>,
    /// Call `annotate` also for predictions without annotations.
    annotate_empty: bool,
}

/// The options and the declaration of one `init` call.
pub struct Setup {
    options: Vec<(String, String)>,
    declaration: Declaration,
}

/// The prediction of an `annotate` call, by the number of the call.
pub struct Prediction {
    call: u64,
}

/// An annotation of an `annotate` call, by the number of the call and its index.
pub struct Annotation {
    call: u64,
    index: usize,
}

/// The input and output of the current `annotate` call, lent to the plugin.
struct Batch {
    call: u64,
    variants: Vec<VcfVariant>,
    txs: Vec<TxAnnotation>,
}

struct HostState {
    table: ResourceTable,
    wasi: WasiCtx,
    run: RunInfo,
    /// As declared by `init`.
    columns: Vec<Column>,
    /// Number of `annotate` calls so far.
    calls: u64,
    batch: Option<Batch>,
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
        let mut data = None;
        let mut plugin_options = Vec::new();
        for (key, value) in options {
            if key == DATA_OPTION {
                if data.replace(PathBuf::from(value)).is_some() {
                    return Err(plugin_message(path, "option `data` given more than once"));
                }
            } else if key.starts_with(RESERVED_OPTION_PREFIX) {
                return Err(plugin_message(
                    path,
                    &format!("option `{key}` is reserved for mehari"),
                ));
            } else {
                plugin_options.push((key.clone(), value.clone()));
            }
        }

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
            options: plugin_options,
            data,
            run,
            engine,
            pre,
            declaration: Declaration::default(),
            instances: ThreadLocal::new(),
        };
        let (instance, declaration) = outside_runtime(path, || plugin.instantiate())?;
        plugin.declaration = declaration;
        plugin.instances.get_or_default().replace(Some(instance));
        Ok(plugin)
    }

    /// A new instance and what its `init` declared.
    fn instantiate(&self) -> Result<(Instance, Declaration), SeqvarsError> {
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
            columns: Vec::new(),
            calls: 0,
            batch: None,
        };
        let mut store = Store::new(&self.engine, state);
        let plugin = self.pre.instantiate(&mut store).map_err(error)?;

        let setup = Setup {
            options: self.options.clone(),
            declaration: Declaration::default(),
        };
        let handle = store
            .data_mut()
            .table
            .push(setup)
            .map_err(|e| error(e.into()))?;
        let result = plugin
            .mehari_plugin_transcript_plugin()
            .call_init(&mut store, Resource::new_borrow(handle.rep()));
        let Setup { declaration, .. } = store
            .data_mut()
            .table
            .delete(handle)
            .map_err(|e| error(e.into()))?;
        result
            .map_err(error)?
            .map_err(|message| plugin_message(&self.path, &message))?;

        store.data_mut().columns = declaration.columns.clone();
        Ok(((store, plugin), declaration))
    }

    /// A new instance for the calling thread.
    fn new_instance(&self) -> Result<Instance, SeqvarsError> {
        let (instance, declaration) = outside_runtime(&self.path, || self.instantiate())?;
        if declaration != self.declaration {
            return Err(plugin_message(
                &self.path,
                "init declared something else than before",
            ));
        }
        Ok(instance)
    }

    /// Lend the annotations to the plugin for the call and take them back afterwards.
    fn call(
        &self,
        (store, plugin): &mut Instance,
        variants: &[VcfVariant],
        txs: &mut Vec<TxAnnotation>,
    ) -> Result<(), SeqvarsError> {
        let error = |e: wasmtime::Error| plugin_error(&self.path, e);
        let state = store.data_mut();
        state.calls += 1;
        let call = state.calls;
        let handle = state
            .table
            .push(Prediction { call })
            .map_err(|e| error(e.into()))?;
        state.batch = Some(Batch {
            call,
            variants: variants.to_vec(),
            txs: std::mem::take(txs),
        });
        let result = plugin
            .mehari_plugin_transcript_plugin()
            .call_annotate(&mut *store, Resource::new_borrow(handle.rep()));
        let state = store.data_mut();
        if let Some(batch) = state.batch.take() {
            *txs = batch.txs;
        }
        state.table.delete(handle).map_err(|e| error(e.into()))?;

        result
            .map_err(error)?
            .map_err(|message| plugin_message(&self.path, &message))
    }
}

impl TranscriptPlugin for WasmPlugin {
    fn columns(&self) -> Vec<String> {
        let columns = &self.declaration.columns;
        columns.iter().map(|(name, _)| name.clone()).collect()
    }

    fn annotate(
        &self,
        variants: &[VcfVariant],
        txs: &mut Vec<TxAnnotation>,
    ) -> Result<(), SeqvarsError> {
        if txs.is_empty() && !self.declaration.annotate_empty {
            return Ok(());
        }
        let slot = self.instances.get_or_default();
        let mut instance = match slot.take() {
            Some(instance) => instance,
            None => self.new_instance()?,
        };
        // An instance that trapped is not put back.
        outside_runtime(&self.path, || self.call(&mut instance, variants, txs))?;
        slot.replace(Some(instance));
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
    fn add_column(
        &mut self,
        setup: &Resource<Setup>,
        name: String,
        kind: ColumnType,
    ) -> wasmtime::Result<()> {
        let valid = |c: char| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.');
        ensure!(
            !name.is_empty() && name.chars().all(valid),
            "column name {name:?} has characters other than letters, digits, `_`, `-` and `.`"
        );
        let builtin = AnnField::ann_field_names(&[]);
        ensure!(
            !builtin
                .iter()
                .any(|field| field.eq_ignore_ascii_case(&name)),
            "column name {name:?} is taken by a field of ANN"
        );
        let columns = &mut self.table.get_mut(setup)?.declaration.columns;
        ensure!(
            columns.iter().all(|(other, _)| *other != name),
            "column {name:?} declared twice"
        );
        columns.push((name, kind));
        Ok(())
    }

    fn batch(&self, call: u64) -> wasmtime::Result<&Batch> {
        let batch = self.batch.as_ref().filter(|batch| batch.call == call);
        batch.ok_or_else(|| format_err!("annotation used outside the annotate call it came from"))
    }

    fn tx(&self, annotation: &Resource<Annotation>) -> wasmtime::Result<&TxAnnotation> {
        let &Annotation { call, index } = self.table.get(annotation)?;
        let txs = &self.batch(call)?.txs;
        txs.get(index)
            .ok_or_else(|| format_err!("annotation {index} out of range"))
    }

    fn ann(&self, annotation: &Resource<Annotation>) -> wasmtime::Result<&AnnField> {
        Ok(&self.tx(annotation)?.ann)
    }

    fn context(&self, annotation: &Resource<Annotation>) -> wasmtime::Result<&TxContext> {
        Ok(&self.tx(annotation)?.context)
    }

    fn alignment(
        &self,
        annotation: &Resource<Annotation>,
    ) -> wasmtime::Result<Option<&GenomeAlignment>> {
        let tx = self.context(annotation)?.tx();
        Ok(tx.and_then(|tx| tx.genome_alignments.first()))
    }

    fn set(
        &mut self,
        annotation: &Resource<Annotation>,
        column: String,
        kind: ColumnType,
        value: Option<String>,
    ) -> wasmtime::Result<()> {
        match self.columns.iter().find(|(name, _)| *name == column) {
            None => bail!("column {column:?} was not declared in init"),
            Some((_, declared)) if *declared != kind => {
                bail!("column {column:?} holds {declared} values, not {kind} values")
            }
            Some(_) => {}
        }
        let &Annotation { call, index } = self.table.get(annotation)?;
        let batch = self.batch.as_mut().filter(|batch| batch.call == call);
        let batch = batch
            .ok_or_else(|| format_err!("annotation used outside the annotate call it came from"))?;
        let tx = batch
            .txs
            .get_mut(index)
            .ok_or_else(|| format_err!("annotation {index} out of range"))?;
        tx.ann.custom_fields.insert(column, value);
        Ok(())
    }
}

impl types::HostSetup for HostState {
    fn options(&mut self, setup: Resource<Setup>) -> wasmtime::Result<Vec<(String, String)>> {
        Ok(self.table.get(&setup)?.options.clone())
    }

    fn add_string_column(&mut self, setup: Resource<Setup>, name: String) -> wasmtime::Result<()> {
        self.add_column(&setup, name, ColumnType::single(ValueType::String))
    }

    fn add_integer_column(&mut self, setup: Resource<Setup>, name: String) -> wasmtime::Result<()> {
        self.add_column(&setup, name, ColumnType::single(ValueType::Integer))
    }

    fn add_float_column(&mut self, setup: Resource<Setup>, name: String) -> wasmtime::Result<()> {
        self.add_column(&setup, name, ColumnType::single(ValueType::Float))
    }

    fn add_boolean_column(&mut self, setup: Resource<Setup>, name: String) -> wasmtime::Result<()> {
        self.add_column(&setup, name, ColumnType::single(ValueType::Boolean))
    }

    fn add_string_list_column(
        &mut self,
        setup: Resource<Setup>,
        name: String,
    ) -> wasmtime::Result<()> {
        self.add_column(&setup, name, ColumnType::list(ValueType::String))
    }

    fn add_integer_list_column(
        &mut self,
        setup: Resource<Setup>,
        name: String,
    ) -> wasmtime::Result<()> {
        self.add_column(&setup, name, ColumnType::list(ValueType::Integer))
    }

    fn add_float_list_column(
        &mut self,
        setup: Resource<Setup>,
        name: String,
    ) -> wasmtime::Result<()> {
        self.add_column(&setup, name, ColumnType::list(ValueType::Float))
    }

    fn add_boolean_list_column(
        &mut self,
        setup: Resource<Setup>,
        name: String,
    ) -> wasmtime::Result<()> {
        self.add_column(&setup, name, ColumnType::list(ValueType::Boolean))
    }

    fn annotate_empty_predictions(&mut self, setup: Resource<Setup>) -> wasmtime::Result<()> {
        self.table.get_mut(&setup)?.declaration.annotate_empty = true;
        Ok(())
    }

    fn drop(&mut self, setup: Resource<Setup>) -> wasmtime::Result<()> {
        self.table.delete(setup)?;
        Ok(())
    }
}

impl types::HostPrediction for HostState {
    fn variants(&mut self, p: Resource<Prediction>) -> wasmtime::Result<Vec<types::VcfVariant>> {
        let call = self.table.get(&p)?.call;
        Ok(self
            .batch(call)?
            .variants
            .iter()
            .map(|var| types::VcfVariant {
                chromosome: var.chromosome.clone(),
                position: var.position,
                reference: var.reference.clone(),
                alternative: var.alternative.clone(),
            })
            .collect())
    }

    fn annotations(
        &mut self,
        p: Resource<Prediction>,
    ) -> wasmtime::Result<Vec<Resource<Annotation>>> {
        let call = self.table.get(&p)?.call;
        let count = self.batch(call)?.txs.len();
        (0..count)
            .map(|index| Ok(self.table.push(Annotation { call, index })?))
            .collect()
    }

    fn drop(&mut self, prediction: Resource<Prediction>) -> wasmtime::Result<()> {
        self.table.delete(prediction)?;
        Ok(())
    }
}

impl types::HostAnnotation for HostState {
    fn feature_id(&mut self, a: Resource<Annotation>) -> wasmtime::Result<String> {
        Ok(self.ann(&a)?.feature_id.clone())
    }

    fn feature_type(&mut self, a: Resource<Annotation>) -> wasmtime::Result<String> {
        Ok(self.ann(&a)?.feature_type.to_string())
    }

    fn consequences(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Vec<String>> {
        Ok(strings(&self.ann(&a)?.consequences))
    }

    fn impact(&mut self, a: Resource<Annotation>) -> wasmtime::Result<String> {
        Ok(self.ann(&a)?.putative_impact.to_string())
    }

    fn gene_symbol(&mut self, a: Resource<Annotation>) -> wasmtime::Result<String> {
        Ok(self.ann(&a)?.gene_symbol.clone())
    }

    fn gene_id(&mut self, a: Resource<Annotation>) -> wasmtime::Result<String> {
        Ok(self.ann(&a)?.gene_id.clone())
    }

    fn biotypes(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Vec<String>> {
        Ok(strings(&self.ann(&a)?.feature_biotype))
    }

    fn tags(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Vec<String>> {
        Ok(strings(&self.ann(&a)?.feature_tags))
    }

    fn rank(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<types::Rank>> {
        Ok(self.ann(&a)?.rank.as_ref().map(|rank| types::Rank {
            ord: rank.ord,
            total: rank.total,
        }))
    }

    fn hgvs_g(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<String>> {
        Ok(self.ann(&a)?.hgvs_g.clone())
    }

    fn hgvs_n(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<String>> {
        Ok(self.ann(&a)?.hgvs_n.clone())
    }

    fn hgvs_c(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<String>> {
        Ok(self.ann(&a)?.hgvs_c.clone())
    }

    fn hgvs_p(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<String>> {
        Ok(self.ann(&a)?.hgvs_p.clone())
    }

    fn cdna_pos(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<types::Position>> {
        Ok(self.ann(&a)?.cdna_pos.as_ref().map(position))
    }

    fn cdna_end(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<types::Position>> {
        let tx = self.tx(&a)?;
        Ok(tx.ann.cdna_pos.as_ref().and_then(|start| {
            let (_, end) = tx.context.tx_span()?;
            Some(types::Position {
                ord: end,
                total: start.total,
            })
        }))
    }

    fn cds_pos(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<types::Position>> {
        Ok(self.ann(&a)?.cds_pos.as_ref().map(position))
    }

    fn protein_pos(
        &mut self,
        a: Resource<Annotation>,
    ) -> wasmtime::Result<Option<types::Position>> {
        Ok(self.ann(&a)?.protein_pos.as_ref().map(position))
    }

    fn distance(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<i32>> {
        Ok(self.ann(&a)?.distance)
    }

    fn warnings(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Vec<String>> {
        let messages = &self.ann(&a)?.messages;
        Ok(messages.iter().flatten().map(ToString::to_string).collect())
    }

    fn field(
        &mut self,
        a: Resource<Annotation>,
        column: String,
    ) -> wasmtime::Result<Option<String>> {
        let fields = &self.ann(&a)?.custom_fields;
        Ok(fields.get(&column).cloned().flatten())
    }

    fn protein_id(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<String>> {
        let tx = self.context(&a)?.tx();
        Ok(tx.and_then(|tx| tx.protein.clone()))
    }

    fn contig(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<String>> {
        Ok(self.alignment(&a)?.map(|aln| aln.contig.clone()))
    }

    fn genome_build(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<String>> {
        Ok(self.alignment(&a)?.map(|aln| aln.genome_build.clone()))
    }

    fn strand(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<types::Strand>> {
        Ok(self.alignment(&a)?.map(|aln| {
            if aln.strand == Strand::Minus as i32 {
                types::Strand::Minus
            } else {
                types::Strand::Plus
            }
        }))
    }

    fn exons(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Vec<types::Exon>> {
        let Some(aln) = self.alignment(&a)? else {
            return Ok(Vec::new());
        };
        let exons = aln
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
            .collect::<Option<Vec<_>>>();
        Ok(exons.unwrap_or_default())
    }

    fn genome_cds_start(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<i32>> {
        let aln = self.alignment(&a)?;
        Ok(aln.and_then(|aln| aln.cds_start).map(|start| start + 1))
    }

    fn genome_cds_end(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<i32>> {
        Ok(self.alignment(&a)?.and_then(|aln| aln.cds_end))
    }

    fn reference_transcript(
        &mut self,
        a: Resource<Annotation>,
    ) -> wasmtime::Result<Option<String>> {
        let data = self.context(&a)?.ref_data();
        Ok(data.map(|data| data.transcript_sequence.to_string()))
    }

    fn reference_protein(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<String>> {
        let data = self.context(&a)?.ref_data();
        Ok(data.map(|data| data.aa_sequence.to_string()))
    }

    fn reference_cds_start(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<i32>> {
        Ok(self.context(&a)?.ref_data().map(|data| data.cds_start))
    }

    fn reference_cds_end(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<i32>> {
        Ok(self.context(&a)?.ref_data().map(|data| data.cds_stop))
    }

    fn alternative_transcript(
        &mut self,
        a: Resource<Annotation>,
    ) -> wasmtime::Result<Option<String>> {
        let data = self.context(&a)?.alt_data();
        Ok(data.map(|data| data.transcript_sequence.clone()))
    }

    fn alternative_protein(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<String>> {
        let data = self.context(&a)?.alt_data();
        Ok(data.map(|data| data.aa_sequence.clone()))
    }

    fn alternative_cds_start(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<i32>> {
        Ok(self.context(&a)?.alt_data().map(|data| data.cds_start))
    }

    fn translation_table(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<String>> {
        let data = self.context(&a)?.ref_data();
        Ok(data.map(|data| {
            match data.translation_table {
                TranslationTable::Standard => "standard",
                TranslationTable::Selenocysteine => "selenocysteine",
                TranslationTable::VertebrateMitochondrial => "vertebrate-mitochondrial",
            }
            .into()
        }))
    }

    fn translation_exceptions(
        &mut self,
        a: Resource<Annotation>,
    ) -> wasmtime::Result<Vec<types::TranslationException>> {
        let data = self.context(&a)?.ref_data();
        let exceptions = data.iter().flat_map(|data| &data.translation_exceptions);
        Ok(exceptions
            .map(|exception| types::TranslationException {
                position: exception.position,
                amino_acid: exception.amino_acid.to_string(),
            })
            .collect())
    }

    fn frameshift(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<bool>> {
        let data = self.context(&a)?.alt_data();
        Ok(data.map(|data| data.is_frameshift))
    }

    fn first_changed_codon(&mut self, a: Resource<Annotation>) -> wasmtime::Result<Option<i32>> {
        let data = self.context(&a)?.alt_data();
        Ok(data.and_then(|data| data.variant_start_aa))
    }

    fn set_string(
        &mut self,
        a: Resource<Annotation>,
        column: String,
        value: String,
    ) -> wasmtime::Result<()> {
        check_string(&column, &value)?;
        let kind = ColumnType::single(ValueType::String);
        self.set(&a, column, kind, Some(value))
    }

    fn set_integer(
        &mut self,
        a: Resource<Annotation>,
        column: String,
        value: i64,
    ) -> wasmtime::Result<()> {
        let kind = ColumnType::single(ValueType::Integer);
        self.set(&a, column, kind, Some(value.to_string()))
    }

    fn set_float(
        &mut self,
        a: Resource<Annotation>,
        column: String,
        value: f64,
    ) -> wasmtime::Result<()> {
        check_float(&column, value)?;
        let kind = ColumnType::single(ValueType::Float);
        self.set(&a, column, kind, Some(value.to_string()))
    }

    fn set_boolean(
        &mut self,
        a: Resource<Annotation>,
        column: String,
        value: bool,
    ) -> wasmtime::Result<()> {
        let kind = ColumnType::single(ValueType::Boolean);
        self.set(&a, column, kind, Some(value.to_string()))
    }

    fn set_string_list(
        &mut self,
        a: Resource<Annotation>,
        column: String,
        values: Vec<String>,
    ) -> wasmtime::Result<()> {
        for value in &values {
            check_string(&column, value)?;
        }
        let kind = ColumnType::list(ValueType::String);
        self.set(&a, column, kind, join(values))
    }

    fn set_integer_list(
        &mut self,
        a: Resource<Annotation>,
        column: String,
        values: Vec<i64>,
    ) -> wasmtime::Result<()> {
        let kind = ColumnType::list(ValueType::Integer);
        self.set(&a, column, kind, join(values))
    }

    fn set_float_list(
        &mut self,
        a: Resource<Annotation>,
        column: String,
        values: Vec<f64>,
    ) -> wasmtime::Result<()> {
        for value in &values {
            check_float(&column, *value)?;
        }
        let kind = ColumnType::list(ValueType::Float);
        self.set(&a, column, kind, join(values))
    }

    fn set_boolean_list(
        &mut self,
        a: Resource<Annotation>,
        column: String,
        values: Vec<bool>,
    ) -> wasmtime::Result<()> {
        let kind = ColumnType::list(ValueType::Boolean);
        self.set(&a, column, kind, join(values))
    }

    fn drop(&mut self, annotation: Resource<Annotation>) -> wasmtime::Result<()> {
        self.table.delete(annotation)?;
        Ok(())
    }
}

/// Strings must not contain what separates values in VCF and `ANN`.
fn check_string(column: &str, value: &str) -> wasmtime::Result<()> {
    let separator = |c: char| c.is_whitespace() || matches!(c, '|' | ',' | ';' | '=' | '&');
    if let Some(c) = value.chars().find(|c| separator(*c)) {
        bail!("value {value:?} of column {column:?} contains {c:?}, which separates values in VCF");
    }
    Ok(())
}

fn check_float(column: &str, value: f64) -> wasmtime::Result<()> {
    ensure!(
        value.is_finite(),
        "value {value} of column {column:?} is not finite"
    );
    Ok(())
}

/// A list as written, joined by `&`; `None` if empty.
fn join(values: Vec<impl ToString>) -> Option<String> {
    (!values.is_empty()).then(|| values.iter().map(ToString::to_string).join("&"))
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
        spdi.parse().unwrap()
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
            .split('~')
            .map(|pair| {
                let (key, value) = pair.split_once('^').unwrap();
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
        // Indels on both sides of the last exon-exon junction.
        let group = ["17:41199711:GG:G", "17:41197803:A:ATT"].map(variant);

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

        let nmd_group = |predictor: &ConsequencePredictor| {
            predictor
                .predict_multiple(&group)
                .unwrap()
                .unwrap()
                .into_iter()
                .map(|ann| (ann.feature_id, ann.custom_fields.get(ANN_NMD).cloned()))
                .collect_vec()
        };
        assert_eq!(nmd_group(&wasm), nmd_group(&builtin));
        Ok(())
    }

    /// Every function of the interface, for the missense variant c.5123C>A (p.Ala1708Glu).
    #[test]
    fn probe_reports_everything() -> Result<(), anyhow::Error> {
        let (report, ann) = probe("17:41215920:G:T", &[])?;
        let join = |items: Vec<String>| items.join("+");
        let show = |value: Option<String>| value.unwrap_or_default();
        let expected = [
            ("assembly", "grch37".to_string()),
            ("version", crate::built_info::PKG_VERSION.to_string()),
            ("vep", "true".into()),
            ("variants", "1".into()),
            ("variant-position", "41215920".into()),
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
            ("cdna-end", "5236/7088".into()),
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
            ("protein-id", "NP_009225.1".into()),
            ("contig", "NC_000017.10".into()),
            ("genome-build", "grch37".into()),
            ("strand", "Strand::Minus".into()),
            ("exons", "23".into()),
            ("first-exon", "1-94:41277288-41277381:94M".into()),
            // c.1 is n.114 and chr17:41276113, the stop codon ends at c.5592.
            ("genome-cds", "41197695-41276113".into()),
            ("reference-transcript", "7091".into()),
            ("reference-protein", "1864".into()),
            ("reference-cds", "114-5705".into()),
            ("alternative-transcript", "7091".into()),
            ("alternative-protein", "1864".into()),
            ("alternative-cds-start", "114".into()),
            ("translation-table", "standard".into()),
            ("translation-exceptions", String::new()),
            ("frameshift", "false".into()),
            ("first-changed-codon", "1708".into()),
        ]
        .map(|(key, value)| (key.to_string(), value))
        .into_iter()
        .collect::<BTreeMap<_, _>>();

        assert_eq!(report, expected);
        assert_eq!(report["hgvs-p"], "p.Ala1708Glu");
        Ok(())
    }

    /// Values of each column type, as written into `ANN`.
    #[test]
    fn probe_fills_typed_columns() -> Result<(), anyhow::Error> {
        let (_, ann) = probe("17:41215920:G:T", &[])?;
        let value = |column: &str| ann.custom_fields[column].clone();
        let columns = load(PROBE_WASM, &[])?.columns();
        assert_eq!(
            columns,
            [
                "probe",
                "probe_int",
                "probe_float",
                "probe_bool",
                "probe_strings",
                "probe_ints",
                "probe_floats",
                "probe_bools"
            ]
        );
        // NM_007294.4 sorts first among the BRCA1 transcripts.
        assert_eq!(value("probe_int").as_deref(), Some("-1"));
        assert_eq!(value("probe_float").as_deref(), Some("0.25"));
        assert_eq!(value("probe_bool").as_deref(), Some("true"));
        assert_eq!(value("probe_strings").as_deref(), Some("a&b"));
        assert_eq!(value("probe_ints").as_deref(), Some("1&-2&3"));
        assert_eq!(value("probe_floats").as_deref(), Some("0.5&1.5"));
        assert_eq!(value("probe_bools").as_deref(), Some("true&false"));
        Ok(())
    }

    #[test]
    fn probe_empty_list_leaves_field_empty() -> Result<(), anyhow::Error> {
        let (_, ann) = probe("17:41215920:G:T", &[("empty-list", "")])?;
        assert_eq!(ann.custom_fields["probe_strings"], None);
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
    #[case(&[("mehari-limit", "1")], "reserved for mehari")]
    #[case(&[("column", "a|b")], "characters other than")]
    #[case(&[("column", "probe")], "declared twice")]
    #[case(&[("column", "HGVS.c")], "taken by a field of ANN")]
    #[case(&[("column", "allele")], "taken by a field of ANN")]
    fn load_fails(#[case] options: &[(&str, &str)], #[case] message: &str) {
        let error = load(PROBE_WASM, options).unwrap_err().to_string();
        assert!(error.contains(message), "{error}");
    }

    /// An annotation kept from one call cannot be used in the next.
    #[test]
    fn stale_annotation_fails() -> Result<(), anyhow::Error> {
        let predictor = predictor(vec![Arc::new(load(PROBE_WASM, &[("stale", "")])?)])?;
        predictor.predict(&variant("17:41215920:G:T"))?;
        let error = predictor
            .predict(&variant("17:41276049:A:T"))
            .unwrap_err()
            .to_string();
        assert!(error.contains("outside the annotate call"), "{error}");
        Ok(())
    }

    /// Predictions without annotations reach a plugin only on request.
    #[rstest::rstest]
    #[case(&[("fail", "called")], false)]
    #[case(&[("fail", "called"), ("annotate-empty", "")], true)]
    fn empty_predictions(
        #[case] options: &[(&str, &str)],
        #[case] called: bool,
    ) -> Result<(), anyhow::Error> {
        let predictor = predictor(vec![Arc::new(load(PROBE_WASM, options)?)])?;
        // Far from the genes of the test database.
        let result = predictor.predict(&variant("17:1000000:A:T"));
        assert!(
            matches!(result, Ok(Some(ref anns)) if anns.is_empty()) != called,
            "{result:?}"
        );
        Ok(())
    }

    #[rstest::rstest]
    #[case(&[("fail", "annotate went wrong")], "annotate went wrong")]
    #[case(&[("set", "separator")], "contains '|'")]
    #[case(&[("set", "ampersand")], "contains '&'")]
    #[case(&[("set", "list-element")], "contains '|'")]
    #[case(&[("set", "list-as-single")], "holds string list values")]
    #[case(&[("set", "nan-list")], "not finite")]
    #[case(&[("set", "undeclared")], "not declared")]
    #[case(&[("set", "wrong-type")], "holds string values")]
    #[case(&[("set", "nan")], "not finite")]
    fn annotate_fails(#[case] options: &[(&str, &str)], #[case] message: &str) {
        let error = probe("17:41215920:G:T", options).unwrap_err().to_string();
        assert!(error.contains(message), "{error}");
    }
}
