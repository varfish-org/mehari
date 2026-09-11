use arrow::array::{Array, Int32Array, RecordBatch, StringArray};
use arrow::compute::cast;
use arrow::datatypes::{DataType, Field, FieldRef};
use arrow::pyarrow::{FromPyArrow, ToPyArrow};
use mehari::annotate::seqvars::consequence::load_tx_db;
use mehari::annotate::seqvars::consequence::logic::ConsequencePredictor;
use mehari::annotate::seqvars::consequence::terms::{
    ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_TX_SEQ_ALT, ANN_TX_SEQ_REF, AnnField, Consequence,
    FeatureBiotype, Pos, PutativeImpact, Rank,
};
use mehari::annotate::seqvars::consequence::{ConfigBuilder, SequenceReporting, VcfVariant};
use mehari::annotate::seqvars::provider::{
    ConfigBuilder as ProviderConfigBuilder, Provider as MehariProvider,
};
use mehari::common::progress::{Progress, ProgressBar, Unit, hidden_bar};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pythonize::pythonize;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_arrow::schema::{SchemaLike, TracingOptions};
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Duration;
use strum::IntoEnumIterator;

/// Clears pyo3-log's cache of Python log levels. Set once, at module import.
static LOG_RESET_HANDLE: OnceLock<pyo3_log::ResetHandle> = OnceLock::new();

/// Makes pyo3-log reread the Python log levels, so recent `logging` configuration applies.
fn reread_log_levels() {
    if let Some(handle) = LOG_RESET_HANDLE.get() {
        handle.reset();
    }
}

#[pyfunction]
fn consequence_variants() -> Vec<String> {
    Consequence::iter()
        .map(|c: Consequence| c.to_string())
        .collect()
}

#[pyfunction]
fn putative_impact_variants() -> Vec<String> {
    PutativeImpact::iter()
        .map(|i: PutativeImpact| i.to_string())
        .collect()
}

#[pyfunction]
fn feature_biotype_variants() -> Vec<String> {
    FeatureBiotype::iter().map(|b| b.to_string()).collect()
}

// _TraceAnnField is the same as ArrowAnnField but without the custom_fields field,
// to enable static tracing of the schema.
#[derive(Deserialize, Serialize)]
struct _TraceAnnField {
    pub allele: String,
    pub consequences: Vec<String>,
    pub putative_impact: String,
    pub gene_symbol: String,
    pub gene_id: String,
    pub feature_type: String,
    pub feature_id: String,
    pub feature_biotype: Vec<String>,
    pub feature_tags: Vec<String>,
    pub rank: Option<Rank>,
    pub hgvs_g: Option<String>,
    pub hgvs_n: Option<String>,
    pub hgvs_c: Option<String>,
    pub hgvs_p: Option<String>,
    pub cdna_pos: Option<Pos>,
    pub cds_pos: Option<Pos>,
    pub protein_pos: Option<Pos>,
    pub distance: Option<i32>,
    pub strand: i32,
    pub messages: Option<Vec<String>>,
}

#[derive(Deserialize, Serialize)]
struct ArrowAnnField {
    pub allele: String,
    pub consequences: Vec<String>,
    pub putative_impact: String,
    pub gene_symbol: String,
    pub gene_id: String,
    pub feature_type: String,
    pub feature_id: String,
    pub feature_biotype: Vec<String>,
    pub feature_tags: Vec<String>,
    pub rank: Option<Rank>,
    pub hgvs_g: Option<String>,
    pub hgvs_n: Option<String>,
    pub hgvs_c: Option<String>,
    pub hgvs_p: Option<String>,
    pub cdna_pos: Option<Pos>,
    pub cds_pos: Option<Pos>,
    pub protein_pos: Option<Pos>,
    pub distance: Option<i32>,
    pub strand: i32,
    pub messages: Option<Vec<String>>,
    pub custom_fields: Option<std::collections::BTreeMap<String, Option<String>>>,
}

impl ArrowAnnField {
    fn from_ann_field(f: AnnField, custom_columns: &[String]) -> Self {
        let custom_fields = if custom_columns.is_empty() {
            None
        } else {
            let mut map = f.custom_fields;
            for col in custom_columns {
                map.entry(col.clone()).or_insert(None);
            }
            Some(map)
        };

        Self {
            allele: f.allele.to_string(),
            consequences: f.consequences.iter().map(|c| c.to_string()).collect(),
            putative_impact: f.putative_impact.to_string(),
            gene_symbol: f.gene_symbol,
            gene_id: f.gene_id,
            feature_type: f.feature_type.to_string(),
            feature_id: f.feature_id,
            feature_biotype: f.feature_biotype.iter().map(|b| b.to_string()).collect(),
            feature_tags: f.feature_tags.iter().map(|t| t.to_string()).collect(),
            rank: f.rank,
            hgvs_g: f.hgvs_g,
            hgvs_n: f.hgvs_n,
            hgvs_c: f.hgvs_c,
            hgvs_p: f.hgvs_p,
            cdna_pos: f.cdna_pos,
            cds_pos: f.cds_pos,
            protein_pos: f.protein_pos,
            distance: f.distance,
            strand: f.strand,
            messages: f
                .messages
                .map(|msgs| msgs.iter().map(|m| m.to_string()).collect()),
            custom_fields,
        }
    }
}

#[derive(Deserialize, Serialize)]
struct ArrowResult {
    pub annotation: Vec<ArrowAnnField>,
}

#[pyclass(name = "SeqvarsAnnotator")]
pub struct PySeqvarsAnnotator {
    predictor: ConsequencePredictor,
    fields: Vec<FieldRef>,
    custom_columns: Vec<String>,
}

#[pymethods]
impl PySeqvarsAnnotator {
    #[new]
    #[pyo3(signature = (transcript_db_paths, reference_path=None, report_cdna_sequence="none", report_protein_sequence="none"))]
    fn new(
        transcript_db_paths: Vec<String>,
        reference_path: Option<String>,
        report_cdna_sequence: &str,
        report_protein_sequence: &str,
    ) -> PyResult<Self> {
        reread_log_levels();

        let mut tx_dbs = Vec::new();
        for path in transcript_db_paths {
            let db = load_tx_db(&path).map_err(|e| {
                pyo3::exceptions::PyIOError::new_err(format!(
                    "Failed to load tx_db {}: {}",
                    path, e
                ))
            })?;
            tx_dbs.push(db);
        }

        let merged_tx_db = mehari::db::transcripts::merge::merge_transcript_databases(tx_dbs)
            .map_err(|e| {
                pyo3::exceptions::PyValueError::new_err(format!("Failed to merge databases: {}", e))
            })?;

        let provider_config = ProviderConfigBuilder::default()
            .build()
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let provider = Arc::new(MehariProvider::new(
            merged_tx_db,
            reference_path.map(PathBuf::from),
            false,
            provider_config,
        ));

        let report_cdna = SequenceReporting::from_str(report_cdna_sequence).map_err(|_| {
            pyo3::exceptions::PyValueError::new_err("Invalid cdna sequence reporting option")
        })?;
        let report_protein =
            SequenceReporting::from_str(report_protein_sequence).map_err(|_| {
                pyo3::exceptions::PyValueError::new_err("Invalid protein sequence reporting option")
            })?;

        let mut custom_columns = Vec::new();
        if report_cdna.includes_ref() {
            custom_columns.push(ANN_TX_SEQ_REF.to_string());
        }
        if report_cdna.includes_alt() {
            custom_columns.push(ANN_TX_SEQ_ALT.to_string());
        }
        if report_protein.includes_ref() {
            custom_columns.push(ANN_AA_SEQ_REF.to_string());
        }
        if report_protein.includes_alt() {
            custom_columns.push(ANN_AA_SEQ_ALT.to_string());
        }

        let config = ConfigBuilder::default()
            .report_cdna_sequence(report_cdna)
            .report_protein_sequence(report_protein)
            .custom_columns(custom_columns.clone())
            .build()
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let options = TracingOptions::default().allow_null_fields(true);
        let mut ann_field_inner_fields = Vec::<FieldRef>::from_type::<_TraceAnnField>(options)
            .map_err(|e| {
                pyo3::exceptions::PyRuntimeError::new_err(format!("Schema trace error: {}", e))
            })?;
        let custom_fields_data_type = DataType::Struct(
            custom_columns
                .iter()
                .map(|name| Field::new(name, DataType::Utf8, true))
                .collect::<Vec<_>>()
                .into(),
        );
        ann_field_inner_fields.push(Arc::new(Field::new(
            "custom_fields",
            custom_fields_data_type,
            true,
        )));

        let ann_field_struct = DataType::Struct(ann_field_inner_fields.into());
        let annotation_list_field = Field::new("item", ann_field_struct, true);
        let fields = vec![Arc::new(Field::new(
            "annotation",
            DataType::List(Arc::new(annotation_list_field)),
            true,
        ))];

        let predictor = ConsequencePredictor::new(provider, config);

        Ok(Self {
            predictor,
            fields,
            custom_columns,
        })
    }

    /// Annotate a single variant. Returns a Python dictionary.
    #[pyo3(signature = (chromosome, position, reference, alternative))]
    fn annotate<'py>(
        &self,
        py: Python<'py>,
        chromosome: &str,
        position: i32,
        reference: &str,
        alternative: &str,
    ) -> PyResult<Bound<'py, PyAny>> {
        let variant = VcfVariant {
            chromosome: chromosome.to_string(),
            position,
            reference: reference.to_string(),
            alternative: alternative.to_string(),
        };

        let ann_fields_opt = self.predictor.predict(&variant).map_err(|e| match e {
            mehari::errors::SeqvarsError::UnknownChromosomeAccession
            | mehari::errors::SeqvarsError::InvalidCoordinates(_, _) => {
                pyo3::exceptions::PyValueError::new_err(e.to_string())
            }
            _ => pyo3::exceptions::PyRuntimeError::new_err(e.to_string()),
        })?;

        let arrow_anns: Vec<ArrowAnnField> = ann_fields_opt
            .unwrap_or_default()
            .into_iter()
            .map(|f| ArrowAnnField::from_ann_field(f, &self.custom_columns))
            .collect();

        #[derive(Serialize)]
        struct SingleResult {
            annotation: Vec<ArrowAnnField>,
        }

        let py_dict = pythonize(
            py,
            &SingleResult {
                annotation: arrow_anns,
            },
        )
        .map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("Serialization error: {}", e))
        })?;

        Ok(py_dict)
    }

    /// Batch annotation via arrow (e.g., for use with polars).
    /// Expects an RecordBatch with columns: 'chromosome', 'position', 'reference', 'alternative'
    #[pyo3(signature = (batch))]
    fn annotate_batch<'py>(
        &self,
        py: Python<'py>,
        batch: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let record_batch = RecordBatch::from_pyarrow_bound(batch).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("Invalid Arrow batch: {}", e))
        })?;

        let get_string_col = |name: &str| -> PyResult<StringArray> {
            let col = record_batch.column_by_name(name).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!("Missing column '{}'", name))
            })?;

            let cast_col = cast(col, &DataType::Utf8).map_err(|e| {
                pyo3::exceptions::PyTypeError::new_err(format!(
                    "Cannot cast '{}' to string: {}",
                    name, e
                ))
            })?;

            let string_arr = cast_col
                .as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| {
                    pyo3::exceptions::PyTypeError::new_err(format!(
                        "Downcast failed for string column '{}'",
                        name
                    ))
                })?;

            Ok(string_arr.clone())
        };

        let get_i32_col = |name: &str| -> PyResult<Int32Array> {
            let col = record_batch.column_by_name(name).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!("Missing column '{}'", name))
            })?;

            let cast_col = cast(col, &DataType::Int32).map_err(|e| {
                pyo3::exceptions::PyTypeError::new_err(format!(
                    "Cannot cast '{}' to Int32: {}",
                    name, e
                ))
            })?;

            let int_arr = cast_col
                .as_any()
                .downcast_ref::<Int32Array>()
                .ok_or_else(|| {
                    pyo3::exceptions::PyTypeError::new_err(format!(
                        "Downcast failed for int column '{}'",
                        name
                    ))
                })?;

            Ok(int_arr.clone())
        };

        let chrom_arr = get_string_col("chromosome")?;
        let pos_arr = get_i32_col("position")?;
        let ref_arr = get_string_col("reference")?;
        let alt_arr = get_string_col("alternative")?;

        let num_rows = record_batch.num_rows();
        let indices: Vec<usize> = (0..num_rows).collect();

        // Release the GIL: the rayon workers log through pyo3-log, which needs the GIL.
        let results: Result<Vec<ArrowResult>, anyhow::Error> = py.detach(|| {
            indices
                .par_iter()
                .map(|&i| {
                    let variant = VcfVariant {
                        chromosome: chrom_arr.value(i).to_string(),
                        position: pos_arr.value(i),
                        reference: ref_arr.value(i).to_string(),
                        alternative: alt_arr.value(i).to_string(),
                    };

                    let ann_fields = self.predictor.predict(&variant)?.unwrap_or_default();

                    let arrow_anns: Vec<ArrowAnnField> = ann_fields
                        .into_iter()
                        .map(|f| ArrowAnnField::from_ann_field(f, &self.custom_columns))
                        .collect();

                    Ok(ArrowResult {
                        annotation: arrow_anns,
                    })
                })
                .collect()
        });

        let results = results.map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!(
                "Prediction failed during batch processing: {}",
                e
            ))
        })?;

        let out_batch = serde_arrow::to_record_batch(&self.fields, &results).map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("Serialization error: {}", e))
        })?;

        out_batch.to_pyarrow(py)
    }

    /// Annotate a group of phased variants. Returns a Python dictionary.
    #[pyo3(signature = (variants))]
    fn annotate_multiple<'py>(
        &self,
        py: Python<'py>,
        variants: Vec<(String, i32, String, String)>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let vcf_variants: Vec<VcfVariant> = variants
            .into_iter()
            .map(|(chrom, pos, ref_seq, alt_seq)| VcfVariant {
                chromosome: chrom,
                position: pos,
                reference: ref_seq,
                alternative: alt_seq,
            })
            .collect();

        let ann_fields_opt =
            self.predictor
                .predict_multiple(&vcf_variants)
                .map_err(|e| match e {
                    mehari::errors::SeqvarsError::GroupValidation(_)
                    | mehari::errors::SeqvarsError::UnknownChromosomeAccession
                    | mehari::errors::SeqvarsError::InvalidCoordinates(_, _) => {
                        pyo3::exceptions::PyValueError::new_err(e.to_string())
                    }
                    _ => pyo3::exceptions::PyRuntimeError::new_err(e.to_string()),
                })?;

        let arrow_anns: Vec<ArrowAnnField> = ann_fields_opt
            .unwrap_or_default()
            .into_iter()
            .map(|f| ArrowAnnField::from_ann_field(f, &self.custom_columns))
            .collect();

        #[derive(Serialize)]
        struct SingleResult {
            annotation: Vec<ArrowAnnField>,
        }

        let py_dict = pythonize(
            py,
            &SingleResult {
                annotation: arrow_anns,
            },
        )
        .map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("Serialization error: {}", e))
        })?;

        Ok(py_dict)
    }
}

/// A progress bar that the Rust side created for one step.
#[derive(Clone)]
struct StepBar {
    step: String,
    total: u64,
    unit: Unit,
    bar: ProgressBar,
}

/// Creates hidden progress bars and keeps them, so that Python can show them.
#[derive(Default)]
struct StepBars(Mutex<Vec<StepBar>>);

impl StepBars {
    fn snapshot(&self) -> Vec<StepBar> {
        self.0.lock().map(|bars| bars.clone()).unwrap_or_default()
    }
}

impl Progress for StepBars {
    fn bar(&self, step: &str, total: u64, unit: Unit) -> ProgressBar {
        let bar = hidden_bar(total);
        if let Ok(mut bars) = self.0.lock() {
            bars.push(StepBar {
                step: step.to_string(),
                total,
                unit,
                bar: bar.clone(),
            });
        }
        bar
    }
}

/// Mirrors the Rust progress bars onto tqdm-compatible Python bars, one per step.
///
/// Exceptions raised by a Python bar go to `sys.unraisablehook`, so they cannot abort the operation.
struct TqdmBars {
    /// Creates a Python bar, called like `tqdm(total=..., desc=..., unit=...)`.
    factory: Py<PyAny>,
    /// Per step: the open Python bar, and the position that it shows.
    shown: Vec<(Option<Py<PyAny>>, u64)>,
}

impl TqdmBars {
    /// Updates the Python bars, and closes the finished ones. `done` closes all of them.
    fn sync(&mut self, py: Python<'_>, steps: &[StepBar], done: bool) {
        for (i, step) in steps.iter().enumerate() {
            if i == self.shown.len() {
                let created = self.create(py, step);
                self.shown.push((created, 0));
            }
            let (open, shown_position) = &mut self.shown[i];
            let Some(py_bar) = open.as_ref() else {
                continue;
            };
            let position = step.bar.position();
            let close = done || step.bar.is_finished();
            let n = position.saturating_sub(*shown_position);
            let result = advance_py_bar(py, py_bar, n, close);
            *shown_position = position;
            if let Err(err) = result {
                err.write_unraisable(py, Some(py_bar.bind(py)));
                *open = None;
            } else if close {
                *open = None;
            }
        }
    }

    fn create(&self, py: Python<'_>, step: &StepBar) -> Option<Py<PyAny>> {
        let (unit, unit_scale, unit_divisor) = match step.unit {
            Unit::Bytes => ("B", true, 1024),
            Unit::Transcripts => ("tx", false, 1000),
        };
        let kwargs = PyDict::new(py);
        let created = kwargs
            .set_item("total", step.total)
            .and_then(|_| kwargs.set_item("desc", &step.step))
            .and_then(|_| kwargs.set_item("unit", unit))
            .and_then(|_| kwargs.set_item("unit_scale", unit_scale))
            .and_then(|_| kwargs.set_item("unit_divisor", unit_divisor))
            .and_then(|_| self.factory.call(py, (), Some(&kwargs)));
        match created {
            Ok(py_bar) => Some(py_bar),
            Err(err) => {
                err.write_unraisable(py, Some(self.factory.bind(py)));
                None
            }
        }
    }
}

fn advance_py_bar(py: Python<'_>, py_bar: &Py<PyAny>, n: u64, close: bool) -> PyResult<()> {
    if n > 0 {
        py_bar.call_method1(py, "update", (n,))?;
    }
    if close {
        py_bar.call_method0(py, "close")?;
    }
    Ok(())
}

/// Builds the database on a worker thread and shows its progress on bars from `factory`.
///
/// Only the calling thread touches the Python bars, ten times per second. So the Rust threads
/// never wait for the GIL to report progress.
fn run_with_tqdm(
    py: Python<'_>,
    factory: Py<PyAny>,
    common_args: &mehari::common::Args,
    args: &mehari::db::transcripts::create::cli::Args,
) -> anyhow::Result<()> {
    let steps = StepBars::default();
    let mut tqdm = TqdmBars {
        factory,
        shown: Vec::new(),
    };
    py.detach(|| {
        std::thread::scope(|scope| {
            let worker = scope.spawn(|| {
                mehari::db::transcripts::create::run_with_progress(common_args, args, &steps)
            });
            while !worker.is_finished() {
                std::thread::sleep(Duration::from_millis(100));
                Python::attach(|py| tqdm.sync(py, &steps.snapshot(), false));
            }
            Python::attach(|py| tqdm.sync(py, &steps.snapshot(), true));
            worker
                .join()
                .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
        })
    })
}

#[pyfunction]
#[pyo3(signature = (
    assembly,
    annotation,
    output,
    transcript_source,
    assembly_version=None,
    annotation_version=None,
    transcript_source_version=None,
    seqrepo=None,
    transcript_sequences=None,
    mane_transcripts=None,
    disable_filters=false,
    threads=1,
    compression_level=19,
    progress=None
))]
#[allow(clippy::too_many_arguments)]
fn build_transcript_db(
    py: Python<'_>,
    assembly: String,
    annotation: Vec<PathBuf>,
    output: PathBuf,
    transcript_source: String,
    assembly_version: Option<String>,
    annotation_version: Option<String>,
    transcript_source_version: Option<String>,
    seqrepo: Option<PathBuf>,
    transcript_sequences: Option<PathBuf>,
    mane_transcripts: Option<PathBuf>,
    disable_filters: bool,
    threads: usize,
    compression_level: i32,
    progress: Option<Py<PyAny>>,
) -> PyResult<()> {
    if seqrepo.is_none() && transcript_sequences.is_none() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Either 'seqrepo' or 'transcript_sequences' must be provided",
        ));
    }

    if transcript_source.to_lowercase() == "ensembl" && transcript_source_version.is_none() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "'transcript_source_version' is required when source is 'ensembl'",
        ));
    }

    let args = mehari::db::transcripts::create::cli::Args {
        assembly,
        assembly_version,
        annotation,
        annotation_version,
        transcript_source,
        transcript_source_version,
        seqrepo,
        transcript_sequences,
        mane_transcripts,
        disable_filters,
        threads,
        compression_level,
        output,
    };

    let common_args = mehari::common::Args::default();

    reread_log_levels();
    let result = match progress {
        // Release the GIL: `run` logs from a rayon thread pool, and pyo3-log needs the GIL.
        None => py.detach(|| mehari::db::transcripts::create::run(&common_args, &args)),
        Some(factory) => run_with_tqdm(py, factory, &common_args, &args),
    };
    result.map_err(|e| {
        pyo3::exceptions::PyRuntimeError::new_err(format!(
            "Failed to build transcript database: {}",
            e
        ))
    })?;

    Ok(())
}

#[pymodule(name = "_mehari")]
fn mehari_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Without a tracing subscriber, mehari's `tracing` events become `log` records
    // (tracing's "log" feature). pyo3-log passes these on to Python's `logging`.
    // If another logger is already set, import the module without the logging bridge.
    if let Ok(log_reset_handle) = pyo3_log::try_init() {
        let _ = LOG_RESET_HANDLE.set(log_reset_handle);
    }

    m.add_class::<PySeqvarsAnnotator>()?;
    m.add_function(wrap_pyfunction!(consequence_variants, m)?)?;
    m.add_function(wrap_pyfunction!(putative_impact_variants, m)?)?;
    m.add_function(wrap_pyfunction!(feature_biotype_variants, m)?)?;
    m.add_function(wrap_pyfunction!(build_transcript_db, m)?)?;
    Ok(())
}
