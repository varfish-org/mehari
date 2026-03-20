use arrow::array::{Array, Int32Array, RecordBatch, StringArray};
use arrow::compute::cast;
use arrow::datatypes::{DataType, FieldRef};
use arrow::pyarrow::{FromPyArrow, ToPyArrow};
use mehari::annotate::seqvars::ann::{AnnField, Consequence, Pos, PutativeImpact, Rank};
use mehari::annotate::seqvars::csq::{Config, ConsequencePredictor, VcfVariant};
use mehari::annotate::seqvars::load_tx_db;
use mehari::annotate::seqvars::provider::{
    ConfigBuilder as ProviderConfigBuilder, Provider as MehariProvider,
};
use pyo3::prelude::*;
use pythonize::pythonize;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_arrow::schema::{SchemaLike, TracingOptions};
use std::path::PathBuf;
use std::sync::Arc;
use strum::IntoEnumIterator;

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
}

impl From<AnnField> for ArrowAnnField {
    fn from(f: AnnField) -> Self {
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
        }
    }
}

#[derive(Deserialize, Serialize)]
struct ArrowResult {
    pub chromosome: String,
    pub position: i32,
    pub reference: String,
    pub alternative: String,
    pub consequences: Vec<ArrowAnnField>,
}

#[pyclass(name = "SeqvarsAnnotator")]
pub struct PySeqvarsAnnotator {
    predictor: ConsequencePredictor,
}

#[pymethods]
impl PySeqvarsAnnotator {
    #[new]
    #[pyo3(signature = (transcript_db_paths, reference_path=None))]
    fn new(transcript_db_paths: Vec<String>, reference_path: Option<String>) -> PyResult<Self> {
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

        let merged_tx_db = mehari::db::merge::merge_transcript_databases(tx_dbs).map_err(|e| {
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

        let predictor = ConsequencePredictor::new(provider, Config::default());

        Ok(Self { predictor })
    }

    /// Annotate a single variant. Returns a Python dictionary.
    #[pyo3(signature = (chrom, position, reference, alternative))]
    fn annotate<'py>(
        &self,
        py: Python<'py>,
        chrom: &str,
        position: i32,
        reference: &str,
        alternative: &str,
    ) -> PyResult<Bound<'py, PyAny>> {
        let variant = VcfVariant {
            chromosome: chrom.to_string(),
            position,
            reference: reference.to_string(),
            alternative: alternative.to_string(),
        };

        let ann_fields_opt = self
            .predictor
            .predict(&variant)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let arrow_anns: Vec<ArrowAnnField> = ann_fields_opt
            .unwrap_or_default()
            .into_iter()
            .map(ArrowAnnField::from)
            .collect();

        #[derive(Serialize)]
        struct SingleResult {
            consequences: Vec<ArrowAnnField>,
        }

        let py_dict = pythonize(
            py,
            &SingleResult {
                consequences: arrow_anns,
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

        let results: Vec<ArrowResult> = indices
            .par_iter()
            .map(|&i| {
                let variant = VcfVariant {
                    chromosome: chrom_arr.value(i).to_string(),
                    position: pos_arr.value(i),
                    reference: ref_arr.value(i).to_string(),
                    alternative: alt_arr.value(i).to_string(),
                };

                let ann_fields = self
                    .predictor
                    .predict(&variant)
                    .unwrap_or(None)
                    .unwrap_or_default();

                let arrow_anns: Vec<ArrowAnnField> =
                    ann_fields.into_iter().map(ArrowAnnField::from).collect();

                ArrowResult {
                    chromosome: variant.chromosome,
                    position: variant.position,
                    reference: variant.reference,
                    alternative: variant.alternative,
                    consequences: arrow_anns,
                }
            })
            .collect();

        let options = TracingOptions::default()
            .allow_null_fields(true)
            .string_dictionary_encoding(true);

        let fields = Vec::<FieldRef>::from_type::<ArrowResult>(options).map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("Schema error: {}", e))
        })?;

        let out_batch = serde_arrow::to_record_batch(&fields, &results).map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("Serialization error: {}", e))
        })?;

        Ok(out_batch.to_pyarrow(py)?)
    }
}

#[pymodule(name = "_mehari")]
fn mehari_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySeqvarsAnnotator>()?;
    m.add_function(wrap_pyfunction!(consequence_variants, m)?)?;
    m.add_function(wrap_pyfunction!(putative_impact_variants, m)?)?;
    Ok(())
}
