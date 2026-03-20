use mehari::annotate::seqvars::VariantAnnotation;
use mehari::annotate::seqvars::csq::{Config, ConsequencePredictor, VcfVariant};
use mehari::annotate::seqvars::load_tx_db;
use mehari::annotate::seqvars::provider::{
    ConfigBuilder as ProviderConfigBuilder, Provider as MehariProvider,
};
use mehari::db::merge::merge_transcript_databases;
use pyo3::prelude::*;
use pythonize::pythonize;
use std::path::PathBuf;
use std::sync::Arc;

#[pyclass(name = "SeqvarsAnnotator")]
pub struct PySeqvarsAnnotator {
    predictor: ConsequencePredictor,
}

#[pymethods]
impl PySeqvarsAnnotator {
    /// Initialize the annotator with transcript databases and an optional FASTA reference.
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

        let merged_tx_db = merge_transcript_databases(tx_dbs).map_err(|e| {
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

        let annotation = VariantAnnotation {
            consequences: ann_fields_opt.unwrap_or_default(),
            frequencies: None,
            clinvar: None,
        };

        let py_dict = pythonize(py, &annotation).map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("Serialization error: {}", e))
        })?;

        Ok(py_dict)
    }
}

/// The Python module definition
#[pymodule]
fn mehari_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySeqvarsAnnotator>()?;
    Ok(())
}
