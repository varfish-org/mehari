//! Annotation of sequence variants.
use std::env;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::pin::Pin;
use std::sync::Arc;
use std::time::Instant;

use self::consequence::terms::AnnField;
use crate::annotate::cli::{PredictorSettings, Sources};
use crate::annotate::seqvars::consequence::logic::ConsequencePredictor;
use crate::common::contig::ContigManager;
use crate::common::noodles::{NoodlesVariantReader, open_variant_reader, open_variant_writer};
use crate::db::keys::Var;
use crate::db::transcripts::merge::merge_transcript_databases;
use crate::errors::SeqvarsError;
use crate::pbs;
use crate::pbs::txs::TxSeqDatabase;
use anyhow::Error;
use biocommons_bioutils::assemblies::Assembly;
use clap::Parser;
use clinvar::{ClinvarAnnotator, ClinvarResult};
use consequence::{ConsequenceAnnotator, VcfVariant};
use frequency::{FreqResult, FrequencyAnnotator};
use indexmap::IndexMap;
use itertools::Itertools;
use noodles::vcf::header::FileFormat;
use noodles::vcf::header::record::value::map::format::Number as FormatNumber;
use noodles::vcf::header::record::value::map::format::Type as FormatType;
use noodles::vcf::header::record::value::map::info::Number;
use noodles::vcf::header::record::value::map::{Info, info::Type as InfoType};
use noodles::vcf::io::writer::Writer as VcfWriter;
use noodles::vcf::variant::RecordBuf as VcfRecord;
use noodles::vcf::variant::io::Write as VcfWrite;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record_buf::info::field;
use noodles::vcf::{Header as VcfHeader, header::record::value::map::Map};
use serde::{Deserialize, Serialize};
use strum::Display;
use thousands::Separable;
use tokio::io::{AsyncWrite, AsyncWriteExt};

pub mod binning;
pub mod cadd;
mod compound;
pub mod provider;
pub(crate) mod reference;
pub mod spliceai;

pub(crate) mod clinvar;
pub mod consequence;
pub mod dbsnp;
pub(crate) mod frequency;

#[derive(Debug, Clone, clap::ValueEnum, PartialEq, Eq)]
pub enum OutputFormat {
    Vcf,
    Jsonl,
}

/// Command line arguments for `annotate seqvars` sub command.
#[derive(Parser, Debug)]
#[command(about = "Annotate sequence variant VCF files", long_about = None)]
pub struct Args {
    /// Assembly to use.
    #[arg(long)]
    pub assembly: Option<String>,

    /// Reference genome FASTA file (with accompanying index), strongly recommended.
    #[arg(long)]
    pub reference: Option<PathBuf>,

    /// Read the reference genome into memory.
    #[arg(long, requires = "reference")]
    pub in_memory_reference: bool,

    /// Path to the input VCF file.
    ///
    /// Use '-' to read from stdin.
    #[arg(short = 'i', long, required = true)]
    pub input: String,

    /// Path to the output file. Defaults to stdout.
    #[arg(short = 'o', long, default_value = "-")]
    pub output: String,

    #[arg(long, value_enum, default_value_t = OutputFormat::Vcf)]
    pub output_format: OutputFormat,

    /// For debug purposes, maximal number of variants to annotate.
    #[arg(long)]
    pub max_var_count: Option<usize>,

    /// What to annotate and which source to use.
    #[command(flatten)]
    pub sources: Sources,

    /// Transcript annotation related settings
    #[command(flatten)]
    pub predictor_settings: PredictorSettings,

    /// Number of threads to use for annotation.
    ///
    /// A good value is around 16.
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
}

#[derive(Debug, Display, Copy, Clone, clap::ValueEnum, PartialEq, Eq, parse_display::FromStr)]
pub enum AnnotationOption {
    Transcripts,
    Frequencies,
    ClinVar,
}

fn build_header(
    header_in: &VcfHeader,
    with_annotations: bool,
    with_frequencies: bool,
    with_clinvar: bool,
    annotators: &[AnnotatorEnum],
    additional_records: &[(String, String)],
    custom_columns: &[String],
) -> VcfHeader {
    let mut header_out = header_in.clone();
    *header_out.file_format_mut() = FileFormat::default();

    if with_frequencies {
        header_out.infos_mut().insert(
            "gnomad_exomes_an".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of alleles in gnomAD exomes",
            ),
        );
        header_out.infos_mut().insert(
            "gnomad_exomes_hom".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of hom. alt. carriers in gnomAD exomes",
            ),
        );
        header_out.infos_mut().insert(
            "gnomad_exomes_het".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of het. alt. carriers in gnomAD exomes",
            ),
        );
        header_out.infos_mut().insert(
            "gnomad_exomes_hemi".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of hemi. alt. carriers in gnomAD exomes",
            ),
        );

        header_out.infos_mut().insert(
            "gnomad_genomes_an".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of alleles in gnomAD genomes",
            ),
        );
        header_out.infos_mut().insert(
            "gnomad_genomes_hom".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of hom. alt. carriers in gnomAD genomes",
            ),
        );
        header_out.infos_mut().insert(
            "gnomad_genomes_het".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of het. alt. carriers in gnomAD genomes",
            ),
        );
        header_out.infos_mut().insert(
            "gnomad_genomes_hemi".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of hemi. alt. carriers in gnomAD genomes",
            ),
        );

        header_out.infos_mut().insert(
            "helix_an".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of alleles in HelixMtDb",
            ),
        );
        header_out.infos_mut().insert(
            "helix_hom".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of hom. alt. carriers in HelixMtDb",
            ),
        );
        header_out.infos_mut().insert(
            "helix_het".into(),
            Map::<Info>::new(
                Number::Count(1),
                InfoType::Integer,
                "Number of het. alt. carriers in HelixMtDb",
            ),
        );
    }

    if with_annotations {
        let fields = AnnField::ann_field_names(custom_columns).join(" | ");
        header_out.infos_mut().insert(
            "ANN".into(),
            Map::<Info>::new(
                Number::Unknown,
                InfoType::String,
                format!("Functional annotations: '{fields}'"),
            ),
        );
    }

    if with_clinvar {
        header_out.infos_mut().insert(
            "clinvar_germline_classification".into(),
            Map::<Info>::new(
                Number::Unknown,
                InfoType::String,
                "ClinVar clinical significance",
            ),
        );
        header_out.infos_mut().insert(
            "clinvar_vcv".into(),
            Map::<Info>::new(Number::Count(1), InfoType::String, "ClinVar VCV accession"),
        );
    }

    for annotator in annotators {
        if let Err(e) = annotator.register_vcf_headers(&mut header_out) {
            tracing::error!("Failed to register VCF headers for annotator: {}", e);
        }
    }

    for (key, value) in additional_records {
        header_out
            .insert(
                key.parse().expect("invalid key"),
                noodles::vcf::header::record::Value::from(value.as_ref()),
            )
            .unwrap();
    }

    header_out
}

/// Return path component for the assembly.
pub fn path_component(assembly: Assembly) -> &'static str {
    match assembly {
        Assembly::Grch37 | Assembly::Grch37p10 => "grch37",
        Assembly::Grch38 => "grch38",
    }
}

/// Trait for writing out annotated VCF records as VCF or VarFish TSV.
///
/// TODO: use async_trait crate
pub trait AsyncAnnotatedVariantWriter {
    #[allow(async_fn_in_trait)]
    async fn write_noodles_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error>;

    #[allow(async_fn_in_trait)]
    async fn write_annotated_record(
        &mut self,
        header: &VcfHeader,
        record: AnnotatedVariant,
    ) -> Result<(), anyhow::Error>;

    #[allow(async_fn_in_trait)]
    async fn shutdown(&mut self) -> Result<(), anyhow::Error>;

    fn set_assembly(&mut self, _assembly: String) {
        // nop
    }

    fn set_custom_columns(&mut self, _custom_columns: Vec<String>) {}

    /// Returns the conversion that `write_annotated_record` applies to a record before writing.
    ///
    /// Callers may apply it beforehand, e.g., on worker threads; writing a converted record
    /// gives the same output.
    fn record_preparer(
        &self,
    ) -> impl Fn(AnnotatedVariant) -> AnnotatedVariant + Send + Sync + 'static {
        |record| record
    }
}

pub(crate) fn prepare_vcf_record(
    record: AnnotatedVariant,
    custom_columns: &[String],
) -> noodles::vcf::variant::RecordBuf {
    let mut out_record = record.vcf;

    macro_rules! insert_field {
        ($key:literal, scalar $variant:ident, $val:expr) => {
            out_record
                .info_mut()
                .insert($key.into(), Some(field::Value::$variant($val)));
        };
        ($key:literal, opt_scalar $variant:ident, $val:expr) => {
            if let Some(v) = $val {
                out_record
                    .info_mut()
                    .insert($key.into(), Some(field::Value::$variant(v)));
            }
        };
        ($key:literal, array $array_variant:ident, $iter:expr) => {
            out_record.info_mut().insert(
                $key.into(),
                Some(field::Value::Array(field::value::Array::$array_variant(
                    $iter.map(|x| Some(x.clone())).collect(),
                ))),
            );
        };
    }

    if let Some(freqs) = &record.annotation.frequencies {
        insert_field!("gnomad_exomes_an", scalar Integer, freqs.gnomad_exomes_an);
        insert_field!("gnomad_exomes_hom", scalar Integer, freqs.gnomad_exomes_hom);
        insert_field!("gnomad_exomes_het", scalar Integer, freqs.gnomad_exomes_het);
        insert_field!("gnomad_exomes_hemi", opt_scalar Integer, freqs.gnomad_exomes_hemi);

        insert_field!("gnomad_genomes_an", scalar Integer, freqs.gnomad_genomes_an);
        insert_field!("gnomad_genomes_hom", scalar Integer, freqs.gnomad_genomes_hom);
        insert_field!("gnomad_genomes_het", scalar Integer, freqs.gnomad_genomes_het);
        insert_field!("gnomad_genomes_hemi", opt_scalar Integer, freqs.gnomad_genomes_hemi);

        insert_field!("helix_an", opt_scalar Integer, freqs.helix_an);
        insert_field!("helix_hom", opt_scalar Integer, freqs.helix_hom);
        insert_field!("helix_het", opt_scalar Integer, freqs.helix_het);
    }

    if let Some(clinvar) = &record.annotation.clinvar {
        insert_field!("clinvar_vcv", array String, clinvar.vcv.iter());
        insert_field!("clinvar_germline_classification", array String, clinvar.germline_classification.iter());
    }

    if let Some(cadd) = &record.annotation.cadd {
        insert_field!("CADD_RAW", scalar Float, cadd.raw_score);
        insert_field!("CADD_PHRED", scalar Float, cadd.phred);
    }

    if let Some(spliceai) = &record.annotation.spliceai {
        let preds = &spliceai.predictions;
        insert_field!("SpliceAI_pred_SYMBOL", array String, preds.iter().map(|p| &p.symbol));
        insert_field!("SpliceAI_pred_DS_AG", array Float, preds.iter().map(|p| &p.ds_ag));
        insert_field!("SpliceAI_pred_DS_AL", array Float, preds.iter().map(|p| &p.ds_al));
        insert_field!("SpliceAI_pred_DS_DG", array Float, preds.iter().map(|p| &p.ds_dg));
        insert_field!("SpliceAI_pred_DS_DL", array Float, preds.iter().map(|p| &p.ds_dl));
        insert_field!("SpliceAI_pred_DP_AG", array Integer, preds.iter().map(|p| &p.dp_ag));
        insert_field!("SpliceAI_pred_DP_AL", array Integer, preds.iter().map(|p| &p.dp_al));
        insert_field!("SpliceAI_pred_DP_DG", array Integer, preds.iter().map(|p| &p.dp_dg));
        insert_field!("SpliceAI_pred_DP_DL", array Integer, preds.iter().map(|p| &p.dp_dl));
    }

    if let Some(dbsnp) = &record.annotation.dbsnp {
        out_record.ids_mut().as_mut().insert(dbsnp.rs_id.clone());
    }

    if !record.annotation.consequences.is_empty() {
        let formatted_anns = record
            .annotation
            .consequences
            .iter()
            .map(|ann| ann.format(custom_columns));
        insert_field!("ANN", array String, formatted_anns);
    }

    out_record
}

/// Implement `AnnotatedVcfWriter` for `VcfWriter`.
impl<Inner: Write> AsyncAnnotatedVariantWriter for VcfWriter<Inner> {
    async fn write_noodles_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error> {
        self.write_header(header)
            .map_err(|e| anyhow::anyhow!("Error writing VCF header: {}", e))
    }

    async fn write_annotated_record(
        &mut self,
        header: &VcfHeader,
        record: AnnotatedVariant,
    ) -> Result<(), anyhow::Error> {
        self.write_variant_record(header, &record.vcf)
            .map_err(|e| anyhow::anyhow!("Error writing VCF record: {}", e))
    }

    async fn shutdown(&mut self) -> Result<(), Error> {
        Ok(<VcfWriter<Inner>>::get_mut(self).flush()?)
    }
}

/// Implement `AsyncAnnotatedVariantWriter` for `AsyncWriter`.
impl<Inner: tokio::io::AsyncWrite + Unpin> AsyncAnnotatedVariantWriter
    for noodles::vcf::r#async::io::Writer<Inner>
{
    async fn write_noodles_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error> {
        self.write_header(header)
            .await
            .map_err(|e| anyhow::anyhow!("Error writing VCF header: {}", e))
    }

    async fn write_annotated_record(
        &mut self,
        header: &VcfHeader,
        record: AnnotatedVariant,
    ) -> Result<(), anyhow::Error> {
        self.write_variant_record(header, &record.vcf)
            .await
            .map_err(|e| anyhow::anyhow!("Error writing VCF record: {}", e))
    }

    async fn shutdown(&mut self) -> Result<(), Error> {
        Ok(<noodles::vcf::r#async::io::Writer<Inner>>::get_mut(self)
            .flush()
            .await?)
    }
}

/// Implement `AsyncAnnotatedVariantWriter` for `AsyncWriter`.
impl<Inner: tokio::io::AsyncWrite + Unpin> AsyncAnnotatedVariantWriter
    for noodles::bcf::r#async::io::Writer<Inner>
{
    async fn write_noodles_header(&mut self, header: &VcfHeader) -> Result<(), Error> {
        self.write_header(header)
            .await
            .map_err(|e| anyhow::anyhow!("Error writing VCF header: {}", e))
    }

    async fn write_annotated_record(
        &mut self,
        header: &VcfHeader,
        record: AnnotatedVariant,
    ) -> Result<(), anyhow::Error> {
        noodles::bcf::r#async::io::Writer::write_variant_record(self, header, &record.vcf)
            .await
            .map_err(|e| anyhow::anyhow!("Error writing VCF record: {}", e))
    }

    async fn shutdown(&mut self) -> Result<(), Error> {
        Ok(<noodles::bcf::r#async::io::Writer<Inner>>::get_mut(self)
            .shutdown()
            .await?)
    }
}

pub struct SeqvarsVcfWriter {
    pub inner: crate::common::noodles::VariantWriter,
    pub custom_columns: Vec<String>,
}

impl SeqvarsVcfWriter {
    pub fn new(inner: crate::common::noodles::VariantWriter) -> Self {
        Self {
            inner,
            custom_columns: Vec::new(),
        }
    }

    /// Moves the annotation into the VCF record; preparing the result again changes nothing.
    fn prepare(record: AnnotatedVariant, custom_columns: &[String]) -> AnnotatedVariant {
        AnnotatedVariant {
            vcf: prepare_vcf_record(record, custom_columns),
            annotation: VariantAnnotation {
                consequences: vec![],
                frequencies: None,
                clinvar: None,
                cadd: None,
                spliceai: None,
                dbsnp: None,
            },
        }
    }
}

impl AsyncAnnotatedVariantWriter for SeqvarsVcfWriter {
    async fn write_noodles_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error> {
        self.inner.write_noodles_header(header).await
    }

    async fn write_annotated_record(
        &mut self,
        header: &VcfHeader,
        record: AnnotatedVariant,
    ) -> Result<(), anyhow::Error> {
        self.inner
            .write_annotated_record(header, Self::prepare(record, &self.custom_columns))
            .await
    }

    async fn shutdown(&mut self) -> Result<(), anyhow::Error> {
        self.inner.shutdown().await
    }

    fn set_custom_columns(&mut self, custom_columns: Vec<String>) {
        self.custom_columns = custom_columns;
    }

    fn record_preparer(
        &self,
    ) -> impl Fn(AnnotatedVariant) -> AnnotatedVariant + Send + Sync + 'static {
        let custom_columns = self.custom_columns.clone();
        move |record| Self::prepare(record, &custom_columns)
    }
}

pub struct SeqvarJsonlWriter {
    inner: Pin<Box<dyn AsyncWrite>>,
    custom_columns: Vec<String>,
}

impl SeqvarJsonlWriter {
    pub fn new(inner: Pin<Box<dyn AsyncWrite>>) -> Self {
        Self {
            inner,
            custom_columns: Vec::new(),
        }
    }
}

impl AsyncAnnotatedVariantWriter for SeqvarJsonlWriter {
    async fn write_noodles_header(&mut self, _header: &VcfHeader) -> Result<(), Error> {
        Ok(())
    }

    async fn write_annotated_record(
        &mut self,
        _header: &VcfHeader,
        record: AnnotatedVariant,
    ) -> Result<(), anyhow::Error> {
        let chrom = record.chrom().to_string();
        let pos = record.pos();
        let reference = record.vcf.reference_bases().to_string();
        let alt: Vec<String> = record
            .vcf
            .alternate_bases()
            .as_ref()
            .iter()
            .map(|a| a.to_string())
            .collect();

        // Extract INFO fields minus the heavy ANN strings
        let mut info_map = serde_json::Map::new();
        for (key, value) in record.vcf.info().as_ref() {
            if key == "ANN" {
                continue;
            }
            let json_val = match value {
                Some(field::Value::Integer(i)) => serde_json::json!(i),
                Some(field::Value::Float(f)) => serde_json::json!(f),
                Some(field::Value::Flag) => serde_json::json!(true),
                Some(field::Value::String(s)) => serde_json::json!(s),
                Some(field::Value::Array(field::value::Array::Integer(arr))) => {
                    serde_json::json!(arr)
                }
                Some(field::Value::Array(field::value::Array::Float(arr))) => {
                    serde_json::json!(arr)
                }
                Some(field::Value::Array(field::value::Array::String(arr))) => {
                    serde_json::json!(arr)
                }
                _ => serde_json::Value::Null,
            };
            info_map.insert(key.to_string(), json_val);
        }

        let json_record = serde_json::json!({
            "chromosome": chrom,
            "position": pos,
            "reference": reference,
            "alternative": alt,
            "annotations": {
                "frequencies": record.annotation.frequencies,
                "clinvar": record.annotation.clinvar,
                "per-transcript": record.annotation.consequences,
            },
            "vcfInfo": info_map,
        });

        let mut json_bytes = serde_json::to_vec(&json_record)
            .map_err(|e| anyhow::anyhow!("JSON serialization error: {}", e))?;
        json_bytes.push(b'\n');

        self.inner
            .write_all(&json_bytes)
            .await
            .map_err(|e| anyhow::anyhow!("Error writing JSONL record: {}", e))
    }

    async fn shutdown(&mut self) -> Result<(), Error> {
        Ok(self.inner.flush().await?)
    }
    fn set_custom_columns(&mut self, custom_columns: Vec<String>) {
        self.custom_columns = custom_columns;
    }
}

#[allow(clippy::large_enum_variant)]
pub(crate) enum AnnotatorEnum {
    Frequency(FrequencyAnnotator),
    Clinvar(ClinvarAnnotator),
    Consequence(ConsequenceAnnotator),
    Cadd(cadd::CaddAnnotator),
    SpliceAi(spliceai::SpliceAiAnnotator),
    Dbsnp(dbsnp::DbsnpAnnotator),
}

impl AnnotatorEnum {
    fn annotate(&self, var: &Var, annotation: &mut VariantAnnotation) -> anyhow::Result<()> {
        match self {
            AnnotatorEnum::Frequency(a) => {
                if let Some(freqs) = a.annotate(var)? {
                    if annotation.frequencies.is_some() {
                        anyhow::bail!(
                            "Multiple frequency databases returned results for variant {:?}; \
                             annotating from multiple frequency databases is not supported",
                            var
                        );
                    }
                    annotation.frequencies = Some(freqs);
                }
                Ok(())
            }
            AnnotatorEnum::Clinvar(a) => {
                if let Some(clinvar) = a.annotate(var)? {
                    if annotation.clinvar.is_some() {
                        anyhow::bail!(
                            "Multiple ClinVar databases returned results for variant {:?}; \
                             annotating from multiple ClinVar databases is not supported",
                            var
                        );
                    }
                    annotation.clinvar = Some(clinvar);
                }
                Ok(())
            }
            AnnotatorEnum::Consequence(a) => {
                annotation.consequences = a.annotate(var)?;
                Ok(())
            }
            AnnotatorEnum::Cadd(a) => {
                if let Some(cadd) = a.annotate(var)? {
                    if annotation.cadd.is_some() {
                        anyhow::bail!(
                            "Multiple CADD databases returned results for variant {:?}",
                            var
                        );
                    }
                    annotation.cadd = Some(cadd);
                }
                Ok(())
            }
            AnnotatorEnum::SpliceAi(a) => {
                if let Some(spliceai) = a.annotate(var)? {
                    if annotation.spliceai.is_some() {
                        anyhow::bail!(
                            "Multiple SpliceAI databases returned results for variant {:?}",
                            var
                        );
                    }
                    annotation.spliceai = Some(spliceai);
                }
                Ok(())
            }
            AnnotatorEnum::Dbsnp(a) => {
                if let Some(dbsnp) = a.annotate(var)? {
                    if annotation.dbsnp.is_some() {
                        anyhow::bail!(
                            "Multiple dbSNP databases returned results for variant {:?}",
                            var
                        );
                    }
                    annotation.dbsnp = Some(dbsnp);
                }
                Ok(())
            }
        }
    }

    fn register_vcf_headers(&self, header: &mut VcfHeader) -> Result<(), Error> {
        match self {
            AnnotatorEnum::Frequency(_) => Ok(()),
            AnnotatorEnum::Clinvar(_) => Ok(()),
            AnnotatorEnum::Consequence(_) => Ok(()),
            AnnotatorEnum::Cadd(_) => {
                header.infos_mut().insert(
                    "CADD_RAW".into(),
                    Map::<Info>::new(Number::Count(1), InfoType::Float, "Raw CADD score"),
                );
                header.infos_mut().insert(
                    "CADD_PHRED".into(),
                    Map::<Info>::new(Number::Count(1), InfoType::Float, "PHRED scaled CADD score"),
                );
                Ok(())
            }
            AnnotatorEnum::SpliceAi(_) => {
                let fields = vec![
                    ("SYMBOL", InfoType::String, "gene symbol"),
                    ("DS_AG", InfoType::Float, "delta score for acceptor gain"),
                    ("DS_AL", InfoType::Float, "delta score for acceptor loss"),
                    ("DS_DG", InfoType::Float, "delta score for donor gain"),
                    ("DS_DL", InfoType::Float, "delta score for donor loss"),
                    (
                        "DP_AG",
                        InfoType::Integer,
                        "delta position for acceptor gain",
                    ),
                    (
                        "DP_AL",
                        InfoType::Integer,
                        "delta position for acceptor loss",
                    ),
                    ("DP_DG", InfoType::Integer, "delta position for donor gain"),
                    ("DP_DL", InfoType::Integer, "delta position for donor loss"),
                ];
                for (name, ty, desc) in fields {
                    let field_name = format!("SpliceAI_pred_{}", name);
                    header.infos_mut().insert(
                        field_name,
                        Map::<Info>::new(
                            Number::Unknown,
                            ty,
                            format!("SpliceAI prediction. {}", desc),
                        ),
                    );
                }
                Ok(())
            }
            AnnotatorEnum::Dbsnp(_) => {
                // dbsnp just inserts IDs into the ID column
                Ok(())
            }
        }
    }
}

pub(crate) struct Annotator {
    annotators: Vec<AnnotatorEnum>,
}

impl Annotator {
    fn versions_for_vcf_header(&self) -> Vec<(String, String)> {
        // TODO also extract version information for frequencies and clinvar

        let tx_db_version = self
            .annotators
            .iter()
            .filter_map(|a| match a {
                AnnotatorEnum::Consequence(a) => a.predictor.data_version(),
                _ => None,
            })
            .next();
        tx_db_version
            .map(|v| vec![("mehariTxDbVersion".to_string(), v)])
            .unwrap_or_default()
    }
}

impl Annotator {
    pub(crate) fn new(annotators: Vec<AnnotatorEnum>) -> Self {
        Self { annotators }
    }

    pub fn annotate(&self, vcf_record: &VcfRecord) -> anyhow::Result<VariantAnnotation> {
        let vcf_var = from_vcf_allele(vcf_record, 0);
        let mut annotation = VariantAnnotation {
            consequences: vec![],
            frequencies: None,
            clinvar: None,
            cadd: None,
            spliceai: None,
            dbsnp: None,
        };

        // Skip records with a deletion as alternative allele.
        if vcf_var.alternative == "*" {
            return Ok(annotation);
        }

        for annotator in &self.annotators {
            annotator.annotate(&vcf_var, &mut annotation)?;
        }

        Ok(annotation)
    }

    pub(crate) fn consequence_predictor(&self) -> Option<&ConsequencePredictor> {
        self.annotators.iter().find_map(|a| match a {
            AnnotatorEnum::Consequence(c) => Some(&c.predictor),
            _ => None,
        })
    }
}

/// Main entry point for `annotate seqvars` sub command.
pub async fn run(_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("config = {:#?}", &args);

    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global();

    match args.output_format {
        OutputFormat::Vcf => {
            let writer = open_variant_writer(&args.output).await?;
            let mut seqvars_writer = SeqvarsVcfWriter::new(writer);
            run_with_writer(&mut seqvars_writer, args).await?;
            // shutdown is already called inside run_with_writer
        }
        OutputFormat::Jsonl => {
            let stream = crate::common::io::tokio::open_write_maybe_bgzf(&args.output).await?;
            let mut writer = SeqvarJsonlWriter::new(stream);
            run_with_writer(&mut writer, args).await?;
            // shutdown is already called inside run_with_writer
        }
    }

    Ok(())
}

/// Run the annotation with the given `Write` within the `VcfWriter`.
async fn run_with_writer(
    writer: &mut impl AsyncAnnotatedVariantWriter,
    args: &Args,
) -> Result<(), Error> {
    tracing::info!("Open VCF and read header");
    let mut reader = open_variant_reader(&args.input).await?;

    let mut header_in = reader.read_header().await?;

    // Work around glnexus issue with RNC.
    if let Some(format) = header_in.formats_mut().get_mut("RNC") {
        *format.number_mut() = FormatNumber::Count(1);
        *format.type_mut() = FormatType::String;
    }

    let tx_dbs = match &args.sources.transcripts {
        Some(sources) => {
            tracing::info!("Opening transcript database(s)");
            sources
                .iter()
                .map(consequence::load_tx_db)
                .collect::<Result<Vec<_>, _>>()?
        }
        None => vec![],
    };

    for db in &tx_dbs {
        for sv in &db.source_version {
            tracing::info!("Transcript source version: {:?}", sv);
        }
    }

    let extract_unstructured_first = |key: &str| -> Option<String> {
        header_in
            .other_records()
            .get(key)
            .and_then(|collection| match collection {
                noodles::vcf::header::record::value::Collection::Unstructured(list) => {
                    list.first().cloned()
                }
                noodles::vcf::header::record::value::Collection::Structured(_) => {
                    // If it's surprisingly structured, we can't extract a simple string hint
                    None
                }
            })
    };

    // Priority 1: Check the ##contig lines for an `assembly=` field
    let vcf_contig_assembly = header_in
        .contigs()
        .values()
        .find_map(|contig| contig.other_fields().get("assembly").map(|v| v.to_string()));

    // Priority 2: Extract the file stem from ##reference=...
    let vcf_reference_tag = extract_unstructured_first("reference").map(|v| {
        let path = Path::new(&v);
        path.file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string()
    });

    // Priority 3: Extract the ##assembly tag, handling URLs and plain strings
    let vcf_assembly_tag = extract_unstructured_first("assembly").map(|v| {
        if v.starts_with("http") || v.starts_with("ftp") || v.starts_with("file") {
            let path = Path::new(&v);
            path.file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string()
        } else {
            v
        }
    });

    // Resolve the best available hint
    let vcf_hint = vcf_contig_assembly
        .or(vcf_reference_tag)
        .or(vcf_assembly_tag)
        .unwrap_or_default();

    // 3. Resolve the Assembly (Using the fuzzy matching logic from before)
    let assembly = match &args.assembly {
        Some(asm) => asm.clone(),
        None => {
            let mut unique_assemblies: IndexMap<String, String> = IndexMap::new();
            for db in &tx_dbs {
                for sv in &db.source_version {
                    let mut final_assembly = sv.assembly.clone();
                    #[allow(deprecated)]
                    if final_assembly.is_empty() {
                        final_assembly = _assembly_from_assembly_enum(sv.assembly_enum);
                    }
                    if !final_assembly.is_empty() {
                        unique_assemblies.insert(final_assembly.to_lowercase(), final_assembly);
                    }
                }
            }

            if unique_assemblies.len() == 1 {
                // Only one DB assembly available, just use it
                unique_assemblies.into_values().next().unwrap()
            } else if !unique_assemblies.is_empty() && !vcf_hint.is_empty() {
                // Multiple DBs available. See if the VCF hint matches one of them!
                let hint_lower = vcf_hint.to_lowercase();

                // Do a fuzzy check: Does the VCF hint contain the DB assembly name, or vice versa?
                let matched_assembly =
                    unique_assemblies.iter().find_map(|(lower_asm, orig_asm)| {
                        if hint_lower.contains(lower_asm) || lower_asm.contains(&hint_lower) {
                            Some(orig_asm.clone())
                        } else {
                            None
                        }
                    });

                if let Some(matched) = matched_assembly {
                    tracing::info!(
                        "Auto-detected assembly {:?} from VCF header metadata",
                        matched
                    );
                    matched
                } else {
                    let list: Vec<_> = unique_assemblies.into_values().collect();
                    anyhow::bail!(
                        "Multiple assemblies found in transcript databases ({}), but the VCF header hint ('{}') didn't clearly match any. Please specify --assembly manually.",
                        list.join(", "),
                        vcf_hint
                    );
                }
            } else if unique_assemblies.is_empty() {
                anyhow::bail!(
                    "No valid assemblies found in transcript databases and --assembly omitted."
                );
            } else {
                let list: Vec<_> = unique_assemblies.into_values().collect();
                anyhow::bail!(
                    "Multiple assemblies found in databases ({}), and no hints found in VCF header. Please specify --assembly.",
                    list.join(", ")
                );
            }
        }
    };

    writer.set_assembly(assembly.clone());
    tracing::info!("Using assembly {:?}", &assembly);

    let annotator = Arc::new(setup_seqvars_annotator(
        &args.sources,
        tx_dbs,
        args.reference.as_ref(),
        args.in_memory_reference,
        &args.predictor_settings,
        assembly,
    )?);

    // With a reference genome, the REF of each record is compared with it.
    let reference_provider = annotator
        .consequence_predictor()
        .map(|predictor| predictor.provider.clone())
        .filter(|provider| provider.reference_available());

    let mut additional_header_info = annotator.versions_for_vcf_header();
    additional_header_info.push(("mehariCmd".into(), env::args().join(" ")));
    additional_header_info.push(("mehariVersion".into(), env!("CARGO_PKG_VERSION").into()));

    let with_annotations = args
        .sources
        .transcripts
        .as_ref()
        .is_some_and(|v| !v.is_empty());
    let with_frequencies = args
        .sources
        .frequencies
        .as_ref()
        .is_some_and(|v| !v.is_empty());
    let with_clinvar = args.sources.clinvar.as_ref().is_some_and(|v| !v.is_empty());

    let enable_compound = args
        .predictor_settings
        .compound_settings
        .enable_compound_variants;

    let custom_columns = annotator
        .consequence_predictor()
        .map(|predictor| predictor.custom_columns().to_vec())
        .unwrap_or_default();

    writer.set_custom_columns(custom_columns.clone());

    let header_out = build_header(
        &header_in,
        with_annotations,
        with_frequencies,
        with_clinvar,
        &annotator.annotators,
        &additional_header_info,
        &custom_columns,
    );

    // Perform the VCF annotation. The compound processor moves into each batch's blocking task
    // and back, so that compound predictions do not run on the async runtime.
    let mut processor = enable_compound.then(|| {
        VariantProcessor::new(
            annotator.clone(),
            &args.predictor_settings.compound_settings,
        )
    });

    tracing::info!("Annotating VCF ...");
    let start = Instant::now();
    let mut prev = Instant::now();
    let mut total_read = 0usize;
    let mut total_written = 0usize;

    writer.write_noodles_header(&header_out).await?;

    use futures::TryStreamExt;
    let mut records = reader.records(&header_in).await;

    let worker_threads = rayon::current_num_threads();
    let batch_size = worker_threads * 1024;

    // Preparing records for output is costly, so do it on the worker threads, unless there is
    // only one: then the writer has time to spare.  Compound variant detection still needs the
    // annotations, so the writer prepares those records itself.
    let preparer =
        (!enable_compound && worker_threads > 1).then(|| Arc::new(writer.record_preparer()));

    let mut next_batch = Vec::with_capacity(batch_size);
    let mut processing_handle: Option<
        tokio::task::JoinHandle<anyhow::Result<Vec<AnnotatedVariant>>>,
    > = None;

    loop {
        // fill the `next_batch` buffer asynchronously (happens concurrently while the previous batch is being processed)
        while next_batch.len() < batch_size {
            if let Some(max) = args.max_var_count
                && total_read >= max
            {
                break;
            }

            match records.try_next().await? {
                Some(vcf_record) => {
                    if vcf_record.alternate_bases().len() != 1 {
                        tracing::error!(
                            "Found record with more than one alternate allele.  This is currently not supported. \
                    Please use `bcftools norm` to split multi-allelic records.  Record: {:?}",
                            &vcf_record
                        );
                        anyhow::bail!("multi-allelic records not supported");
                    }
                    if let Some(provider) = &reference_provider
                        && ref_differs(provider, &vcf_record) == Some(true)
                    {
                        tracing::warn!(
                            "REF at {}:{} differs from the reference genome; the annotation uses \
                             the reference",
                            vcf_record.reference_sequence_name(),
                            vcf_record.variant_start().map_or(0, |pos| pos.get())
                        );
                    }
                    next_batch.push(vcf_record);
                    total_read += 1;
                }
                _ => break, // EOF
            }
        }

        // start processing `next_batch` before writing the previous batch, so that writing
        // overlaps with processing
        let previous_handle = processing_handle.take();
        if !next_batch.is_empty() {
            let batch_to_process =
                std::mem::replace(&mut next_batch, Vec::with_capacity(batch_size));
            let annotator_clone = annotator.clone();
            let preparer = preparer.clone();

            processing_handle = Some(tokio::task::spawn_blocking(move || {
                use rayon::prelude::*;

                batch_to_process
                    .into_par_iter()
                    .map(|record| {
                        let annotation = annotator_clone.annotate(&record)?;
                        let annotated = AnnotatedVariant {
                            vcf: record,
                            annotation,
                        };
                        Ok(match &preparer {
                            Some(prepare) => prepare(annotated),
                            None => annotated,
                        })
                    })
                    .collect::<anyhow::Result<Vec<_>>>()
            }));
        }

        // await the previous batch's results and write them out sequentially
        if let Some(handle) = previous_handle {
            let processed_batch = handle.await??;
            // Compound predictions run on the blocking pool too, while the next batch is
            // being annotated.
            let ready_records = match processor.take() {
                Some(mut batch_processor) => {
                    let (ready_records, batch_processor) =
                        tokio::task::spawn_blocking(move || -> anyhow::Result<_> {
                            let mut ready_records = Vec::new();
                            for annotated_record in processed_batch {
                                ready_records.extend(batch_processor.process_annotated_record(
                                    annotated_record.vcf,
                                    annotated_record.annotation,
                                )?);
                            }
                            Ok((ready_records, batch_processor))
                        })
                        .await??;
                    processor = Some(batch_processor);
                    ready_records
                }
                None => processed_batch,
            };

            for annotated_record in ready_records {
                if prev.elapsed().as_secs() >= 60 {
                    tracing::info!("at {:?}", from_vcf_allele(&annotated_record.vcf, 0));
                    prev = Instant::now();
                }
                writer
                    .write_annotated_record(&header_out, annotated_record)
                    .await?;
                total_written += 1;
            }
        }

        // if there are no more records to process, stop.
        if processing_handle.is_none() {
            break;
        }
    }

    if let Some(mut processor) = processor {
        let final_records = tokio::task::spawn_blocking(move || processor.flush()).await??;
        for rec in final_records {
            writer.write_annotated_record(&header_out, rec).await?;
            total_written += 1;
        }
    }

    tracing::info!(
        "... annotated {} records in {:?}",
        total_written.separate_with_commas(),
        start.elapsed()
    );
    if let Some(predictor) = annotator.consequence_predictor() {
        predictor.log_restrict_summary();
    }
    writer.shutdown().await?;
    Ok(())
}

/// Whether the REF of `record` differs from the reference genome; `None` if the reference genome
/// lacks the bases.
fn ref_differs(provider: &provider::Provider, record: &VcfRecord) -> Option<bool> {
    let start = record.variant_start()?.get() as u64 - 1;
    let reference = record.reference_bases();
    let end = start + reference.len() as u64;
    let bases = provider.reference_bases(record.reference_sequence_name(), start, end)?;
    Some(!bases.eq_ignore_ascii_case(reference.as_bytes()))
}

pub(crate) fn setup_seqvars_annotator(
    sources: &Sources,
    preloaded_tx_dbs: Vec<TxSeqDatabase>,
    reference: Option<impl AsRef<Path>>,
    in_memory_reference: bool,
    predictor_settings: &PredictorSettings,
    assembly: String,
) -> Result<Annotator, Error> {
    let contig_manager = Arc::new(ContigManager::new(&assembly));
    let mut annotators = vec![];

    // Add the frequency annotator if requested.
    if let Some(rocksdb_paths) = &sources.frequencies {
        let freq_dbs = initialize_frequency_annotators_for_assembly(
            rocksdb_paths,
            &assembly,
            contig_manager.clone(),
        )?;
        if freq_dbs.is_empty() && !rocksdb_paths.is_empty() {
            anyhow::bail!(
                "Frequency paths were provided but no frequency databases could be initialized for assembly '{}'",
                assembly
            );
        }
        for freq_db in freq_dbs {
            annotators.push(AnnotatorEnum::Frequency(freq_db))
        }
    }

    // Add the ClinVar annotator if requested.
    if let Some(rocksdb_paths) = &sources.clinvar {
        let clinvar_dbs = initialize_clinvar_annotators_for_assembly(
            rocksdb_paths,
            &assembly,
            contig_manager.clone(),
        )?;
        if clinvar_dbs.is_empty() && !rocksdb_paths.is_empty() {
            anyhow::bail!(
                "ClinVar paths were provided but no ClinVar databases could be initialized for assembly '{}'",
                assembly
            );
        }
        for clinvar_db in clinvar_dbs {
            annotators.push(AnnotatorEnum::Clinvar(clinvar_db))
        }
    }

    // Add the CADD annotator if requested.
    if let Some(rocksdb_paths) = &sources.cadd {
        let cadd_dbs = initialize_cadd_annotators_for_assembly(
            rocksdb_paths,
            &assembly,
            contig_manager.clone(),
        )?;
        for cadd_db in cadd_dbs {
            annotators.push(AnnotatorEnum::Cadd(cadd_db))
        }
    }

    // Add the SpliceAI annotator if requested.
    if let Some(rocksdb_paths) = &sources.spliceai {
        let spliceai_dbs = initialize_spliceai_annotators_for_assembly(
            rocksdb_paths,
            &assembly,
            contig_manager.clone(),
        )?;
        for spliceai_db in spliceai_dbs {
            annotators.push(AnnotatorEnum::SpliceAi(spliceai_db))
        }
    }

    // Add the dbSNP annotator if requested.
    if let Some(rocksdb_paths) = &sources.dbsnp {
        let dbsnp_dbs = initialize_dbsnp_annotators_for_assembly(
            rocksdb_paths,
            &assembly,
            contig_manager.clone(),
        )?;
        for dbsnp_db in dbsnp_dbs {
            annotators.push(AnnotatorEnum::Dbsnp(dbsnp_db))
        }
    }

    // Add the consequence annotator if requested.
    if !preloaded_tx_dbs.is_empty() {
        // Filter out any loaded databases that don't match the active assembly
        let databases: Vec<_> = preloaded_tx_dbs
            .into_iter()
            .filter(|db| {
                db.source_version
                    .iter()
                    .any(|sv| sv.assembly.eq_ignore_ascii_case(&assembly))
            })
            .collect();

        if databases.is_empty() {
            tracing::warn!(
                "No suitable transcript databases found for requested assembly {:?}, therefore no consequence prediction will occur.",
                &assembly
            );
        } else {
            let tx_db = merge_transcript_databases(databases)?;
            if let Some(tx_sources) = &sources.transcripts {
                tracing::info!(
                    "Loaded transcript database(s) from {}",
                    &tx_sources.join(", ")
                );
            }
            annotators.push(
                ConsequenceAnnotator::from_db_and_settings(
                    tx_db,
                    reference,
                    in_memory_reference,
                    predictor_settings,
                )
                .map(AnnotatorEnum::Consequence)?,
            );
        }
    } else {
        if let Some(tx_sources) = &sources.transcripts {
            tracing::info!("Opening transcript database(s)");

            let databases = consequence::load_transcript_dbs_for_assembly(tx_sources, &assembly)?;

            if databases.is_empty() {
                tracing::warn!(
                    "No suitable transcript databases found for requested assembly {:?}, therefore no consequence prediction will occur.",
                    &assembly
                );
            } else {
                let tx_db = merge_transcript_databases(databases)?;
                tracing::info!(
                    "Loaded transcript database(s) from {}",
                    &tx_sources.join(", ")
                );
                annotators.push(
                    ConsequenceAnnotator::from_db_and_settings(
                        tx_db,
                        reference,
                        in_memory_reference,
                        predictor_settings,
                    )
                    .map(AnnotatorEnum::Consequence)?,
                );
            }
        }
    }

    let annotator = Annotator::new(annotators);
    Ok(annotator)
}

/// Dynamically reads the database's internal metadata to verify its assembly.
/// Falls back to a loose path-string validation warning if metadata is unreadable.
fn verify_database_assembly(path: &str, expected_assembly: &str) {
    let options = rocksdb::Options::default();

    if let Ok(cfs) = rocksdb::DB::list_cf(&options, path)
        && let Ok(db) = rocksdb::DB::open_cf_for_read_only(&options, path, &cfs, false)
        && let Some(cf_meta) = db.cf_handle("meta")
        && let Ok(Some(bytes)) = db.get_cf(&cf_meta, b"assembly")
        && let Ok(actual_assembly) = std::str::from_utf8(&bytes)
    {
        if !actual_assembly.eq_ignore_ascii_case(expected_assembly) {
            tracing::warn!(
                "Database assembly mismatch! The database at '{}' internal metadata specifies assembly '{}', but the active annotation pipeline is configured for '{}'.",
                path,
                actual_assembly,
                expected_assembly
            );
        }
        return;
    }

    let path_lower = path.to_lowercase();
    let assembly_lower = expected_assembly.to_lowercase();
    if !path_lower.contains(&assembly_lower) {
        tracing::warn!(
            "Database path '{}' does not contain the active assembly name string '{}'. Please ensure it matches your target genome build configuration.",
            path,
            expected_assembly
        );
    }
}
pub fn initialize_clinvar_annotators_for_assembly(
    rocksdb_paths: &[String],
    assembly: &str,
    contig_manager: Arc<ContigManager>,
) -> Result<Vec<ClinvarAnnotator>, Error> {
    let mut annotators = Vec::new();

    for rocksdb_path in rocksdb_paths {
        verify_database_assembly(rocksdb_path, assembly);
        tracing::info!(
            "Loading ClinVar database for assembly {:?} from {}",
            assembly,
            rocksdb_path
        );
        let annotator = ClinvarAnnotator::from_path(rocksdb_path, contig_manager.clone())?;
        annotators.push(annotator);
    }
    Ok(annotators)
}

pub fn initialize_frequency_annotators_for_assembly(
    rocksdb_paths: &[String],
    assembly: &str,
    contig_manager: Arc<ContigManager>,
) -> Result<Vec<FrequencyAnnotator>, Error> {
    let mut annotators = Vec::new();

    for rocksdb_path in rocksdb_paths {
        verify_database_assembly(rocksdb_path, assembly);
        tracing::info!(
            "Loading frequency database for assembly {:?} from {}",
            assembly,
            rocksdb_path
        );
        let annotator = FrequencyAnnotator::from_path(rocksdb_path, contig_manager.clone())?;
        annotators.push(annotator);
    }
    Ok(annotators)
}

pub fn initialize_cadd_annotators_for_assembly(
    rocksdb_paths: &[String],
    assembly: &str,
    contig_manager: Arc<ContigManager>,
) -> Result<Vec<cadd::CaddAnnotator>, Error> {
    let mut annotators = Vec::new();

    for rocksdb_path in rocksdb_paths {
        verify_database_assembly(rocksdb_path, assembly);
        tracing::info!(
            "Loading CADD database for assembly {:?} from {}",
            assembly,
            rocksdb_path
        );
        let annotator = cadd::CaddAnnotator::from_path(rocksdb_path, contig_manager.clone())?;
        annotators.push(annotator);
    }
    Ok(annotators)
}

pub fn initialize_spliceai_annotators_for_assembly(
    rocksdb_paths: &[String],
    assembly: &str,
    contig_manager: Arc<ContigManager>,
) -> Result<Vec<spliceai::SpliceAiAnnotator>, Error> {
    let mut annotators = Vec::new();

    for rocksdb_path in rocksdb_paths {
        verify_database_assembly(rocksdb_path, assembly);
        tracing::info!(
            "Loading SpliceAI database for assembly {:?} from {}",
            assembly,
            rocksdb_path
        );
        let annotator =
            spliceai::SpliceAiAnnotator::from_path(rocksdb_path, contig_manager.clone())?;
        annotators.push(annotator);
    }
    Ok(annotators)
}

pub fn initialize_dbsnp_annotators_for_assembly(
    rocksdb_paths: &[String],
    assembly: &str,
    contig_manager: Arc<ContigManager>,
) -> Result<Vec<dbsnp::DbsnpAnnotator>, Error> {
    let mut annotators = Vec::new();

    for rocksdb_path in rocksdb_paths {
        verify_database_assembly(rocksdb_path, assembly);
        tracing::info!(
            "Loading dbSNP database for assembly {:?} from {}",
            assembly,
            rocksdb_path
        );
        let annotator = dbsnp::DbsnpAnnotator::from_path(rocksdb_path, contig_manager.clone())?;
        annotators.push(annotator);
    }
    Ok(annotators)
}

fn _assembly_from_assembly_enum(val: i32) -> String {
    #[allow(deprecated)]
    match pbs::txs::Assembly::from_i32(val).unwrap_or(pbs::txs::Assembly::Unknown) {
        pbs::txs::Assembly::Unknown => "".into(),
        pbs::txs::Assembly::Grch37 => "GRCh37".into(),
        pbs::txs::Assembly::Grch38 => "GRCh38".into(),
    }
}

fn _genome_build_from_enum(val: i32) -> String {
    #[allow(deprecated)]
    match pbs::txs::GenomeBuild::from_i32(val).unwrap_or(pbs::txs::GenomeBuild::Unknown) {
        pbs::txs::GenomeBuild::Unknown => "".into(),
        pbs::txs::GenomeBuild::Grch37 => "GRCh37".into(),
        pbs::txs::GenomeBuild::Grch38 => "GRCh38".into(),
    }
}

fn _source_from_source_enum(val: i32) -> String {
    #[allow(deprecated)]
    match pbs::txs::Source::from_i32(val).unwrap_or(pbs::txs::Source::Unknown) {
        pbs::txs::Source::Unknown => "".into(),
        pbs::txs::Source::Refseq => "RefSeq".into(),
        pbs::txs::Source::Ensembl => "Ensembl".into(),
    }
}

/// Create for all alternate alleles from the given VCF record.
pub fn from_vcf_allele(value: &noodles::vcf::variant::RecordBuf, allele_no: usize) -> Var {
    let chrom = value.reference_sequence_name().to_string();
    let pos: usize = value
        .variant_start()
        .expect("Telomeric breakends not supported")
        .get();
    let pos = i32::try_from(pos).unwrap();
    let reference = value.reference_bases().to_string();
    Var {
        chrom,
        pos,
        reference,
        alternative: value.alternate_bases().as_ref()[allele_no].to_string(),
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantAnnotation {
    pub consequences: Vec<AnnField>,
    pub frequencies: Option<FreqResult>,
    pub clinvar: Option<ClinvarResult>,
    pub cadd: Option<cadd::CaddResult>,
    pub spliceai: Option<spliceai::SpliceAiResult>,
    pub dbsnp: Option<dbsnp::DbsnpResult>,
}

#[derive(Debug, Clone)]
pub struct AnnotatedVariant {
    pub vcf: VcfRecord,
    pub annotation: VariantAnnotation,
}

impl AnnotatedVariant {
    /// Convenience method to get the chromosome
    pub fn chrom(&self) -> &str {
        self.vcf.reference_sequence_name()
    }

    /// Convenience method to get the 1-based start position
    pub fn pos(&self) -> usize {
        self.vcf.variant_start().map(usize::from).unwrap_or(0)
    }
}

struct VariantProcessor {
    annotator: Arc<Annotator>,
    buffer: crate::annotate::seqvars::compound::VariantBuffer,
    annotations: std::collections::HashMap<String, VariantAnnotation>,
}

impl VariantProcessor {
    pub fn new(
        annotator: Arc<Annotator>,
        compound_settings: &crate::annotate::cli::CompoundSettings,
    ) -> Self {
        Self {
            annotator,
            buffer: compound::VariantBuffer::new(
                compound_settings.phasing_strategy,
                compound_settings.compound_window_size,
            ),
            annotations: Default::default(),
        }
    }

    /// Takes an already-annotated record, buffers it, and returns any records ready to be written.
    pub fn process_annotated_record(
        &mut self,
        record: VcfRecord,
        annotation: VariantAnnotation,
    ) -> anyhow::Result<Vec<AnnotatedVariant>> {
        let mut ready_records = Vec::new();
        let vcf_var = from_vcf_allele(&record, 0);

        if self.buffer.should_flush(&vcf_var.chrom, vcf_var.pos) {
            ready_records.extend(self.process_buffer_flush()?);
        }
        let csq_var = VcfVariant {
            chromosome: vcf_var.chrom,
            position: vcf_var.pos,
            reference: vcf_var.reference,
            alternative: vcf_var.alternative,
        };
        let key = csq_var.to_string();
        if self.annotations.insert(key.clone(), annotation).is_some() {
            tracing::warn!(
                "Duplicate SPDI detected in annotation buffer: {}; \
                 the earlier annotation was overwritten.",
                key
            );
        }

        let (tx_accessions, min_start, max_end) =
            if let Some(predictor) = self.annotator.consequence_predictor() {
                consequence::get_transcript_boundaries(
                    predictor,
                    &csq_var.chromosome,
                    csq_var.position,
                    csq_var.reference.len(),
                )
            } else {
                (std::collections::HashSet::new(), i32::MAX, -1)
            };

        // Note: We don't support multi allelic records at the moment, so this is a constant 1.
        let alt_allele_idx = 1;
        // Currently, phasing is evaluated based on the first sample in the VCF only.
        let sample_idx = 0;

        self.buffer.push(
            csq_var,
            record,
            tx_accessions,
            min_start,
            max_end,
            sample_idx,
            alt_allele_idx,
        );

        Ok(ready_records)
    }

    pub fn flush(&mut self) -> anyhow::Result<Vec<AnnotatedVariant>> {
        self.process_buffer_flush()
    }

    fn process_buffer_flush(&mut self) -> anyhow::Result<Vec<AnnotatedVariant>> {
        let (variants, compound_groups) = self.buffer.flush();
        let predictor = self.annotator.consequence_predictor();

        let mut result_annotations: Vec<VariantAnnotation> = variants
            .iter()
            .map(|b| {
                self.annotations
                    .remove(&b.vcf_var.to_string())
                    .unwrap_or_else(|| VariantAnnotation {
                        consequences: vec![],
                        frequencies: None,
                        clinvar: None,
                        cadd: None,
                        spliceai: None,
                        dbsnp: None,
                    })
            })
            .collect();

        for group_indices in compound_groups {
            if let Some(predictor) = predictor {
                let vcf_vars: Vec<_> = group_indices
                    .iter()
                    .map(|&i| variants[i].vcf_var.clone())
                    .collect();

                // The group's variants, for messages.
                let group = || {
                    vcf_vars
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                };
                match predictor.predict_multiple(&vcf_vars) {
                    Ok(Some(ann_fields)) => {
                        for &i in &group_indices {
                            result_annotations[i]
                                .consequences
                                .extend(ann_fields.clone());
                        }
                    }
                    Ok(None) => {}
                    // An invalid group cannot be combined, but its variants keep their
                    // single-variant annotations.
                    Err(SeqvarsError::GroupValidation(e)) => {
                        tracing::warn!("Skipping compound variant {}: {}", group(), e);
                    }
                    // Other errors, e.g. of a plugin, stop the run, as for single variants.
                    Err(e) => {
                        return Err(Error::from(e)
                            .context(format!("Failed to annotate compound variant {}", group())));
                    }
                }
            }
        }

        Ok(variants
            .into_iter()
            .zip(result_annotations)
            .map(|(b, ann)| AnnotatedVariant {
                vcf: b.record,
                annotation: ann,
            })
            .collect())
    }
}

#[cfg(test)]
mod test {
    use super::binning::bin_from_range;
    use super::consequence::{ConsequenceAnnotator, VcfVariant, load_tx_db};
    use super::{
        AnnotatedVariant, Annotator, AnnotatorEnum, Args, OutputFormat, VariantProcessor,
        VcfRecord, run,
    };
    use crate::annotate::cli::{CompoundSettings, PhasingStrategy};
    use crate::annotate::cli::{PickUnit, PredictorSettings};
    use crate::annotate::cli::{Sources, TranscriptSettings};
    use crate::annotate::seqvars::consequence::plugin::{BuiltinPlugin, PluginSpec};
    use crate::common::noodles::{NoodlesVariantReader, open_variant_reader};
    use clap_verbosity_flag::Verbosity;
    use futures::TryStreamExt;
    use pretty_assertions::assert_eq;
    use std::collections::HashSet;
    use std::path::{Path, PathBuf};
    use std::sync::Arc;
    use temp_testdir::TempDir;

    use noodles::core::Position;
    use noodles::vcf::variant::record_buf::AlternateBases;
    use noodles::vcf::variant::record_buf::info::field::Value;
    use noodles::vcf::variant::record_buf::info::field::value::Array;

    async fn read_vcf(
        path: impl AsRef<Path>,
    ) -> Result<Vec<noodles::vcf::variant::RecordBuf>, anyhow::Error> {
        let mut output_reader = open_variant_reader(path.as_ref()).await?;
        let header = output_reader.read_header().await?;
        let mut record_iter = output_reader.records(&header).await;
        let mut records = Vec::new();
        while let Some(record) = record_iter.try_next().await? {
            records.push(record);
        }
        Ok(records)
    }

    #[tokio::test]
    async fn smoke_test_output_vcf() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.vcf");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db";
        let assembly = "grch37";
        let args = Args {
            threads: 1,
            reference: None,
            in_memory_reference: true,
            assembly: Some(assembly.into()),
            predictor_settings: PredictorSettings {
                transcript_settings: TranscriptSettings::most_severe(PickUnit::Gene),
                ..Default::default()
            },
            input: String::from("tests/data/annotate/seqvars/brca1.examples.vcf"),
            output: path_out.into_os_string().into_string().unwrap(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
                ..Default::default()
            },
        };

        run(&args_common, &args).await?;

        let actual = std::fs::read_to_string(args.output)?;
        // remove vcf header lines starting with ##mehari
        let actual = actual
            .lines()
            .filter(|line| !line.starts_with("##mehari"))
            .collect::<Vec<_>>()
            .join("\n");
        insta::assert_snapshot!(actual);

        Ok(())
    }

    /// The ANN header lists the custom columns of all enabled features, in order.
    #[tokio::test]
    async fn custom_columns_in_ann_header() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.vcf");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let mut predictor_settings = PredictorSettings::default();
        predictor_settings.reporting_settings.plugins = [
            BuiltinPlugin::TxSequenceRef,
            BuiltinPlugin::TxSequenceAlt,
            BuiltinPlugin::AaSequenceRef,
            BuiltinPlugin::AaSequenceAlt,
            BuiltinPlugin::Nmd,
        ]
        .into_iter()
        .map(PluginSpec::Builtin)
        .collect();
        predictor_settings
            .compound_settings
            .enable_compound_variants = true;
        let args = Args {
            threads: 1,
            reference: None,
            in_memory_reference: true,
            assembly: Some("grch37".into()),
            predictor_settings,
            input: String::from("tests/data/annotate/seqvars/brca1.examples.vcf"),
            output: path_out.into_os_string().into_string().unwrap(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            sources: Sources {
                transcripts: Some(vec!["tests/data/annotate/db/grch37/txs.bin.zst".into()]),
                ..Default::default()
            },
        };

        run(&args_common, &args).await?;

        let output = std::fs::read_to_string(&args.output)?;
        let ann_header = output
            .lines()
            .find(|line| line.starts_with("##INFO=<ID=ANN,"))
            .unwrap();
        assert!(
            ann_header.contains(
                "| tx_sequence_ref | tx_sequence_alt | aa_sequence_ref | aa_sequence_alt \
                 | nmd | compound_ids | compound_variants'"
            ),
            "{ann_header}"
        );

        Ok(())
    }

    /// A plugin that fails for a compound group stops the run.
    #[cfg(feature = "wasm")]
    #[tokio::test]
    async fn compound_plugin_error_stops_the_run() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_in = temp.join("input.vcf");
        let path_out = temp.join("output.vcf");

        let header = std::fs::read_to_string("tests/data/annotate/seqvars/brca1.examples.vcf")?
            .lines()
            .filter(|line| line.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n");
        std::fs::write(
            &path_in,
            format!(
                "{header}\n\
                 17\t41197701\t.\tG\tC\t.\t.\t.\tGT\t1/1\n\
                 17\t41215920\t.\tG\tT\t.\t.\t.\tGT\t1/1\n"
            ),
        )?;

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let mut predictor_settings = PredictorSettings::default();
        predictor_settings.reporting_settings.plugins = vec![PluginSpec::Wasm {
            path: "tests/data/plugins/probe.wasm".into(),
            options: vec![("fail-compound".into(), "group went wrong".into())],
        }];
        predictor_settings
            .compound_settings
            .enable_compound_variants = true;
        let args = Args {
            threads: 1,
            reference: None,
            in_memory_reference: true,
            assembly: Some("grch37".into()),
            predictor_settings,
            input: path_in.into_os_string().into_string().unwrap(),
            output: path_out.into_os_string().into_string().unwrap(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            sources: Sources {
                transcripts: Some(vec!["tests/data/annotate/db/grch37/txs.bin.zst".into()]),
                ..Default::default()
            },
        };

        let error = run(&args_common, &args).await.unwrap_err();
        assert!(
            format!("{error:#}").contains("group went wrong"),
            "{error:#}"
        );

        Ok(())
    }

    /// Single variants are labelled `Single`, compound groups `comp_N` with their variants.
    #[tokio::test]
    async fn compound_columns() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_in = temp.join("input.vcf");
        let path_out = temp.join("output.vcf");

        // Two homozygous variants on the same transcripts form one compound group.
        let header = std::fs::read_to_string("tests/data/annotate/seqvars/brca1.examples.vcf")?
            .lines()
            .filter(|line| line.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n");
        std::fs::write(
            &path_in,
            format!(
                "{header}\n\
                 17\t41197701\t.\tG\tC\t.\t.\t.\tGT\t1/1\n\
                 17\t41215920\t.\tG\tT\t.\t.\t.\tGT\t1/1\n"
            ),
        )?;

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let mut predictor_settings = PredictorSettings::default();
        predictor_settings
            .compound_settings
            .enable_compound_variants = true;
        let args = Args {
            threads: 1,
            reference: None,
            in_memory_reference: true,
            assembly: Some("grch37".into()),
            predictor_settings,
            input: path_in.into_os_string().into_string().unwrap(),
            output: path_out.into_os_string().into_string().unwrap(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            sources: Sources {
                transcripts: Some(vec!["tests/data/annotate/db/grch37/txs.bin.zst".into()]),
                ..Default::default()
            },
        };

        run(&args_common, &args).await?;

        let compound_columns = read_vcf(&args.output)
            .await?
            .iter()
            .flat_map(|record| match record.info().get("ANN").flatten() {
                Some(Value::Array(Array::String(anns))) => anns.clone(),
                _ => vec![],
            })
            .flatten()
            .map(|ann| {
                let fields = ann.split('|').collect::<Vec<_>>();
                format!(
                    "{} {} {}",
                    fields[6],
                    fields[fields.len() - 2],
                    fields[fields.len() - 1]
                )
            })
            .collect::<Vec<_>>();
        insta::assert_yaml_snapshot!(compound_columns);

        Ok(())
    }

    #[test]
    fn test_bin_from_range() -> Result<(), anyhow::Error> {
        assert_eq!(bin_from_range(0, 0)?, 585);
        assert_eq!(bin_from_range(1, 0)?, 585);
        assert_eq!(bin_from_range(0, 1)?, 585);
        assert_eq!(bin_from_range(0, 42)?, 585);
        assert_eq!(bin_from_range(42_424_242, 42_424_243)?, 908);
        assert_eq!(bin_from_range(0, 42_424_243)?, 1);

        Ok(())
    }

    #[tokio::test]
    async fn test_badly_formed_vcf_entry() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.vcf");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db";
        let assembly = "grch37";
        let args = Args {
            threads: 1,
            reference: None,
            in_memory_reference: true,
            assembly: Some(assembly.into()),
            predictor_settings: PredictorSettings {
                transcript_settings: TranscriptSettings::most_severe(PickUnit::Gene),
                ..Default::default()
            },
            input: String::from("tests/data/annotate/seqvars/badly_formed_vcf_entry.vcf"),
            output: path_out.into_os_string().into_string().unwrap(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
                ..Default::default()
            },
        };

        run(&args_common, &args).await?;

        let records_written = read_vcf(&args.output).await?;
        let mut snapshot_data = std::collections::BTreeMap::new();

        for record in records_written {
            let key = format!(
                "{}:{}:{}:{}",
                record.reference_sequence_name(),
                record
                    .variant_start()
                    .map_or_else(|| "0".into(), |s| s.to_string()),
                record.reference_bases(),
                record.alternate_bases().as_ref().join(",")
            );

            let ann_field = record.info().get("ANN").flatten().map(|v| match v {
                Value::Array(Array::String(inner)) => inner
                    .iter()
                    .map(|s| s.clone().unwrap_or_default())
                    .collect::<Vec<_>>()
                    .join("\n"),
                _ => "".into(),
            });

            snapshot_data.insert(key, ann_field);
        }

        insta::assert_yaml_snapshot!("badly_formed_vcf_entry_output", snapshot_data);

        Ok(())
    }

    /// Mitochondrial variants called by the DRAGEN v310 germline caller have special format
    /// considerations.
    ///
    /// See: https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/DRAGEN/MitochondrialCalling.htm
    #[tokio::test]
    async fn test_dragen_mitochondrial_variant() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.vcf");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db";
        let assembly = "grch37";
        let args = Args {
            threads: 1,
            reference: None,
            in_memory_reference: true,
            assembly: Some(assembly.into()),
            predictor_settings: PredictorSettings {
                transcript_settings: TranscriptSettings::most_severe(PickUnit::Gene),
                ..Default::default()
            },
            input: String::from("tests/data/annotate/seqvars/mitochondrial_variants.vcf"), // <-- Corrected
            output: path_out.into_os_string().into_string().unwrap(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
                ..Default::default()
            },
        };

        run(&args_common, &args).await?;

        let records_written = read_vcf(&args.output).await?;
        let mut snapshot_data = std::collections::BTreeMap::new();

        for record in records_written {
            let key = format!(
                "{}:{}:{}:{}",
                record.reference_sequence_name(),
                record
                    .variant_start()
                    .map_or_else(|| "0".into(), |s| s.to_string()),
                record.reference_bases(),
                record.alternate_bases().as_ref().join(",")
            );

            let ann_field = record.info().get("ANN").flatten().map(|v| match v {
                Value::Array(Array::String(inner)) => inner
                    .iter()
                    .map(|s| s.clone().unwrap_or_default())
                    .collect::<Vec<_>>()
                    .join("\n"),
                _ => "".into(),
            });

            snapshot_data.insert(key, ann_field);
        }

        insta::assert_yaml_snapshot!("dragen_mitochondrial_variant_output", snapshot_data);

        Ok(())
    }

    /// Variants at the origin of the circular chrM shift across it during normalization.
    #[tokio::test]
    async fn test_mitochondrial_origin() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.vcf");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db/grch38";
        let args = Args {
            threads: 1,
            reference: Some(format!("{prefix}/chrM.fasta").into()),
            in_memory_reference: true,
            assembly: Some("grch38".into()),
            predictor_settings: Default::default(),
            input: String::from("tests/data/annotate/seqvars/mt_origin_variants.vcf"),
            output: path_out.into_os_string().into_string().unwrap(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            sources: Sources {
                transcripts: Some(vec![format!("{prefix}/mt-origin.txs.bin.zst")]),
                ..Default::default()
            },
        };

        run(&args_common, &args).await?;

        let mut snapshot_data = std::collections::BTreeMap::new();
        for record in read_vcf(&args.output).await? {
            let key = format!(
                "{}:{}:{}",
                record
                    .variant_start()
                    .map_or_else(|| "0".into(), |s| s.to_string()),
                record.reference_bases(),
                record.alternate_bases().as_ref().join(",")
            );
            let ann_field = match record.info().get("ANN").flatten() {
                Some(Value::Array(Array::String(inner))) => {
                    inner.iter().flatten().cloned().collect::<Vec<_>>()
                }
                _ => Vec::new(),
            };
            snapshot_data.insert(key, ann_field);
        }

        insta::assert_yaml_snapshot!("mitochondrial_origin_output", snapshot_data);

        Ok(())
    }

    /// With a reference genome, records whose REF differs from it are logged. chrM begins with
    /// GATCA, so the REF of chrM:3 matches and those of chrM:4 and chrM:5 do not.
    #[tokio::test]
    #[tracing_test::traced_test]
    async fn ref_differing_from_the_reference_is_logged() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_in = temp.join("input.vcf");
        std::fs::write(
            &path_in,
            "##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chrM\t3\t.\tT\tC\t.\t.\t.\nchrM\t4\t.\tG\tA\t.\t.\t.\nchrM\t5\t.\tC\tT\t.\t.\t.\n",
        )?;
        let prefix = "tests/data/annotate/db/grch38";
        let args = Args {
            threads: 1,
            reference: Some(format!("{prefix}/chrM.fasta").into()),
            in_memory_reference: true,
            assembly: Some("grch38".into()),
            predictor_settings: Default::default(),
            input: path_in.to_string_lossy().into(),
            output: temp.join("output.vcf").to_string_lossy().into(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            sources: Sources {
                transcripts: Some(vec![format!("{prefix}/mt-origin.txs.bin.zst")]),
                ..Default::default()
            },
        };
        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };

        run(&args_common, &args).await?;

        assert!(!logs_contain("REF at chrM:3 "));
        assert!(logs_contain(
            "REF at chrM:4 differs from the reference genome"
        ));
        assert!(logs_contain(
            "REF at chrM:5 differs from the reference genome"
        ));
        Ok(())
    }

    /// Test some variants originating from Clair3 on ONT data merged via GLNexus.
    ///
    /// Note that we currently re-use the GRCh37 databases in
    /// `tests/data/annotate/db/grch37` for GRCh38 via symlink.  As the below
    /// only is a smoke test, this is sufficient.  However, for other tests,
    /// this will pose a problem.
    #[tokio::test]
    async fn test_clair3_glnexus_variants() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.vcf");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db";
        let assembly = "grch37";
        let args = Args {
            threads: 1,
            reference: None,
            in_memory_reference: true,
            assembly: Some(assembly.into()),
            predictor_settings: PredictorSettings {
                transcript_settings: TranscriptSettings::most_severe(PickUnit::Gene),
                ..Default::default()
            },
            input: String::from("tests/data/annotate/seqvars/clair3-glnexus-min.vcf"), // <-- Corrected
            output: path_out.into_os_string().into_string().unwrap(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
                ..Default::default()
            },
        };

        run(&args_common, &args).await?;

        let records_written = read_vcf(&args.output).await?;
        let mut snapshot_data = std::collections::BTreeMap::new();

        for record in records_written {
            let key = format!(
                "{}:{}:{}:{}",
                record.reference_sequence_name(),
                record
                    .variant_start()
                    .map_or_else(|| "0".into(), |s| s.to_string()),
                record.reference_bases(),
                record.alternate_bases().as_ref().join(",")
            );

            let ann_field = record.info().get("ANN").flatten().map(|v| match v {
                Value::Array(Array::String(inner)) => inner
                    .iter()
                    .map(|s| s.clone().unwrap_or_default())
                    .collect::<Vec<_>>()
                    .join("\n"),
                _ => "".into(),
            });

            snapshot_data.insert(key, ann_field);
        }

        insta::assert_yaml_snapshot!("clair3_glnexus_variants_output", snapshot_data);

        Ok(())
    }

    /// Test corresponding to https://github.com/varfish-org/mehari/issues/409
    ///
    /// Note that we currently re-use the GRCh37 databases in
    /// `tests/data/annotate/db/grch37` for GRCh38 via symlink.  As the below
    /// only is a smoke test, this is sufficient.  However, for other tests,
    /// this will pose a problem.
    #[tokio::test]
    async fn test_brca2_zar1l_affected() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.vcf");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db";
        let assembly = "grch38";
        let args = Args {
            threads: 1,
            reference: None,
            in_memory_reference: true,
            assembly: Some(assembly.into()),
            predictor_settings: PredictorSettings {
                transcript_settings: TranscriptSettings::most_severe(PickUnit::Gene),
                ..Default::default()
            },
            input: String::from("tests/data/annotate/seqvars/brca2_zar1l/brca2_zar1l.vcf"),
            output: path_out.into_os_string().into_string().unwrap(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
                ..Default::default()
            },
        };

        run(&args_common, &args).await?;

        let records_written = read_vcf(&args.output).await?;
        let mut snapshot_data = std::collections::BTreeMap::new();

        for record in records_written {
            let key = format!(
                "{}:{}:{}:{}",
                record.reference_sequence_name(),
                record
                    .variant_start()
                    .map_or_else(|| "0".into(), |s| s.to_string()),
                record.reference_bases(),
                record.alternate_bases().as_ref().join(",")
            );

            let ann_field = record.info().get("ANN").flatten().map(|v| match v {
                Value::Array(Array::String(inner)) => inner
                    .iter()
                    .map(|s| s.clone().unwrap_or_default())
                    .collect::<Vec<_>>()
                    .join("\n"),
                _ => "".into(),
            });

            snapshot_data.insert(key, ann_field);
        }

        insta::assert_yaml_snapshot!("brca2_zar1l_affected_output", snapshot_data);

        Ok(())
    }

    /// Push `(chrom, pos, ref, alt)` variants without genotypes on one transcript
    /// into the compound buffer and flush it. `PhasingStrategy::Ignore` makes them one group.
    fn compound_processor_flush(
        annotator: &Arc<Annotator>,
        variants: &[(&str, i32, &str, &str)],
    ) -> Result<Vec<AnnotatedVariant>, anyhow::Error> {
        let compound_settings = CompoundSettings {
            enable_compound_variants: true,
            phasing_strategy: PhasingStrategy::Ignore,
            ..Default::default()
        };
        let mut processor = VariantProcessor::new(annotator.clone(), &compound_settings);
        for &(chrom, pos, reference, alternative) in variants {
            let record = VcfRecord::builder()
                .set_reference_sequence_name(chrom)
                .set_variant_start(Position::try_from(pos as usize)?)
                .set_reference_bases(reference)
                .set_alternate_bases(AlternateBases::from(vec![alternative.to_string()]))
                .build();
            let vcf_var = VcfVariant {
                chromosome: chrom.into(),
                position: pos,
                reference: reference.into(),
                alternative: alternative.into(),
            };
            let tx_accessions = HashSet::from(["NM_007294.3".to_string()]);
            processor
                .buffer
                .push(vcf_var, record, tx_accessions, 0, 100_000_000, 0, 1);
        }
        processor.flush()
    }

    fn compound_annotator() -> Result<Arc<Annotator>, anyhow::Error> {
        let predictor_settings = PredictorSettings {
            compound_settings: CompoundSettings {
                enable_compound_variants: true,
                ..Default::default()
            },
            ..Default::default()
        };
        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let consequence = ConsequenceAnnotator::from_db_and_settings(
            tx_db,
            None::<PathBuf>,
            false,
            &predictor_settings,
        )?;
        Ok(Arc::new(Annotator::new(vec![AnnotatorEnum::Consequence(
            consequence,
        )])))
    }

    #[test]
    fn compound_prediction_error_is_returned() -> Result<(), anyhow::Error> {
        let annotator = compound_annotator()?;
        // "chrFake" has no accession, so `predict_multiple` fails with an
        // error that is not a group validation error.
        let result = compound_processor_flush(
            &annotator,
            &[("chrFake", 100, "A", "C"), ("chrFake", 200, "A", "C")],
        );

        let msg = format!("{:#}", result.expect_err("expected an error"));
        assert!(msg.contains("chrFake:100:A:C,chrFake:200:A:C"), "{msg}");
        assert!(
            msg.contains("Could not determine chromosome accession"),
            "{msg}"
        );
        Ok(())
    }

    #[test]
    #[tracing_test::traced_test]
    fn compound_group_validation_error_is_logged() -> Result<(), anyhow::Error> {
        let annotator = compound_annotator()?;
        // The deletion covers the SNV, so the group fails validation.
        let annotated = compound_processor_flush(
            &annotator,
            &[("17", 41197700, "AG", "A"), ("17", 41197701, "G", "C")],
        )?;

        assert_eq!(annotated.len(), 2);
        assert!(logs_contain("Overlapping variants detected"));
        Ok(())
    }
}
