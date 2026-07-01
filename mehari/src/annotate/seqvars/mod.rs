//! Annotation of sequence variants.
use std::env;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::pin::Pin;
use std::sync::Arc;
use std::time::Instant;

use self::consequence::terms::AnnField;
use crate::annotate::cli::{PredictorSettings, Sources};
use crate::annotate::seqvars::consequence::ConfigBuilder;
use crate::annotate::seqvars::consequence::logic::ConsequencePredictor;
use crate::annotate::seqvars::consequence::terms::{
    ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_COMPOUND_IDS, ANN_COMPOUND_VARIANTS, ANN_TX_SEQ_ALT,
    ANN_TX_SEQ_REF,
};
use crate::common::contig::ContigManager;
use crate::common::noodles::{NoodlesVariantReader, open_variant_reader, open_variant_writer};
use crate::db::keys::Var;
use crate::db::transcripts::merge::merge_transcript_databases;
use crate::pbs;
use crate::pbs::txs::TxSeqDatabase;
use anyhow::Error;
use biocommons_bioutils::assemblies::Assembly;
use clap::Parser;
use clinvar::{ClinvarAnnotator, ClinvarResult};
use consequence::{Config, ConsequenceAnnotator, VcfVariant};
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
pub mod custom;
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

    /// Reference genome FASTA file (with accompanying index).
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
    /// A sweet spot regarding trade-off between I/O and CPU-bound tasks is around 5.
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
    csq_config: &Config,
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
        let fields = AnnField::ann_field_names(csq_config).join(" | ");
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

    fn set_csq_config(&mut self, _config: Config) {}
}

pub(crate) fn prepare_vcf_record(
    record: AnnotatedVariant,
    csq_config: &Config,
) -> noodles::vcf::variant::RecordBuf {
    let mut out_record = record.vcf;

    if let Some(freqs) = &record.annotation.frequencies {
        let infos = out_record.info_mut();
        infos.insert(
            "gnomad_exomes_an".into(),
            Some(field::Value::Integer(freqs.gnomad_exomes_an)),
        );
        infos.insert(
            "gnomad_exomes_hom".into(),
            Some(field::Value::Integer(freqs.gnomad_exomes_hom)),
        );
        infos.insert(
            "gnomad_exomes_het".into(),
            Some(field::Value::Integer(freqs.gnomad_exomes_het)),
        );
        if let Some(hemi) = freqs.gnomad_exomes_hemi {
            infos.insert(
                "gnomad_exomes_hemi".into(),
                Some(field::Value::Integer(hemi)),
            );
        }

        infos.insert(
            "gnomad_genomes_an".into(),
            Some(field::Value::Integer(freqs.gnomad_genomes_an)),
        );
        infos.insert(
            "gnomad_genomes_hom".into(),
            Some(field::Value::Integer(freqs.gnomad_genomes_hom)),
        );
        infos.insert(
            "gnomad_genomes_het".into(),
            Some(field::Value::Integer(freqs.gnomad_genomes_het)),
        );
        if let Some(hemi) = freqs.gnomad_genomes_hemi {
            infos.insert(
                "gnomad_genomes_hemi".into(),
                Some(field::Value::Integer(hemi)),
            );
        }

        if let Some(an) = freqs.helix_an {
            infos.insert("helix_an".into(), Some(field::Value::Integer(an)));
        }
        if let Some(hom) = freqs.helix_hom {
            infos.insert("helix_hom".into(), Some(field::Value::Integer(hom)));
        }
        if let Some(het) = freqs.helix_het {
            infos.insert("helix_het".into(), Some(field::Value::Integer(het)));
        }
    }

    if let Some(clinvar) = &record.annotation.clinvar {
        let infos = out_record.info_mut();
        infos.insert(
            "clinvar_vcv".into(),
            Some(field::Value::Array(field::value::Array::String(
                clinvar.vcv.iter().map(|s| Some(s.clone())).collect(),
            ))),
        );
        infos.insert(
            "clinvar_germline_classification".into(),
            Some(field::Value::Array(field::value::Array::String(
                clinvar
                    .germline_classification
                    .iter()
                    .map(|s| Some(s.clone()))
                    .collect(),
            ))),
        );
    }

    if let Some(cadd) = &record.annotation.cadd {
        let infos = out_record.info_mut();
        infos.insert(
            "CADD_RAW".into(),
            Some(field::Value::Array(field::value::Array::Float(vec![Some(
                cadd.raw_score,
            )]))),
        );
        infos.insert(
            "CADD_PHRED".into(),
            Some(field::Value::Array(field::value::Array::Float(vec![Some(
                cadd.phred,
            )]))),
        );
    }

    if let Some(spliceai) = &record.annotation.spliceai {
        let infos = out_record.info_mut();
        let mut symbols = Vec::new();
        let mut ds_ags = Vec::new();
        let mut ds_als = Vec::new();
        let mut ds_dgs = Vec::new();
        let mut ds_dls = Vec::new();
        let mut dp_ags = Vec::new();
        let mut dp_als = Vec::new();
        let mut dp_dgs = Vec::new();
        let mut dp_dls = Vec::new();

        for pred in &spliceai.predictions {
            symbols.push(pred.symbol.clone());
            ds_ags.push(format!("{:.2}", pred.ds_ag));
            ds_als.push(format!("{:.2}", pred.ds_al));
            ds_dgs.push(format!("{:.2}", pred.ds_dg));
            ds_dls.push(format!("{:.2}", pred.ds_dl));
            dp_ags.push(format!("{}", pred.dp_ag));
            dp_als.push(format!("{}", pred.dp_al));
            dp_dgs.push(format!("{}", pred.dp_dg));
            dp_dls.push(format!("{}", pred.dp_dl));
        }

        let symbol_str = symbols.join("&");
        let ds_ag_str = ds_ags.join("&");
        let ds_al_str = ds_als.join("&");
        let ds_dg_str = ds_dgs.join("&");
        let ds_dl_str = ds_dls.join("&");
        let dp_ag_str = dp_ags.join("&");
        let dp_al_str = dp_als.join("&");
        let dp_dg_str = dp_dgs.join("&");
        let dp_dl_str = dp_dls.join("&");

        infos.insert(
            "SpliceAI_pred_SYMBOL".into(),
            Some(field::Value::Array(field::value::Array::String(vec![
                Some(symbol_str),
            ]))),
        );
        infos.insert(
            "SpliceAI_pred_DS_AG".into(),
            Some(field::Value::Array(field::value::Array::String(vec![
                Some(ds_ag_str),
            ]))),
        );
        infos.insert(
            "SpliceAI_pred_DS_AL".into(),
            Some(field::Value::Array(field::value::Array::String(vec![
                Some(ds_al_str),
            ]))),
        );
        infos.insert(
            "SpliceAI_pred_DS_DG".into(),
            Some(field::Value::Array(field::value::Array::String(vec![
                Some(ds_dg_str),
            ]))),
        );
        infos.insert(
            "SpliceAI_pred_DS_DL".into(),
            Some(field::Value::Array(field::value::Array::String(vec![
                Some(ds_dl_str),
            ]))),
        );
        infos.insert(
            "SpliceAI_pred_DP_AG".into(),
            Some(field::Value::Array(field::value::Array::String(vec![
                Some(dp_ag_str),
            ]))),
        );
        infos.insert(
            "SpliceAI_pred_DP_AL".into(),
            Some(field::Value::Array(field::value::Array::String(vec![
                Some(dp_al_str),
            ]))),
        );
        infos.insert(
            "SpliceAI_pred_DP_DG".into(),
            Some(field::Value::Array(field::value::Array::String(vec![
                Some(dp_dg_str),
            ]))),
        );
        infos.insert(
            "SpliceAI_pred_DP_DL".into(),
            Some(field::Value::Array(field::value::Array::String(vec![
                Some(dp_dl_str),
            ]))),
        );
    }

    if let Some(dbsnp) = &record.annotation.dbsnp {
        out_record.ids_mut().as_mut().insert(dbsnp.rs_id.clone());
    }

    for custom_db in &record.annotation.custom {
        let infos = out_record.info_mut();
        for (key, val) in &custom_db.fields {
            let field_name = format!("{}_{}", custom_db.db_name, key);
            infos.insert(
                field_name,
                Some(field::Value::Array(field::value::Array::String(vec![
                    Some(val.clone()),
                ]))),
            );
        }
    }

    if !record.annotation.consequences.is_empty() {
        let formatted_anns: Vec<Option<String>> = record
            .annotation
            .consequences
            .iter()
            .map(|ann| Some(ann.format(csq_config)))
            .collect();
        out_record.info_mut().insert(
            "ANN".into(),
            Some(field::Value::Array(field::value::Array::String(
                formatted_anns,
            ))),
        );
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
    pub csq_config: Config,
}

impl SeqvarsVcfWriter {
    pub fn new(inner: crate::common::noodles::VariantWriter) -> Self {
        Self {
            inner,
            csq_config: Config::default(),
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
        let out_record = prepare_vcf_record(record, &self.csq_config);

        self.inner
            .write_annotated_record(
                header,
                AnnotatedVariant {
                    vcf: out_record,
                    annotation: VariantAnnotation {
                        consequences: vec![],
                        frequencies: None,
                        clinvar: None,
                        cadd: None,
                        spliceai: None,
                        dbsnp: None,
                        custom: vec![],
                    },
                },
            )
            .await
    }

    async fn shutdown(&mut self) -> Result<(), anyhow::Error> {
        self.inner.shutdown().await
    }

    fn set_csq_config(&mut self, config: Config) {
        self.csq_config = config;
    }
}

pub struct SeqvarJsonlWriter {
    inner: Pin<Box<dyn AsyncWrite>>,
    csq_config: Config,
}

impl SeqvarJsonlWriter {
    pub fn new(inner: Pin<Box<dyn AsyncWrite>>) -> Self {
        Self {
            inner,
            csq_config: Config::default(),
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
    fn set_csq_config(&mut self, config: Config) {
        self.csq_config = config;
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
    Custom(String, custom::CustomDbAnnotator),
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
            AnnotatorEnum::Custom(name, a) => {
                if let Some(custom_res) = a.annotate(name, var)? {
                    annotation.custom.push(custom_res);
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
                    ("SYMBOL", "gene symbol"),
                    ("DS_AG", "delta score for acceptor gain"),
                    ("DS_AL", "delta score for acceptor loss"),
                    ("DS_DG", "delta score for donor gain"),
                    ("DS_DL", "delta score for donor loss"),
                    ("DP_AG", "delta position for acceptor gain"),
                    ("DP_AL", "delta position for acceptor loss"),
                    ("DP_DG", "delta position for donor gain"),
                    ("DP_DL", "delta position for donor loss"),
                ];
                for (name, desc) in fields {
                    let field_name = format!("SpliceAI_pred_{}", name);
                    header.infos_mut().insert(
                        field_name,
                        Map::<Info>::new(
                            Number::Unknown,
                            InfoType::String,
                            format!("SpliceAI (v1.3.1) prediction. {}", desc),
                        ),
                    );
                }
                Ok(())
            }
            AnnotatorEnum::Dbsnp(_) => {
                // dbsnp just inserts IDs into the ID column
                Ok(())
            }
            AnnotatorEnum::Custom(name, a) => {
                let fields = a.get_fields()?;
                for f in fields {
                    let field_name = format!("{}_{}", name, f);
                    header.infos_mut().insert(
                        field_name,
                        Map::<Info>::new(
                            Number::Unknown,
                            InfoType::String,
                            format!("Custom lookup annotation from database {}", name),
                        ),
                    );
                }
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
            custom: vec![],
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

    let mut custom_columns = Vec::new();
    let c = &args.predictor_settings.reporting_settings;
    if c.report_cdna_sequence.includes_ref() {
        custom_columns.push(ANN_TX_SEQ_REF.to_string());
    }
    if c.report_cdna_sequence.includes_alt() {
        custom_columns.push(ANN_TX_SEQ_ALT.to_string());
    }
    if c.report_protein_sequence.includes_ref() {
        custom_columns.push(ANN_AA_SEQ_REF.to_string());
    }
    if c.report_protein_sequence.includes_alt() {
        custom_columns.push(ANN_AA_SEQ_ALT.to_string());
    }

    let enable_compound = args
        .predictor_settings
        .compound_settings
        .enable_compound_variants;

    if enable_compound {
        custom_columns.push(ANN_COMPOUND_IDS.to_string());
        custom_columns.push(ANN_COMPOUND_VARIANTS.to_string());
    }

    // TODO: manually rebuilding Config here so we can automatically build the VCF ANN header
    //   is not the best way of doing things.
    let csq_config = ConfigBuilder::default()
        .report_cdna_sequence(c.report_cdna_sequence)
        .report_protein_sequence(c.report_protein_sequence)
        .custom_columns(custom_columns)
        .build()?;

    writer.set_csq_config(csq_config.clone());

    let header_out = build_header(
        &header_in,
        with_annotations,
        with_frequencies,
        with_clinvar,
        &annotator.annotators,
        &additional_header_info,
        &csq_config,
    );

    // Perform the VCF annotation.
    let mut processor =
        VariantProcessor::new(&annotator, &args.predictor_settings.compound_settings);

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
                    next_batch.push(vcf_record);
                    total_read += 1;
                }
                _ => break, // EOF
            }
        }

        // await the previous batch's results and write them out sequentially
        if let Some(handle) = processing_handle.take() {
            let processed_batch = handle.await??;

            for annotated_record in processed_batch {
                if prev.elapsed().as_secs() >= 60 {
                    tracing::info!("at {:?}", from_vcf_allele(&annotated_record.vcf, 0));
                    prev = Instant::now();
                }

                if enable_compound {
                    let ready_records = processor.process_annotated_record(
                        annotated_record.vcf,
                        annotated_record.annotation,
                    )?;
                    for rec in ready_records {
                        writer.write_annotated_record(&header_out, rec).await?;
                    }
                } else {
                    writer
                        .write_annotated_record(&header_out, annotated_record)
                        .await?;
                }
                total_written += 1;
            }
        }

        // if there are no more records to process, stop.
        if next_batch.is_empty() {
            break;
        }

        // process `next_batch` in parallel
        let batch_to_process = std::mem::replace(&mut next_batch, Vec::with_capacity(batch_size));
        let annotator_clone = annotator.clone();

        processing_handle = Some(tokio::task::spawn_blocking(move || {
            use rayon::prelude::*;

            batch_to_process
                .into_par_iter()
                .map(|record| {
                    let annotation = annotator_clone.annotate(&record)?;
                    Ok(AnnotatedVariant {
                        vcf: record,
                        annotation,
                    })
                })
                .collect::<anyhow::Result<Vec<_>>>()
        }));
    }

    if enable_compound {
        let final_records = processor.flush()?;
        for rec in final_records {
            writer.write_annotated_record(&header_out, rec).await?;
        }
    }

    tracing::info!(
        "... annotated {} records in {:?}",
        total_written.separate_with_commas(),
        start.elapsed()
    );
    writer.shutdown().await?;
    Ok(())
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

    // Add the custom DB annotator if requested.
    if let Some(custom_db_specs) = &sources.custom_db {
        let custom_dbs = initialize_custom_db_annotators_for_assembly(
            custom_db_specs,
            &assembly,
            contig_manager.clone(),
        )?;
        for (name, custom_db) in custom_dbs {
            annotators.push(AnnotatorEnum::Custom(name, custom_db))
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

pub fn initialize_custom_db_annotators_for_assembly(
    custom_db_specs: &[String],
    assembly: &str,
    contig_manager: Arc<ContigManager>,
) -> Result<Vec<(String, custom::CustomDbAnnotator)>, Error> {
    let mut annotators = Vec::new();

    for spec in custom_db_specs {
        let parts: Vec<&str> = spec.splitn(2, '=').collect();
        if parts.len() != 2 {
            anyhow::bail!(
                "Invalid custom DB specification: '{}'. Must be in name=path format.",
                spec
            );
        }
        let name = parts[0].to_string();
        let rocksdb_path = parts[1];

        verify_database_assembly(rocksdb_path, assembly);
        tracing::info!(
            "Loading custom database '{}' for assembly {:?} from {}",
            name,
            assembly,
            rocksdb_path
        );
        let annotator = custom::CustomDbAnnotator::from_path(rocksdb_path, contig_manager.clone())?;
        annotators.push((name, annotator));
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
    pub custom: Vec<custom::CustomDbResult>,
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

struct VariantProcessor<'a> {
    annotator: &'a Annotator,
    buffer: crate::annotate::seqvars::compound::VariantBuffer,
    next_group_id: usize,
    annotations: std::collections::HashMap<String, VariantAnnotation>,
}

impl<'a> VariantProcessor<'a> {
    pub fn new(
        annotator: &'a Annotator,
        compound_settings: &crate::annotate::cli::CompoundSettings,
    ) -> Self {
        Self {
            annotator,
            buffer: compound::VariantBuffer::new(compound_settings.phasing_strategy),
            next_group_id: 0,
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
        let key = format!(
            "{}:{}:{}:{}",
            vcf_var.chrom, vcf_var.pos, vcf_var.reference, vcf_var.alternative
        );
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
                    &vcf_var.chrom,
                    vcf_var.pos,
                    vcf_var.reference.len(),
                )
            } else {
                (std::collections::HashSet::new(), i32::MAX, -1)
            };

        let csq_var = VcfVariant {
            chromosome: vcf_var.chrom,
            position: vcf_var.pos,
            reference: vcf_var.reference,
            alternative: vcf_var.alternative,
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
                let vcf_var = from_vcf_allele(&b.record, 0);
                let key = format!(
                    "{}:{}:{}:{}",
                    vcf_var.chrom, vcf_var.pos, vcf_var.reference, vcf_var.alternative
                );
                self.annotations
                    .remove(&key)
                    .unwrap_or_else(|| VariantAnnotation {
                        consequences: vec![],
                        frequencies: None,
                        clinvar: None,
                        cadd: None,
                        spliceai: None,
                        dbsnp: None,
                        custom: vec![],
                    })
            })
            .collect();

        for group_indices in compound_groups {
            if let Some(predictor) = predictor {
                let vcf_vars: Vec<_> = group_indices
                    .iter()
                    .map(|&i| variants[i].vcf_var.clone())
                    .collect();

                if let Ok(Some(mut ann_fields)) = predictor.predict_multiple(&vcf_vars) {
                    if ann_fields.is_empty() {
                        continue;
                    }

                    self.next_group_id += 1;
                    let group_id_str = format!("comp_{}", self.next_group_id);
                    let constituent_vars_str = vcf_vars
                        .iter()
                        .map(|v| {
                            format!(
                                "{}:{}:{}:{}",
                                v.chromosome, v.position, v.reference, v.alternative
                            )
                        })
                        .collect::<Vec<_>>()
                        .join(",");

                    for ann in &mut ann_fields {
                        ann.custom_fields
                            .insert(ANN_COMPOUND_IDS.to_string(), Some(group_id_str.clone()));
                        ann.custom_fields.insert(
                            ANN_COMPOUND_VARIANTS.to_string(),
                            Some(constituent_vars_str.clone()),
                        );
                    }

                    for &i in &group_indices {
                        result_annotations[i]
                            .consequences
                            .extend(ann_fields.clone());
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
    use super::{Args, OutputFormat, run};
    use crate::annotate::cli::{ConsequenceBy, PredictorSettings};
    use crate::annotate::cli::{Sources, TranscriptSettings};
    use crate::common::noodles::{NoodlesVariantReader, open_variant_reader};
    use clap_verbosity_flag::Verbosity;
    use futures::TryStreamExt;
    use pretty_assertions::assert_eq;
    use std::path::Path;
    use temp_testdir::TempDir;

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
                transcript_settings: TranscriptSettings {
                    report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                    ..Default::default()
                },
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
                transcript_settings: TranscriptSettings {
                    report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                    ..Default::default()
                },
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
                transcript_settings: TranscriptSettings {
                    report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                    ..Default::default()
                },
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
                transcript_settings: TranscriptSettings {
                    report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                    ..Default::default()
                },
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
                transcript_settings: TranscriptSettings {
                    report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                    ..Default::default()
                },
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
}
