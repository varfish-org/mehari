//! Annotation of sequence variants.
use std::env;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::pin::Pin;
use std::sync::Arc;
use std::time::Instant;

use self::consequence::terms::{ANN_KEY, Allele, AnnField, Message};
use crate::annotate::cli::{PredictorSettings, Sources};
use crate::annotate::seqvars::consequence::logic::{ConsequencePredictor, is_sequence_allele};
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
use consequence::event::BreakendEvent;
use consequence::{
    ConfidenceInterval, ConsequenceAnnotator, SvClaim, SvInfo, Variant, VariantKind, VcfVariant,
    round_copy_number,
};
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
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use strum::Display;
use thousands::Separable;
use tokio::io::{AsyncWrite, AsyncWriteExt};

pub mod binning;
pub mod cadd;
mod compound;
mod mates;
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

/// Accepts an INFO key that matches the VCF specification's
/// `^([A-Za-z_][0-9A-Za-z_.]*|1000G)$`.
fn parse_ann_key(key: &str) -> Result<String, String> {
    let mut chars = key.chars();
    let valid = key == "1000G"
        || (chars
            .next()
            .is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
            && chars.all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.')));
    if valid {
        Ok(key.to_string())
    } else {
        Err("an INFO key must match ^([A-Za-z_][0-9A-Za-z_.]*|1000G)$".to_string())
    }
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

    /// INFO key to write the annotations to.
    #[arg(long, default_value = ANN_KEY, value_parser = parse_ann_key)]
    pub ann_key: String,

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
    ann_key: Option<&str>,
    with_frequencies: bool,
    with_clinvar: bool,
    annotators: &[AnnotatorEnum],
    additional_records: &[(String, String)],
    custom_columns: &[String],
) -> Result<VcfHeader, anyhow::Error> {
    if let Some(ann_key) = ann_key {
        // The annotations would replace the values of such a field, or its header their names.
        let other_fields = build_header(
            &VcfHeader::default(),
            None,
            with_frequencies,
            with_clinvar,
            annotators,
            &[],
            &[],
        )?;
        anyhow::ensure!(
            !other_fields.infos().contains_key(ann_key),
            "--ann-key {ann_key} is also an INFO field that this run writes"
        );
    }

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

    if let Some(ann_key) = ann_key {
        let fields = AnnField::ann_field_names(custom_columns).join(" | ");
        // `Source` and `Version` are the VCF specification's attributes to identify the producer
        // of an INFO field; consumers can tell `ANN` from mehari apart from SnpEff's.
        let mut ann = Map::<Info>::new(
            Number::Unknown,
            InfoType::String,
            format!("Functional annotations: '{fields}'"),
        );
        ann.other_fields_mut()
            .insert("Source".parse()?, "mehari".into());
        ann.other_fields_mut()
            .insert("Version".parse()?, env!("CARGO_PKG_VERSION").into());
        header_out.infos_mut().insert(ann_key.into(), ann);
        header_out.infos_mut().insert(
            "FUSION".into(),
            Map::<Info>::new(
                Number::Unknown,
                InfoType::String,
                "Gene fusions of a structural variant or breakend as GENE1::GENE2: the 5' gene \
                 first, or for a bidirectional fusion the gene before the junction",
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

    Ok(header_out)
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

    fn set_ann_key(&mut self, _ann_key: String) {}

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
    ann_key: &str,
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
        ($key:expr, array $array_variant:ident, $iter:expr) => {
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

    if !record.annotation.fusions.is_empty() {
        insert_field!("FUSION", array String, record.annotation.fusions.iter());
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
        insert_field!(ann_key, array String, formatted_anns);
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
    pub ann_key: String,
}

impl SeqvarsVcfWriter {
    pub fn new(inner: crate::common::noodles::VariantWriter) -> Self {
        Self {
            inner,
            custom_columns: Vec::new(),
            ann_key: ANN_KEY.to_string(),
        }
    }

    /// Moves the annotation into the VCF record; preparing the result again changes nothing.
    fn prepare(
        record: AnnotatedVariant,
        custom_columns: &[String],
        ann_key: &str,
    ) -> AnnotatedVariant {
        AnnotatedVariant {
            vcf: prepare_vcf_record(record, custom_columns, ann_key),
            annotation: VariantAnnotation {
                consequences: vec![],
                fusions: vec![],
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
            .write_annotated_record(
                header,
                Self::prepare(record, &self.custom_columns, &self.ann_key),
            )
            .await
    }

    async fn shutdown(&mut self) -> Result<(), anyhow::Error> {
        self.inner.shutdown().await
    }

    fn set_custom_columns(&mut self, custom_columns: Vec<String>) {
        self.custom_columns = custom_columns;
    }

    fn set_ann_key(&mut self, ann_key: String) {
        self.ann_key = ann_key;
    }

    fn record_preparer(
        &self,
    ) -> impl Fn(AnnotatedVariant) -> AnnotatedVariant + Send + Sync + 'static {
        let custom_columns = self.custom_columns.clone();
        let ann_key = self.ann_key.clone();
        move |record| Self::prepare(record, &custom_columns, &ann_key)
    }
}

pub struct SeqvarJsonlWriter {
    inner: Pin<Box<dyn AsyncWrite>>,
    custom_columns: Vec<String>,
    ann_key: String,
}

impl SeqvarJsonlWriter {
    pub fn new(inner: Pin<Box<dyn AsyncWrite>>) -> Self {
        Self {
            inner,
            custom_columns: Vec::new(),
            ann_key: ANN_KEY.to_string(),
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
            if key == &self.ann_key {
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
                "fusions": record.annotation.fusions,
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

    fn set_ann_key(&mut self, ann_key: String) {
        self.ann_key = ann_key;
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

    /// The annotation of `vcf_record` with its allele `variant`, see `classify`. The last record
    /// of a group of breakends also gets the rows and fusions of `event`, the event of its group,
    /// see `mates::MateGrouper`.
    pub fn annotate(
        &self,
        vcf_record: &VcfRecord,
        variant: &Classified,
        event: Option<&BreakendEvent>,
    ) -> anyhow::Result<VariantAnnotation> {
        let vcf_var = from_vcf_allele(vcf_record, 0);
        let mut annotation = VariantAnnotation {
            consequences: vec![],
            fusions: vec![],
            frequencies: None,
            clinvar: None,
            cadd: None,
            spliceai: None,
            dbsnp: None,
        };

        // A structural variant or a breakend gets consequences only: the databases hold no
        // symbolic alleles.
        if let Ok(Some(variant)) = variant
            && !matches!(variant.kind, VariantKind::Sequence)
        {
            if let Some(predictor) = self.consequence_predictor() {
                let annotations = predictor.predict_variant(variant)?.unwrap_or_default();
                annotation.consequences = annotations.anns;
                annotation.fusions = annotations.fusions;
                // Every row of a record has its allele, also a row of its event.
                if let Some(event) = event {
                    let event = predictor.predict_event(event).map_err(|e| {
                        anyhow::anyhow!(
                            "failed to annotate the event of breakends of {}: {e}",
                            variant.var
                        )
                    })?;
                    let allele = Allele::Alt {
                        alternative: variant.var.alternative.clone(),
                    };
                    let rows = (event.iter().flat_map(|event| &event.anns)).map(|ann| AnnField {
                        allele: allele.clone(),
                        ..ann.clone()
                    });
                    let settings = &predictor.config.transcript_settings;
                    let pick = settings.pick.is_some() && !settings.flag_pick;
                    add_event_rows(&mut annotation.consequences, rows, pick);
                    (annotation.fusions)
                        .extend(event.iter().flat_map(|event| event.fusions.iter().cloned()));
                    annotation.fusions.sort();
                    annotation.fusions.dedup();
                }
            }
            return Ok(annotation);
        }
        if !is_annotatable(vcf_record) {
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

/// Adds `rows`, those of the event of breakends, to `anns`, those of the last record of its group.
/// With `pick`, both were picked, and a transcript keeps one row: that of higher impact, the
/// record's own on a tie.
fn add_event_rows(anns: &mut Vec<AnnField>, rows: impl Iterator<Item = AnnField>, pick: bool) {
    for row in rows {
        let own = pick
            .then(|| anns.iter().position(|ann| ann.feature_id == row.feature_id))
            .flatten();
        match own {
            None => anns.push(row),
            Some(i) if row.putative_impact < anns[i].putative_impact => anns[i] = row,
            Some(_) => {}
        }
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
    warn_on_hg19_chrm(&header_in, reference_provider.as_deref());

    let mut additional_header_info = annotator.versions_for_vcf_header();
    additional_header_info.push(("mehariCmd".into(), env::args().join(" ")));
    additional_header_info.push(("mehariVersion".into(), env!("CARGO_PKG_VERSION").into()));

    let with_annotations = args
        .sources
        .transcripts
        .as_ref()
        .is_some_and(|v| !v.is_empty());
    if with_annotations && header_in.infos().contains_key(&args.ann_key) {
        tracing::warn!(
            "The input already has the INFO field {0}. It is removed from all records, \
            and {0} in the output only holds mehari's annotations.",
            args.ann_key
        );
    }
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
    writer.set_ann_key(args.ann_key.clone());

    let header_out = build_header(
        &header_in,
        with_annotations.then_some(args.ann_key.as_str()),
        with_frequencies,
        with_clinvar,
        &annotator.annotators,
        &additional_header_info,
        &custom_columns,
    )?;

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
    let mut total_unsupported = 0usize;
    let mut total_unknown_claims = 0usize;
    // Holds the last record of a breakend group back until the group is closed.
    let mut grouper = mates::MateGrouper::new(header_in.contigs().keys().cloned());

    writer.write_noodles_header(&header_out).await?;

    use futures::TryStreamExt;
    let mut records = reader.records(&header_in).await;
    let file_format = header_in.file_format();

    let worker_threads = rayon::current_num_threads();
    let batch_size = worker_threads * 1024;

    // Preparing records for output is costly, so do it on the worker threads, unless there is
    // only one: then the writer has time to spare.  Compound variant detection still needs the
    // annotations, so the writer prepares those records itself.
    let preparer =
        (!enable_compound && worker_threads > 1).then(|| Arc::new(writer.record_preparer()));

    let mut next_batch = Vec::with_capacity(batch_size);
    // Each record with its allele, see `classify`.
    type Batch = Vec<(AnnotatedVariant, Classified)>;
    let mut processing_handle: Option<tokio::task::JoinHandle<anyhow::Result<Batch>>> = None;

    loop {
        // fill the `next_batch` buffer asynchronously (happens concurrently while the previous batch is being processed)
        while next_batch.len() < batch_size {
            if let Some(max) = args.max_var_count
                && total_read >= max
            {
                next_batch.extend(grouper.finish());
                break;
            }

            match records.try_next().await? {
                Some(vcf_record) => {
                    // A record without an alternate allele (`ALT` `.`) passes through.
                    if vcf_record.alternate_bases().len() > 1 {
                        tracing::error!(
                            "Found record with more than one alternate allele.  This is currently not supported. \
                    Please split multi-allelic records with `bcftools norm -m-any`, without `-f`, \
                    which rejects `N` as `REF`.  Record: {:?}",
                            &vcf_record
                        );
                        anyhow::bail!("multi-allelic records not supported");
                    }
                    let variant = classify(&vcf_record, file_format);
                    if is_event(&vcf_record) && variant.is_err() {
                        total_unsupported += 1;
                    }
                    if is_event(&vcf_record)
                        && compound::string(&vcf_record, "SVCLAIM")
                            .is_some_and(|claim| claim.parse::<SvClaim>().is_err())
                    {
                        total_unknown_claims += 1;
                    }
                    if is_annotatable(&vcf_record)
                        && let Some(provider) = &reference_provider
                        && ref_differs(provider, &vcf_record) == Some(true)
                    {
                        tracing::warn!(
                            "REF at {}:{} differs from the reference genome; the annotation uses \
                             the reference",
                            vcf_record.reference_sequence_name(),
                            vcf_record.variant_start().map_or(0, |pos| pos.get())
                        );
                    }
                    next_batch.extend(grouper.push(vcf_record, variant));
                    total_read += 1;
                }
                // EOF
                _ => {
                    next_batch.extend(grouper.finish());
                    break;
                }
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
            let ann_key = args.ann_key.clone();

            processing_handle = Some(tokio::task::spawn_blocking(move || {
                use rayon::prelude::*;

                batch_to_process
                    .into_par_iter()
                    .map(|(mut record, variant, event)| {
                        let annotation =
                            annotator_clone.annotate(&record, &variant, event.as_ref())?;
                        if with_annotations {
                            // The header now describes mehari's annotations under this key, so
                            // those from the input must go, also for records that get no
                            // consequences.
                            record.info_mut().as_mut().shift_remove(&ann_key);
                        }
                        let annotated = AnnotatedVariant {
                            vcf: record,
                            annotation,
                        };
                        let annotated = match &preparer {
                            Some(prepare) => prepare(annotated),
                            None => annotated,
                        };
                        Ok((annotated, variant))
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
                            for (annotated_record, variant) in processed_batch {
                                ready_records.extend(batch_processor.process_annotated_record(
                                    annotated_record.vcf,
                                    annotated_record.annotation,
                                    &variant,
                                )?);
                            }
                            Ok((ready_records, batch_processor))
                        })
                        .await??;
                    processor = Some(batch_processor);
                    ready_records
                }
                None => (processed_batch.into_iter())
                    .map(|(annotated, _)| annotated)
                    .collect(),
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
    if total_unsupported > 0 {
        tracing::warn!(
            "{} record(s) stay without annotations: malformed breakends, other records at \
             position 0, symbolic alleles without a length and unknown symbolic alleles are not \
             supported",
            total_unsupported.separate_with_commas()
        );
    }
    if total_unknown_claims > 0 {
        tracing::warn!(
            "{} record(s) have an SVCLAIM other than D, J and DJ, which counts as none",
            total_unknown_claims.separate_with_commas()
        );
    }
    if grouper.unresolved > 0 {
        tracing::warn!(
            "{} breakend record(s) whose MATEID or EVENT mates never arrived; they keep their own \
             annotations",
            grouper.unresolved.separate_with_commas()
        );
    }
    writer.shutdown().await?;
    Ok(())
}

/// The length of chrM in hg19, the Yoruba sequence NC_001807.4. GRCh37, GRCh38 and the
/// transcript databases use the rCRS, NC_012920.1, with 16,569 bp.
const HG19_CHRM_LENGTH: usize = 16_571;

/// Warn if the VCF header or the reference genome has the chrM of hg19, whose positions differ
/// from those of the rCRS.
fn warn_on_hg19_chrm(header: &VcfHeader, reference: Option<&provider::Provider>) {
    let source = if header.contigs().iter().any(|(name, contig)| {
        ContigManager::canonicalize(name) == "MT" && contig.length() == Some(HG19_CHRM_LENGTH)
    }) {
        "VCF header"
    } else if reference.and_then(|provider| provider.reference_length("MT"))
        == Some(HG19_CHRM_LENGTH)
    {
        "reference genome"
    } else {
        return;
    };
    tracing::warn!(
        "chrM has {HG19_CHRM_LENGTH} bp in the {source}, so it is NC_001807.4 of hg19, but the \
         transcripts are on the rCRS, NC_012920.1: chrM variants are misannotated"
    );
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

/// Whether `record` can be annotated. Symbolic alleles, breakends and `*` have no bases, `ALT`
/// `.` has no allele, and position 0 lies in front of the first base of the contig.
fn is_annotatable(record: &VcfRecord) -> bool {
    let alternate_bases = record.alternate_bases().as_ref();
    record.variant_start().is_some()
        && !alternate_bases.is_empty()
        && (alternate_bases.iter()).all(|alt| is_sequence_allele(alt))
}

/// Whether `record` is a variant without bases, which is not evaluated together with others,
/// unlike `*`, `<*>` and `<NON_REF>`, which are no variants of their own.
fn is_event(record: &VcfRecord) -> bool {
    !is_annotatable(record)
        && !(record.alternate_bases().as_ref().iter())
            .all(|alt| matches!(alt.as_str(), "*" | "<*>" | "<NON_REF>"))
}

/// The allele of a record, see `classify`: `Err` for one that is not supported.
type Classified = Result<Option<Variant>, SeqvarsError>;

/// The allele of `record`, classified once and passed on with it, see `Variant::from_info`.
/// `file_format` decides what `CN` means, see `copy_number`. An `SVCLAIM` other than `D`, `J`
/// and `DJ` counts as none, and at position 0 only a breakend is supported.
fn classify(record: &VcfRecord, file_format: FileFormat) -> Classified {
    let Var {
        chrom,
        pos,
        reference,
        alternative,
    } = from_vcf_allele(record, 0);
    let var = VcfVariant {
        chromosome: chrom,
        position: pos,
        reference,
        alternative,
    };
    let variant = Variant::from_info(var, &sv_info(record, file_format))?;
    // Position 0 lies in front of the first base of the contig, see `is_annotatable`.
    if pos == 0
        && (variant.as_ref())
            .is_some_and(|variant| !matches!(variant.kind, VariantKind::Breakend(_)))
    {
        return Err(SeqvarsError::InvalidVariant(
            "only a breakend can lie at position 0".into(),
        ));
    }
    Ok(variant)
}

/// The INFO fields of `record` that decide its allele, see `classify`. `file_format` decides
/// what `CN` means, see `copy_number`.
fn sv_info(record: &VcfRecord, file_format: FileFormat) -> SvInfo {
    use noodles::vcf::variant::record::info::field::key::{
        END_CONFIDENCE_INTERVALS, END_POSITION, IS_IMPRECISE, POSITION_CONFIDENCE_INTERVALS,
        SV_TYPE,
    };

    let info = record.info();
    SvInfo {
        end: compound::integer(record, END_POSITION),
        sv_length: compound::sv_length(record),
        sv_type: compound::string(record, SV_TYPE).and_then(|sv_type| sv_type.parse().ok()),
        // An `SVCLAIM` other than those of VCF 4.5 counts as none.
        claim: compound::string(record, "SVCLAIM").and_then(|claim| claim.parse().ok()),
        copy_number: copy_number(record, file_format),
        cipos: confidence_interval(record, POSITION_CONFIDENCE_INTERVALS),
        ciend: confidence_interval(record, END_CONFIDENCE_INTERVALS),
        imprecise: info.get(IS_IMPRECISE).is_some(),
    }
}

/// The copy number of the allele of `record` from `CN`, rounded. Since VCF 4.4, `CN` is that of
/// the allele, before it was the total of a diploid sample, so 0 and 1 copies give 0, 2 give 1,
/// and each more one more. Files without a header line for `CN` give a string below VCF 4.3.
fn copy_number(record: &VcfRecord, file_format: FileFormat) -> Option<i32> {
    use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
    use noodles::vcf::variant::record_buf::info::field::value::Array;

    let copies = match record.info().get("CN").flatten()? {
        InfoValue::Integer(copies) => f64::from(*copies),
        InfoValue::Float(copies) => f64::from(*copies),
        InfoValue::String(copies) => copies.parse().ok()?,
        InfoValue::Array(Array::Integer(copies)) => match copies[..] {
            [Some(copies)] => f64::from(copies),
            _ => return None,
        },
        InfoValue::Array(Array::Float(copies)) => match copies[..] {
            [Some(copies)] => f64::from(copies),
            _ => return None,
        },
        _ => return None,
    };
    let copies = round_copy_number(copies).ok()?;
    Some(if file_format >= FileFormat::new(4, 4) {
        copies
    } else {
        (copies - 1).max(0)
    })
}

/// The confidence interval `key` of `record`, e.g. `CIPOS`, as offsets `(min, max)`: the outer
/// bounds of its pairs, each in order, as a split record keeps the pairs of all its alleles.
/// Files without a header line for `key` give a string below VCF 4.3.
fn confidence_interval(record: &VcfRecord, key: &str) -> Option<ConfidenceInterval> {
    use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
    use noodles::vcf::variant::record_buf::info::field::value::Array;

    let values: Vec<Option<i32>> = match record.info().get(key).flatten()? {
        InfoValue::Array(Array::Integer(values)) => values.clone(),
        InfoValue::String(values) => values.split(',').map(|value| value.parse().ok()).collect(),
        _ => return None,
    };
    let pairs = values.chunks(2).filter_map(|pair| match pair {
        [Some(a), Some(b)] => Some((*a.min(b), *a.max(b))),
        _ => None,
    });
    pairs.reduce(|(min_a, max_a), (min_b, max_b)| (min_a.min(min_b), max_a.max(max_b)))
}

/// The key of the annotation of `var` from `record` in `VariantProcessor`. Two events at one
/// position with one allele can differ in their end.
fn annotation_key(var: &VcfVariant, record: &VcfRecord) -> String {
    if is_event(record) {
        format!("{var}:{}", compound::event_end(record, var.position))
    } else {
        var.to_string()
    }
}

/// Create for all alternate alleles from the given VCF record. A record without alternate
/// allele gives `.`, as `ALT` has it.
pub fn from_vcf_allele(value: &noodles::vcf::variant::RecordBuf, allele_no: usize) -> Var {
    let chrom = value.reference_sequence_name().to_string();
    // Position 0 lies in front of the first base of the contig, e.g., at a telomeric breakend.
    let pos: usize = value.variant_start().map_or(0, |pos| pos.get());
    let pos = i32::try_from(pos).unwrap();
    let reference = value.reference_bases().to_string();
    Var {
        chrom,
        pos,
        reference,
        alternative: (value.alternate_bases().as_ref().get(allele_no))
            .map_or_else(|| ".".to_string(), String::clone),
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantAnnotation {
    pub consequences: Vec<AnnField>,
    /// The gene fusions of a structural variant or breakend, see `Fusion::label`.
    pub fusions: Vec<String>,
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
    /// The annotations of the variants in `buffer` by `annotation_key`, in input order, as
    /// records at one position can have the same key.
    annotations: FxHashMap<String, std::collections::VecDeque<VariantAnnotation>>,
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
    /// `variant` is its allele, see `classify`.
    pub fn process_annotated_record(
        &mut self,
        record: VcfRecord,
        annotation: VariantAnnotation,
        variant: &Classified,
    ) -> anyhow::Result<Vec<AnnotatedVariant>> {
        let mut ready_records = Vec::new();
        let vcf_var = from_vcf_allele(&record, 0);

        if self.buffer.should_flush(&vcf_var.chrom, vcf_var.pos) {
            ready_records.extend(self.process_buffer_flush()?);
            self.buffer.forget_events(&vcf_var.chrom, vcf_var.pos);
        }
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

        // An event reaches its transcripts by its length, as it has no bases. A breakend reaches
        // up to its mate behind it on the same contig, as the junction changes the bases in
        // between, whether or not its group derives an event. A copy number of the reference is
        // no event.
        let is_event = is_event(&record) && !matches!(variant, Ok(None));
        let event_end = is_event.then(|| {
            let end = compound::event_end(&record, csq_var.position);
            match variant {
                Ok(Some(Variant {
                    kind: VariantKind::Breakend(breakend),
                    ..
                })) => (breakend.mate.as_ref())
                    .filter(|mate| mate.chromosome == csq_var.chromosome)
                    .map_or(end, |mate| end.max(mate.position)),
                _ => end,
            }
        });
        // An event is never evaluated together with the variants on its haplotype, which all
        // its rows say, see `VariantBuffer::push_event`.
        let mut annotation = annotation;
        if event_end.is_some() && self.buffer.is_phased(&record, sample_idx, alt_allele_idx) {
            for ann in &mut annotation.consequences {
                (ann.messages.get_or_insert_with(Vec::new))
                    .push(Message::WarningCompoundIncomplete);
            }
        }
        (self.annotations)
            .entry(annotation_key(&csq_var, &record))
            .or_default()
            .push_back(annotation);

        let len = event_end.map_or(csq_var.reference.len(), |end| {
            (end - csq_var.position + 1) as usize
        });
        let (tx_accessions, min_start, max_end) = if let Some(predictor) =
            self.annotator.consequence_predictor()
            && (event_end.is_some() || is_annotatable(&record))
        {
            consequence::get_transcript_boundaries(
                predictor,
                &csq_var.chromosome,
                csq_var.position,
                len,
            )
        } else {
            (std::collections::HashSet::new(), i32::MAX, -1)
        };

        if let Some(end) = event_end {
            self.buffer.push_event(
                csq_var,
                record,
                tx_accessions,
                max_end,
                end,
                sample_idx,
                alt_allele_idx,
            );
        } else {
            self.buffer.push(
                csq_var,
                record,
                tx_accessions,
                min_start,
                max_end,
                sample_idx,
                alt_allele_idx,
            );
        }

        Ok(ready_records)
    }

    pub fn flush(&mut self) -> anyhow::Result<Vec<AnnotatedVariant>> {
        self.process_buffer_flush()
    }

    fn process_buffer_flush(&mut self) -> anyhow::Result<Vec<AnnotatedVariant>> {
        let mut incomplete = self.buffer.transcripts_with_events();
        let (variants, compound_groups) = self.buffer.flush();
        let predictor = self.annotator.consequence_predictor();

        let mut annotations = std::mem::take(&mut self.annotations);
        let mut result_annotations: Vec<VariantAnnotation> = variants
            .iter()
            .map(|b| {
                (annotations.get_mut(&annotation_key(&b.vcf_var, &b.record)))
                    .and_then(std::collections::VecDeque::pop_front)
                    .unwrap_or_else(|| VariantAnnotation {
                        consequences: vec![],
                        fusions: vec![],
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
                        // The group forms on all transcripts that its variants share.
                        for &i in &group_indices {
                            let shared = (variants[i].tx_accessions.iter()).filter(|tx| {
                                (group_indices.iter())
                                    .all(|&j| variants[j].tx_accessions.contains(*tx))
                            });
                            incomplete[i].extend(shared.cloned());
                        }
                    }
                    // Other errors, e.g. of a plugin, stop the run, as for single variants.
                    Err(e) => {
                        return Err(Error::from(e)
                            .context(format!("Failed to annotate compound variant {}", group())));
                    }
                }
            }
        }

        for (annotation, transcripts) in result_annotations.iter_mut().zip(&incomplete) {
            for ann in &mut annotation.consequences {
                if transcripts.contains(&ann.feature_id) {
                    (ann.messages.get_or_insert_with(Vec::new))
                        .push(Message::WarningCompoundIncomplete);
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
    use super::build_header;
    use super::classify;
    use super::consequence::{
        Breakend, ConfidenceInterval, ConsequenceAnnotator, StructuralVariant, SvKind, Variant,
        VariantKind, VcfVariant, load_tx_db,
    };
    use super::{
        AnnotatedVariant, Annotator, AnnotatorEnum, Args, OutputFormat, VariantProcessor,
        VcfRecord, run,
    };
    use crate::annotate::cli::{CompoundSettings, PhasingStrategy};
    use crate::annotate::cli::{PickUnit, PredictorSettings};
    use crate::annotate::cli::{Sources, TranscriptSettings};
    use crate::annotate::seqvars::consequence::plugin::{BuiltinPlugin, PluginSpec};
    use crate::annotate::seqvars::consequence::terms::{
        ANN_COMPOUND_IDS, ANN_COMPOUND_VARIANTS, AnnField,
    };
    use crate::common::noodles::{NoodlesVariantReader, open_variant_reader};
    use clap::Parser;
    use clap_verbosity_flag::Verbosity;
    use futures::TryStreamExt;
    use noodles::vcf::Header as VcfHeader;
    use noodles::vcf::header::FileFormat;
    use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
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
            ann_key: "ANN".into(),
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
                ..Default::default()
            },
        };

        run(&args_common, &args).await?;

        let actual = std::fs::read_to_string(args.output)?;
        // remove vcf header lines starting with ##mehari, and normalize the version of mehari
        // that the `ANN` header line carries, as it changes with every release
        let version = format!("Version=\"{}\"", env!("CARGO_PKG_VERSION"));
        let actual = actual
            .lines()
            .filter(|line| !line.starts_with("##mehari"))
            .map(|line| line.replace(&version, "Version=\"[VERSION]\""))
            .collect::<Vec<_>>()
            .join("\n");
        insta::assert_snapshot!(actual);

        Ok(())
    }

    /// Only the fields that the run writes are taken, not those of sources it does not use.
    #[rstest::rstest]
    #[case::frequencies("gnomad_exomes_an", true, false, false)]
    #[case::clinvar("clinvar_vcv", false, true, false)]
    #[case::source_not_used("clinvar_vcv", false, false, true)]
    #[case::other_source_used("clinvar_vcv", true, false, true)]
    #[case::default_key("ANN", true, true, true)]
    fn ann_key_must_differ_from_other_fields_of_the_run(
        #[case] ann_key: &str,
        #[case] with_frequencies: bool,
        #[case] with_clinvar: bool,
        #[case] valid: bool,
    ) {
        let header = build_header(
            &VcfHeader::default(),
            Some(ann_key),
            with_frequencies,
            with_clinvar,
            &[],
            &[],
            &[],
        );
        assert_eq!(header.is_ok(), valid, "{ann_key}");
    }

    #[tokio::test]
    async fn ann_key_of_another_field_of_the_run_stops_it() {
        let sources = Sources {
            frequencies: Some(vec!["tests/data/annotate/db/grch37/seqvars/freqs".into()]),
            ..transcripts()
        };
        let error =
            annotate_input_with_field(sources, "ANN", "gnomad_exomes_an", OutputFormat::Vcf)
                .await
                .unwrap_err();
        assert!(
            error.to_string().contains("--ann-key gnomad_exomes_an"),
            "{error}"
        );
    }

    /// The `ANN` header names mehari as its source, and a reader still finds the field names.
    #[test]
    fn ann_header_has_source_and_version() -> Result<(), anyhow::Error> {
        let header = build_header(
            &VcfHeader::default(),
            Some("ANN"),
            false,
            false,
            &[],
            &[],
            &[],
        )?;

        let mut buf = Vec::new();
        noodles::vcf::io::Writer::new(&mut buf).write_header(&header)?;
        let text = String::from_utf8(buf)?;
        let line = text
            .lines()
            .find(|line| line.starts_with("##INFO=<ID=ANN,"))
            .unwrap();
        assert!(line.contains(",Source=\"mehari\""), "{line}");
        assert!(
            line.ends_with(&format!(",Version=\"{}\">", env!("CARGO_PKG_VERSION"))),
            "{line}"
        );

        let read_back: VcfHeader = text.parse()?;
        assert_eq!(
            AnnField::ann_field_names_from_header(&read_back, "ANN"),
            Some(AnnField::ann_field_names(&[]))
        );
        let info = read_back.infos().get("ANN").unwrap();
        assert_eq!(
            info.other_fields().get("Source").map(|v| v.as_str()),
            Some("mehari")
        );

        Ok(())
    }

    /// Annotates a VCF whose field `input_key` comes from another tool, and returns the output
    /// lines.
    ///
    /// The records have no transcript nearby, a variant in BRCA1, and a deletion allele.
    async fn annotate_input_with_field(
        sources: Sources,
        input_key: &str,
        ann_key: &str,
        output_format: OutputFormat,
    ) -> Result<Vec<String>, anyhow::Error> {
        let temp = TempDir::default();
        let path_in = temp.join("input.vcf");
        let path_out = temp.join("output.vcf");
        std::fs::write(
            &path_in,
            format!(
                "##fileformat=VCFv4.2\n\
                 ##INFO=<ID={input_key},Number=.,Type=String,Description=\"Functional annotations: \
                 'Allele | Annotation'\">\n\
                 ##contig=<ID=17,length=81195210>\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 17\t1000\t.\tA\tT\t.\t.\t{input_key}=T|from_input\n\
                 17\t41215920\t.\tG\tT\t.\t.\t{input_key}=T|from_input\n\
                 17\t41215921\t.\tG\t*\t.\t.\t{input_key}=*|from_input\n"
            ),
        )?;

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            threads: 1,
            reference: None,
            in_memory_reference: true,
            assembly: Some("grch37".into()),
            predictor_settings: PredictorSettings::default(),
            input: path_in.into_os_string().into_string().unwrap(),
            output: path_out.into_os_string().into_string().unwrap(),
            output_format,
            max_var_count: None,
            ann_key: ann_key.into(),
            sources,
        };
        run(&args_common, &args).await?;

        Ok(std::fs::read_to_string(&args.output)?
            .lines()
            .map(String::from)
            .collect())
    }

    /// The message on rows of variants that were not evaluated with all their phased partners.
    const INCOMPLETE: &str = "WARNING_COMPOUND_INCOMPLETE";

    /// Two SNVs in BRCA1 that form a compound group, and records without bases.
    const RECORDS_WITHOUT_BASES: &str = "\
        17\t0\ttelomere\tN\t.[13:32900000[\t.\t.\t.\tGT\t1/1\n\
        17\t0\ttelomere_ins\tN\tNA\t.\t.\t.\tGT\t1/1\n\
        17\t41197701\tsnv1\tG\tC\t.\t.\t.\tGT\t1/1\n\
        17\t41197710\tdel\tN\t<DEL>\t.\t.\t.\tGT\t1/1\n\
        17\t41215920\tsnv2\tG\tT\t.\t.\t.\tGT\t1/1\n\
        17\t41215930\tbnd\tG\tG]13:32900000]\t.\t.\t.\tGT\t1/1\n\
        17\t41215940\tstar\tG\t*\t.\t.\t.\tGT\t1/1\n\
        17\t41215950\tany\tG\t<*>\t.\t.\t.\tGT\t1/1\n\
        17\t41215960\tnon_ref\tG\t<NON_REF>\t.\t.\t.\tGT\t1/1\n\
        17\t41215970\tno_alt\tG\t.\t.\t.\t.\tGT\t0/0\n\
        17\t41215980\tcn2\tG\t<CNV>\t.\t.\tEND=41216000;CN=2\tGT\t0/1\n";

    /// Annotates `records`, and returns the ID and the `ANN` rows of each record.
    async fn annotate_records(
        records: &str,
        compound: bool,
    ) -> Result<Vec<(String, Vec<String>)>, anyhow::Error> {
        let mut predictor_settings = PredictorSettings::default();
        predictor_settings
            .compound_settings
            .enable_compound_variants = compound;
        predictor_settings.compound_settings.phasing_strategy = PhasingStrategy::Strict;
        annotate_records_with(records, predictor_settings).await
    }

    /// `annotate_records` with `predictor_settings`.
    async fn annotate_records_with(
        records: &str,
        predictor_settings: PredictorSettings,
    ) -> Result<Vec<(String, Vec<String>)>, anyhow::Error> {
        let temp = TempDir::default();
        let path_in = temp.join("input.vcf");
        let path_out = temp.join("output.vcf");
        std::fs::write(
            &path_in,
            format!(
                "##fileformat=VCFv4.2\n\
                 ##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
                 ##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length\">\n\
                 ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise\">\n\
                 ##INFO=<ID=CIPOS,Number=.,Type=Integer,Description=\"Interval\">\n\
                 ##INFO=<ID=MATEID,Number=A,Type=String,Description=\"ID of mate breakend\">\n\
                 ##INFO=<ID=EVENT,Number=A,Type=String,Description=\"ID of associated event\">\n\
                 ##INFO=<ID=EVENTTYPE,Number=A,Type=String,Description=\"Type of associated event\">\n\
                 ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                 ##contig=<ID=13,length=115169878>\n\
                 ##contig=<ID=17,length=81195210>\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n\
                 {records}"
            ),
        )?;

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
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
            ann_key: "ANN".into(),
            sources: transcripts(),
        };
        run(&args_common, &args).await?;

        Ok(std::fs::read_to_string(&args.output)?
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let fields = line.split('\t').collect::<Vec<_>>();
                let rows = (fields[7].split(';'))
                    .find_map(|info| info.strip_prefix("ANN="))
                    .map_or_else(Vec::new, |anns| anns.split(',').map(String::from).collect());
                (fields[2].to_string(), rows)
            })
            .collect())
    }

    /// Records that cannot be annotated stay in the output without annotations, also at position
    /// 0, where a breakend touches no transcript. Mehari counts them in one warning, but not that
    /// breakend, nor `*`, `<*>`, `<NON_REF>` and `ALT` `.`, which are no variants of their own,
    /// nor a `<CNV>` with the total copy number 2 of VCF 4.2.
    #[tokio::test]
    #[tracing_test::traced_test]
    async fn records_without_bases_pass_through() -> Result<(), anyhow::Error> {
        let records = annotate_records(RECORDS_WITHOUT_BASES, false).await?;

        let annotated = records
            .iter()
            .map(|(id, rows)| (id.as_str(), !rows.is_empty()))
            .collect::<Vec<_>>();
        assert_eq!(
            annotated,
            [
                ("telomere", false),
                ("telomere_ins", false),
                ("snv1", true),
                ("del", false),
                ("snv2", true),
                ("bnd", true),
                ("star", false),
                ("any", false),
                ("non_ref", false),
                ("no_alt", false),
                ("cn2", false),
            ]
        );
        assert!(logs_contain(
            "2 record(s) stay without annotations: malformed breakends"
        ));
        // Without compound mode, nothing is evaluated together.
        let rows = records.iter().flat_map(|(_, rows)| rows);
        assert!(rows.into_iter().all(|row| !row.contains(INCOMPLETE)));
        Ok(())
    }

    /// A record with more than one alternate allele aborts the run, with the command that splits
    /// it.
    #[tokio::test]
    #[tracing_test::traced_test]
    async fn multi_allelic_record_names_split_command() {
        let records = "17\t41197701\tsnv\tG\tC,T\t.\t.\t.\tGT\t1/2\n";
        assert!(annotate_records(records, false).await.is_err());
        assert!(logs_contain(
            "split multi-allelic records with `bcftools norm -m-any`, without `-f`, which rejects \
             `N` as `REF`"
        ));
    }

    /// An `SVCLAIM` other than `D`, `J` and `DJ` counts as none, so this `<DUP>` is no tandem
    /// dup. Mehari counts such records in one warning.
    #[tokio::test]
    #[tracing_test::traced_test]
    async fn unknown_claim_counts_as_none() -> Result<(), anyhow::Error> {
        let records = "17\t41197710\tdup\tN\t<DUP>\t.\t.\tSVLEN=50;SVCLAIM=JD\tGT\t0/1\n";
        let records = annotate_records(records, false).await?;
        assert!(
            (records[0].1.iter()).any(|row| row.starts_with("<DUP>|copy_number_increase|")),
            "{records:?}"
        );
        assert!(logs_contain(
            "1 record(s) have an SVCLAIM other than D, J and DJ, which counts as none"
        ));
        Ok(())
    }

    /// An event is not evaluated with the variants phased with it, but all their rows say so, as
    /// do the rows of the breakend. The other variants are still evaluated together.
    #[tokio::test]
    async fn phased_event_marks_rows_of_partners() -> Result<(), anyhow::Error> {
        for (id, rows) in annotate_records(RECORDS_WITHOUT_BASES, true).await? {
            if id.starts_with("snv") {
                assert!(rows.iter().any(|row| row.contains("|comp_1|")), "{rows:?}");
                assert!(rows.iter().all(|row| row.contains(INCOMPLETE)), "{rows:?}");
            } else if id == "bnd" {
                assert!(!rows.is_empty());
                assert!(rows.iter().all(|row| row.contains(INCOMPLETE)), "{rows:?}");
            } else {
                assert!(rows.is_empty(), "{id}: {rows:?}");
            }
        }
        Ok(())
    }

    /// A record at 1000 with the allele `alternative` and the INFO fields `info`.
    fn sv_record(alternative: &str, info: Vec<(&str, InfoValue)>) -> VcfRecord {
        use noodles::core::Position;
        use noodles::vcf::variant::record_buf::AlternateBases;

        VcfRecord::builder()
            .set_variant_start(Position::try_from(1000).unwrap())
            .set_alternate_bases(AlternateBases::from(vec![alternative.to_string()]))
            .set_info(
                (info.into_iter())
                    .map(|(key, value)| (key.to_string(), Some(value)))
                    .collect(),
            )
            .build()
    }

    type Intervals = (Option<ConfidenceInterval>, Option<ConfidenceInterval>);

    /// The structural variant of `record`, see `classify`.
    fn structural_variant(
        record: &VcfRecord,
        file_format: FileFormat,
    ) -> Option<StructuralVariant> {
        match classify(record, file_format) {
            Ok(Some(Variant {
                kind: VariantKind::Structural(sv),
                ..
            })) => Some(sv),
            _ => None,
        }
    }

    /// The breakend of `record`, see `classify`.
    fn breakend(record: &VcfRecord) -> Option<Breakend> {
        match classify(record, FileFormat::default()) {
            Ok(Some(Variant {
                kind: VariantKind::Breakend(breakend),
                ..
            })) => Some(breakend),
            _ => None,
        }
    }

    /// The structural variant of a record, see `structural_variant`. `SVLEN` as `Number=A` and
    /// `SVCLAIM` as `Number=A` come as arrays of one value.
    #[rstest::rstest]
    #[case::deletion("<DEL>", vec![("SVLEN", InfoValue::Integer(-50))], Some((SvKind::Deletion, 1050, false, (None, None))))]
    #[case::deletion_by_end("<DEL:ME:ALU>", vec![("END", InfoValue::Integer(1500))], Some((SvKind::Deletion, 1500, false, (None, None))))]
    #[case::deletion_without_length("<DEL>", vec![], None)]
    #[case::abundance_claim("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("SVCLAIM", InfoValue::Array(Array::String(vec![Some("D".into())])))], Some((SvKind::Deletion, 1050, true, (None, None))))]
    #[case::imprecise_flag("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("IMPRECISE", InfoValue::Flag)], Some((SvKind::Deletion, 1050, true, (None, None))))]
    #[case::zero_intervals("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("CIPOS", InfoValue::Array(Array::Integer(vec![Some(0), Some(0)])))], Some((SvKind::Deletion, 1050, false, (Some((0, 0)), Some((0, 0))))))]
    #[case::intervals("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("CIPOS", InfoValue::Array(Array::Integer(vec![Some(-10), Some(10)]))), ("CIEND", InfoValue::Array(Array::Integer(vec![Some(-20), Some(20)])))], Some((SvKind::Deletion, 1050, false, (Some((-10, 10)), Some((-20, 20))))))]
    #[case::imprecise_zero_intervals("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("IMPRECISE", InfoValue::Flag), ("CIPOS", InfoValue::Array(Array::Integer(vec![Some(0), Some(0)])))], Some((SvKind::Deletion, 1050, false, (Some((0, 0)), Some((0, 0))))))]
    #[case::end_interval_only("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("CIEND", InfoValue::Array(Array::Integer(vec![Some(-20), Some(20)])))], Some((SvKind::Deletion, 1050, false, (None, Some((-20, 20))))))]
    #[case::split_record_outer_intervals("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("CIPOS", InfoValue::Array(Array::Integer(vec![Some(-5), Some(5), Some(-20), Some(20)])))], Some((SvKind::Deletion, 1050, false, (Some((-20, 20)), Some((-20, 20))))))]
    #[case::imprecise_intervals("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("IMPRECISE", InfoValue::Flag), ("CIPOS", InfoValue::Array(Array::Integer(vec![Some(-10), Some(10)]))), ("CIEND", InfoValue::Array(Array::Integer(vec![Some(-20), Some(20)])))], Some((SvKind::Deletion, 1050, true, (Some((-10, 10)), Some((-20, 20))))))]
    #[case::imprecise_end_interval_only("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("IMPRECISE", InfoValue::Flag), ("CIEND", InfoValue::Array(Array::Integer(vec![Some(-20), Some(20)])))], Some((SvKind::Deletion, 1050, true, (None, Some((-20, 20))))))]
    #[case::imprecise_zero_end_interval_only("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("IMPRECISE", InfoValue::Flag), ("CIEND", InfoValue::Array(Array::Integer(vec![Some(0), Some(0)])))], Some((SvKind::Deletion, 1050, true, (None, Some((0, 0))))))]
    #[case::imprecise_reversed_interval("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("IMPRECISE", InfoValue::Flag), ("CIPOS", InfoValue::Array(Array::Integer(vec![Some(10), Some(-10)])))], Some((SvKind::Deletion, 1050, true, (Some((-10, 10)), Some((-10, 10))))))]
    #[case::interval_strings_without_header("<DEL>", vec![("SVLEN", InfoValue::Integer(-50)), ("IMPRECISE", InfoValue::Flag), ("CIPOS", InfoValue::String("-10,10".into())), ("CIEND", InfoValue::String("-20,20".into()))], Some((SvKind::Deletion, 1050, true, (Some((-10, 10)), Some((-20, 20))))))]
    #[case::dup("<DUP>", vec![("SVLEN", InfoValue::Integer(50))], Some((SvKind::Duplication { tandem: false }, 1050, false, (None, None))))]
    #[case::dup_tandem("<DUP:TANDEM>", vec![("SVLEN", InfoValue::Integer(50))], Some((SvKind::Duplication { tandem: true }, 1050, false, (None, None))))]
    #[case::dup_junction_claim("<DUP>", vec![("SVLEN", InfoValue::Integer(50)), ("SVCLAIM", InfoValue::String("DJ".into()))], Some((SvKind::Duplication { tandem: true }, 1050, false, (None, None))))]
    #[case::dup_tandem_abundance_claim("<DUP:TANDEM>", vec![("SVLEN", InfoValue::Integer(50)), ("SVCLAIM", InfoValue::String("D".into()))], Some((SvKind::Duplication { tandem: false }, 1050, true, (None, None))))]
    #[case::unknown_claim("<DUP>", vec![("SVLEN", InfoValue::Integer(50)), ("SVCLAIM", InfoValue::String("JD".into()))], Some((SvKind::Duplication { tandem: false }, 1050, false, (None, None))))]
    #[case::inversion("<INV>", vec![("END", InfoValue::Integer(1500))], Some((SvKind::Inversion, 1500, false, (None, None))))]
    #[case::inversion_length_difference("<INV>", vec![("SVLEN", InfoValue::Integer(0)), ("END", InfoValue::Integer(1500))], Some((SvKind::Inversion, 1500, false, (None, None))))]
    #[case::strings_without_header("<DEL>", vec![("SVLEN", InfoValue::String("-50".into())), ("END", InfoValue::String("1050".into()))], Some((SvKind::Deletion, 1050, false, (None, None))))]
    #[case::end_string_without_header("<DEL>", vec![("END", InfoValue::String("1050".into()))], Some((SvKind::Deletion, 1050, false, (None, None))))]
    #[case::insertion("<INS>", vec![("SVLEN", InfoValue::Integer(300))], Some((SvKind::Insertion { length: Some(300) }, 1000, false, (None, None))))]
    #[case::insertion_unknown_length("<INS:ME:ALU>", vec![], Some((SvKind::Insertion { length: None }, 1000, false, (None, None))))]
    #[case::insertion_length_0("<INS>", vec![("SVLEN", InfoValue::Integer(0))], Some((SvKind::Insertion { length: None }, 1000, false, (None, None))))]
    #[case::cnv("<CNV>", vec![("SVLEN", InfoValue::Integer(50))], Some((SvKind::CopyNumberChange, 1050, true, (None, None))))]
    #[case::cnv_lost("<CNV>", vec![("SVLEN", InfoValue::Integer(50)), ("CN", InfoValue::Integer(0))], Some((SvKind::Deletion, 1050, true, (None, None))))]
    #[case::cnv_gained("<CNV>", vec![("SVLEN", InfoValue::Integer(50)), ("CN", InfoValue::Array(Array::Integer(vec![Some(3)])))], Some((SvKind::Duplication { tandem: false }, 1050, true, (None, None))))]
    #[case::cn0("<CN0>", vec![("END", InfoValue::Integer(1050))], Some((SvKind::Deletion, 1050, true, (None, None))))]
    #[case::cn1("<CN1>", vec![("END", InfoValue::Integer(1050))], None)]
    #[case::cn2("<CN2>", vec![("END", InfoValue::Integer(1050))], Some((SvKind::Duplication { tandem: false }, 1050, true, (None, None))))]
    #[case::cn_with_sign("<CN+2>", vec![("END", InfoValue::Integer(1050))], None)]
    #[case::cn_negative("<CN-1>", vec![("END", InfoValue::Integer(1050))], None)]
    #[case::tandem_repeat_copy_number("<CNV:TR>", vec![("END", InfoValue::Integer(1050)), ("CN", InfoValue::Float(2.5))], None)]
    #[case::complex("<CPX>", vec![("END", InfoValue::Integer(1050))], None)]
    #[case::custom_by_sv_type("<MYDEL>", vec![("END", InfoValue::Integer(1050)), ("SVTYPE", InfoValue::String("DEL".into()))], Some((SvKind::Deletion, 1050, false, (None, None))))]
    #[case::custom_by_other_sv_type("<MYDEL>", vec![("END", InfoValue::Integer(1050)), ("SVTYPE", InfoValue::String("CN0".into()))], None)]
    #[case::no_variant_with_sv_type("<*>", vec![("END", InfoValue::Integer(1050)), ("SVTYPE", InfoValue::String("DEL".into()))], None)]
    #[case::breakend("G]13:100]", vec![("END", InfoValue::Integer(1050))], None)]
    fn structural_variant_of_record(
        #[case] alternative: &str,
        #[case] info: Vec<(&str, InfoValue)>,
        #[case] expected: Option<(SvKind, i32, bool, Intervals)>,
    ) {
        let sv = structural_variant(&sv_record(alternative, info), FileFormat::default());
        let fields = sv.map(|sv| (sv.kind, sv.end, sv.imprecise, sv.intervals));
        assert_eq!(fields, expected);
    }

    /// `CN` of a `<CNV>` is the copy number of the allele since VCF 4.4, as a float, and the
    /// total of a diploid sample before. With the copy number of the reference, there is no
    /// variant, and the record is supported.
    #[rstest::rstest]
    #[case::lost((4, 5), InfoValue::Array(Array::Float(vec![Some(0.0)])), Some(SvKind::Deletion))]
    #[case::kept((4, 5), InfoValue::Array(Array::Float(vec![Some(1.0)])), None)]
    #[case::gained((4, 4), InfoValue::Float(1.6), Some(SvKind::Duplication { tandem: false }))]
    #[case::total_lost_one((4, 3), InfoValue::Integer(1), Some(SvKind::Deletion))]
    #[case::total_kept((4, 3), InfoValue::Integer(2), None)]
    #[case::total_gained((4, 3), InfoValue::Integer(3), Some(SvKind::Duplication { tandem: false }))]
    #[case::total_string_without_header((4, 2), InfoValue::String("0".into()), Some(SvKind::Deletion))]
    fn copy_number_by_file_format(
        #[case] (major, minor): (u32, u32),
        #[case] copy_number: InfoValue,
        #[case] expected: Option<SvKind>,
    ) {
        let record = sv_record(
            "<CNV>",
            vec![("END", InfoValue::Integer(1050)), ("CN", copy_number)],
        );
        let file_format = FileFormat::new(major, minor);
        let sv = structural_variant(&record, file_format);
        assert_eq!(sv.map(|sv| sv.kind), expected);
        assert_eq!(
            matches!(classify(&record, file_format), Ok(None)),
            expected.is_none()
        );
    }

    /// The uncertainty of a breakend from the INFO of its record, see `Breakend::imprecise`.
    /// Without `IMPRECISE`, `CIPOS` is the micro-homology of a precise call. A split record keeps
    /// the pairs of all its alleles.
    #[rstest::rstest]
    #[case::precise(vec![], (false, None))]
    #[case::imprecise_flag(vec![("IMPRECISE", InfoValue::Flag)], (true, None))]
    #[case::interval_without_flag(vec![("CIPOS", InfoValue::Array(Array::Integer(vec![Some(-500), Some(500)])))], (false, Some((-500, 500))))]
    #[case::imprecise_interval(vec![("IMPRECISE", InfoValue::Flag), ("CIPOS", InfoValue::Array(Array::Integer(vec![Some(-500), Some(500)])))], (true, Some((-500, 500))))]
    #[case::imprecise_zero_interval(vec![("IMPRECISE", InfoValue::Flag), ("CIPOS", InfoValue::Array(Array::Integer(vec![Some(0), Some(0)])))], (false, Some((0, 0))))]
    #[case::split_record(vec![("IMPRECISE", InfoValue::Flag), ("CIPOS", InfoValue::Array(Array::Integer(vec![Some(-10), Some(10), Some(0), Some(0)])))], (true, Some((-10, 10))))]
    fn breakend_of_record(
        #[case] info: Vec<(&str, InfoValue)>,
        #[case] expected: (bool, Option<(i32, i32)>),
    ) {
        let breakend = breakend(&sv_record("G]13:32900000]", info));
        let fields = breakend.map(|breakend| (breakend.imprecise, breakend.cipos));
        assert_eq!(fields, Some(expected));
    }

    /// A `<DEL>` with its length is annotated as the deletion of its bases. In compound mode, its
    /// rows say that it was not evaluated with the variants on its haplotype, as the rows of
    /// those variants do.
    #[rstest::rstest]
    #[case::phased("0|1", true)]
    #[case::unphased("0/1", false)]
    #[tokio::test]
    async fn symbolic_deletion_is_annotated(
        #[case] genotype: &str,
        #[case] marked: bool,
    ) -> Result<(), anyhow::Error> {
        let records = format!(
            "17\t41197701\tsnv\tG\tC\t.\t.\t.\tGT\t0|1\n\
             17\t41197710\tdel\tN\t<DEL>\t.\t.\tSVLEN=-50\tGT\t{genotype}\n"
        );
        let records = annotate_records(&records, true).await?;
        let (snv, del) = (&records[0].1, &records[1].1);
        assert!(
            del.iter().any(
                |row| row.starts_with("<DEL>|3_prime_UTR_exon_variant|MODIFIER|BRCA1|")
                    && row.contains("|g.41197711_41197760del|")
            ),
            "{del:?}"
        );
        assert!(
            del.iter().all(|row| row.contains(INCOMPLETE) == marked),
            "{del:?}"
        );
        assert!(
            snv.iter().all(|row| row.contains(INCOMPLETE) == marked),
            "{snv:?}"
        );
        Ok(())
    }

    /// Two `<DEL>` of one length at one position keep their own rows in compound mode: only the
    /// phased one was not evaluated with the variants on its haplotype.
    #[tokio::test]
    async fn symbolic_deletions_at_one_position() -> Result<(), anyhow::Error> {
        let records = "17\t41197701\tsnv\tG\tC\t.\t.\t.\tGT\t0|1\n\
                       17\t41197710\tphased\tN\t<DEL>\t.\t.\tSVLEN=-50\tGT\t0|1\n\
                       17\t41197710\tunphased\tN\t<DEL>\t.\t.\tSVLEN=-50\tGT\t0/1\n";
        let records = annotate_records(records, true).await?;
        for (id, rows) in &records[1..] {
            let marked = id == "phased";
            assert!(!rows.is_empty(), "{id}");
            assert!(
                rows.iter().all(|row| row.contains(INCOMPLETE) == marked),
                "{id}: {rows:?}"
            );
        }
        Ok(())
    }

    /// A `<DEL>` with the largest END ends at the end of the contig, also in compound mode.
    #[rstest::rstest]
    #[tokio::test]
    async fn symbolic_deletion_with_largest_end(
        #[values(false, true)] compound: bool,
    ) -> Result<(), anyhow::Error> {
        let records = "17\t41197710\tdel\tN\t<DEL>\t.\t.\tEND=2147483647\tGT\t0|1\n";
        let rows = &annotate_records(records, compound).await?[0].1;
        assert!(
            rows.iter()
                .any(|row| row.contains("|g.41197711_81195210del|")),
            "{rows:?}"
        );
        Ok(())
    }

    /// A breakend is annotated from its allele alone. It lies in BRCA1, on the minus strand, and
    /// keeps the bases in front of it, so the transcripts lose their 5' ends. In compound mode,
    /// its rows say that it was not evaluated with the variants on its haplotype, as the rows of
    /// those variants do.
    #[rstest::rstest]
    #[case::phased("0|1", true)]
    #[case::unphased("0/1", false)]
    #[tokio::test]
    async fn breakend_is_annotated(
        #[case] genotype: &str,
        #[case] marked: bool,
    ) -> Result<(), anyhow::Error> {
        let records = format!(
            "17\t41197701\tsnv\tG\tC\t.\t.\t.\tGT\t0|1\n\
             17\t41215930\tbnd\tG\tG]13:32900000]\t.\t.\t.\tGT\t{genotype}\n"
        );
        let records = annotate_records(&records, true).await?;
        let (snv, bnd) = (&records[0].1, &records[1].1);
        assert!(
            bnd.iter().any(|row| row.starts_with("G]13:32900000]|")
                && row.contains("start_lost")
                && row.contains("feature_truncation")
                && row.contains("|g.41215931_qterdelins[NC_000013.10:g.pter_32900000inv]|")),
            "{bnd:?}"
        );
        assert!(
            bnd.iter().all(|row| row.contains(INCOMPLETE) == marked),
            "{bnd:?}"
        );
        assert!(
            snv.iter().all(|row| row.contains(INCOMPLETE) == marked),
            "{snv:?}"
        );
        Ok(())
    }

    /// Whether `row` is one of a breakend, whose `HGVS.g` joins a piece of the mate contig.
    fn is_breakend_row(row: &str) -> bool {
        row.contains("delins[")
    }

    /// The breakends of a deletion, a tandem duplication, a double minute or an inversion of all
    /// BRCA1 transcripts, which lie between 41196312 and 41277500, give their last record the
    /// row of the event, one for the gene, with the allele of the record, next to its own rows:
    /// those of a breakend downstream or upstream of the transcripts. The other records keep
    /// their own rows. `EVENT` groups them as `MATEID` does, an `EVENTTYPE` of their event keeps
    /// it, and that of a double minute makes it a copy elsewhere.
    #[rstest::rstest]
    #[case::deletion(
        "17\t41195000\tbnd1\tG\tG[17:41280000[\t.\t.\tMATEID=bnd2\tGT\t0/1\n\
         17\t41280000\tbnd2\tT\t]17:41195000]T\t.\t.\tMATEID=bnd1\tGT\t0/1\n",
        "transcript_ablation|HIGH|BRCA1|",
        "|g.41195001_41279999del|"
    )]
    #[case::deletion_by_event(
        "17\t41195000\tbnd1\tG\tG[17:41280000[\t.\t.\tEVENT=del\tGT\t0/1\n\
         17\t41280000\tbnd2\tT\t]17:41195000]T\t.\t.\tEVENT=del\tGT\t0/1\n",
        "transcript_ablation|HIGH|BRCA1|",
        "|g.41195001_41279999del|"
    )]
    #[case::tandem_duplication(
        "17\t41195000\tbnd1\tG\t]17:41280000]G\t.\t.\tMATEID=bnd2\tGT\t0/1\n\
         17\t41280000\tbnd2\tT\tT[17:41195000[\t.\t.\tMATEID=bnd1\tGT\t0/1\n",
        "|BRCA1|",
        "|g.41195000_41280000dup|"
    )]
    #[case::tandem_duplication_by_event_type(
        "17\t41195000\tbnd1\tG\t]17:41280000]G\t.\t.\tMATEID=bnd2;EVENTTYPE=DUP:TANDEM\tGT\t0/1\n\
         17\t41280000\tbnd2\tT\tT[17:41195000[\t.\t.\tMATEID=bnd1;EVENTTYPE=DUP:TANDEM\tGT\t0/1\n",
        "|BRCA1|",
        "|g.41195000_41280000dup|"
    )]
    #[case::double_minute(
        "17\t41195000\tbnd1\tG\t]17:41280000]G\t.\t.\tMATEID=bnd2;EVENTTYPE=DOUBLEMINUTE\tGT\t0/1\n\
         17\t41280000\tbnd2\tT\tT[17:41195000[\t.\t.\tMATEID=bnd1;EVENTTYPE=DOUBLEMINUTE\tGT\t0/1\n",
        "transcript_amplification|HIGH|BRCA1|",
        "|g.?_?ins[NC_000017.10:g.41195000_41280000]|"
    )]
    #[case::inversion(
        "17\t41195000\tbnd1\tG\tG]17:41280000]\t.\t.\tMATEID=bnd3;EVENT=inv\tGT\t0/1\n\
         17\t41195001\tbnd2\tA\t[17:41280001[A\t.\t.\tMATEID=bnd4;EVENT=inv\tGT\t0/1\n\
         17\t41280000\tbnd3\tT\tT]17:41195000]\t.\t.\tMATEID=bnd1;EVENT=inv\tGT\t0/1\n\
         17\t41280001\tbnd4\tC\t[17:41195001[C\t.\t.\tMATEID=bnd2;EVENT=inv\tGT\t0/1\n",
        "|BRCA1|",
        "|g.41195001_41280000inv|"
    )]
    #[case::inversion_by_event(
        "17\t41195000\tbnd1\tG\tG]17:41280000]\t.\t.\tEVENT=inv\tGT\t0/1\n\
         17\t41195001\tbnd2\tA\t[17:41280001[A\t.\t.\tEVENT=inv\tGT\t0/1\n\
         17\t41280000\tbnd3\tT\tT]17:41195000]\t.\t.\tEVENT=inv\tGT\t0/1\n\
         17\t41280001\tbnd4\tC\t[17:41195001[C\t.\t.\tEVENT=inv\tGT\t0/1\n",
        "|BRCA1|",
        "|g.41195001_41280000inv|"
    )]
    #[tokio::test]
    async fn breakend_group_gives_event_rows(
        #[case] records: &str,
        #[case] term: &str,
        #[case] hgvs_g: &str,
    ) -> Result<(), anyhow::Error> {
        let alts = records.lines().map(|line| line.split('\t').nth(4).unwrap());
        let annotated = annotate_records(records, false).await?;
        assert_eq!(annotated.len(), alts.clone().count());
        for (i, ((id, rows), alt)) in annotated.iter().zip(alts).enumerate() {
            assert!(
                rows.iter().any(|row| is_breakend_row(row)),
                "{id}: {rows:?}"
            );
            let event = (rows.iter()).filter(|row| row.contains(hgvs_g) && row.contains(term));
            let event = event.collect::<Vec<_>>();
            let last = i + 1 == annotated.len();
            assert_eq!(event.len(), usize::from(last), "{id}: {rows:?}");
            if last {
                assert!(
                    event[0].ends_with("|INFO_ONE_ROW_PER_GENE"),
                    "{id}: {rows:?}"
                );
            }
            assert!(
                rows.iter().all(|row| row.starts_with(&format!("{alt}|"))),
                "{id}: {rows:?}"
            );
        }
        Ok(())
    }

    /// A group without an event keeps the rows of its breakends only, also a pair whose
    /// `EVENTTYPE` is not the event of its shape. Mehari counts the records
    /// whose mates never arrived in one warning, also at a contig change, but not those of a
    /// translocation, whose group completes on the contig of its last record, nor a head-to-head
    /// pair without its tail-to-tail pair.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case::missing_mate(
        "17\t41195000\tbnd1\tG\tG[17:41280000[\t.\t.\tMATEID=gone\tGT\t0/1\n",
        1
    )]
    #[case::missing_mate_at_contig_change(
        "13\t32900000\tbnd1\tG\tG[13:32950000[\t.\t.\tMATEID=gone\tGT\t0/1\n\
         17\t41197701\tsnv\tG\tC\t.\t.\t.\tGT\t0/1\n",
        1
    )]
    #[case::translocation(
        "13\t32900000\tbnd2\tT\t]17:41195000]T\t.\t.\tMATEID=bnd1\tGT\t0/1\n\
         17\t41195000\tbnd1\tG\tG[13:32900000[\t.\t.\tMATEID=bnd2\tGT\t0/1\n",
        0
    )]
    #[case::head_to_head_pair(
        "17\t41195000\tbnd1\tG\tG]17:41280000]\t.\t.\tMATEID=bnd2\tGT\t0/1\n\
         17\t41280000\tbnd2\tT\tT]17:41195000]\t.\t.\tMATEID=bnd1\tGT\t0/1\n",
        0
    )]
    #[case::deletion_as_insertion(
        "17\t41195000\tbnd1\tG\tG[17:41280000[\t.\t.\tMATEID=bnd2;EVENTTYPE=INS\tGT\t0/1\n\
         17\t41280000\tbnd2\tT\t]17:41195000]T\t.\t.\tMATEID=bnd1;EVENTTYPE=INS\tGT\t0/1\n",
        0
    )]
    #[tokio::test]
    async fn breakend_group_without_event(
        #[case] records: &str,
        #[case] count: usize,
    ) -> Result<(), anyhow::Error> {
        let annotated = annotate_records(records, false).await?;
        let ids = records.lines().map(|line| line.split('\t').nth(2).unwrap());
        assert!(annotated.iter().map(|(id, _)| id.as_str()).eq(ids));
        for (id, rows) in annotated.iter().filter(|(id, _)| id.starts_with("bnd")) {
            assert!(!rows.is_empty(), "{id}");
            assert!(
                rows.iter().all(|row| is_breakend_row(row)),
                "{id}: {rows:?}"
            );
        }
        // The log buffer is shared by the cases, so only the presence of a warning is checked
        // here; `mates::test::unresolved_records` checks the count.
        if count > 0 {
            assert!(logs_contain(&format!(
                "{count} breakend record(s) whose MATEID or EVENT mates never arrived"
            )));
        }
        Ok(())
    }

    /// The breakends of an inversion around BRCA2 that loses a few bases at both breaks, and of
    /// a deletion of BRCA1 that takes in a fragment of chromosome 13, give the last record of
    /// each group the rows of the segments between their breaks: the inversion and the two lost
    /// pieces, and the deletion, each with one row per gene for the transcripts inside. The
    /// fragment, which stays in one piece, gives none, and the other records keep their own rows.
    #[tokio::test]
    async fn breakend_group_walks_its_segments() -> Result<(), anyhow::Error> {
        let records = "13\t32800000\tf1\tT\t]17:41195000]T\t.\t.\tEVENT=del\tGT\t0/1\n\
                       13\t32800100\tf2\tG\tG[17:41280000[\t.\t.\tEVENT=del\tGT\t0/1\n\
                       13\t32885000\ti1\tG\tG]13:32980000]\t.\t.\tEVENT=inv\tGT\t0/1\n\
                       13\t32885010\ti2\tA\t[13:32980005[A\t.\t.\tEVENT=inv\tGT\t0/1\n\
                       13\t32980000\ti3\tT\tT]13:32885000]\t.\t.\tEVENT=inv\tGT\t0/1\n\
                       13\t32980005\ti4\tC\t[13:32885010[C\t.\t.\tEVENT=inv\tGT\t0/1\n\
                       17\t41195000\td1\tG\tG[13:32800000[\t.\t.\tEVENT=del\tGT\t0/1\n\
                       17\t41280000\td2\tT\t]13:32800100]T\t.\t.\tEVENT=del\tGT\t0/1\n";
        let inv = "|g.32885010_32980000inv|";
        let lost = ["|g.32885001_32885009del|", "|g.32980001_32980004del|"];
        let del = "|g.41195001_41279999del|";
        let annotated = annotate_records(records, false).await?;
        assert_eq!(annotated.len(), 8);
        for (id, rows) in &annotated {
            let derived = (rows.iter()).filter(|row| !is_breakend_row(row));
            let derived = derived.collect::<Vec<_>>();
            let gene_row = |hgvs_g: &str, gene: &str| {
                (derived.iter()).any(|row| {
                    row.contains(hgvs_g)
                        && row.contains(gene)
                        && row.ends_with("|INFO_ONE_ROW_PER_GENE")
                })
            };
            match id.as_str() {
                "i4" => {
                    assert!(gene_row(inv, "|BRCA2|"), "{id}: {rows:?}");
                    assert!(
                        (derived.iter()).any(|row| row.contains(lost[0])),
                        "{id}: {rows:?}"
                    );
                    assert!(
                        (derived.iter()).all(|row| [inv, lost[0], lost[1]]
                            .iter()
                            .any(|hgvs_g| row.contains(hgvs_g))),
                        "{id}: {rows:?}"
                    );
                }
                "d2" => {
                    assert!(
                        gene_row(del, "transcript_ablation|HIGH|BRCA1|"),
                        "{id}: {rows:?}"
                    );
                    assert!(
                        derived.iter().all(|row| row.contains(del)),
                        "{id}: {rows:?}"
                    );
                }
                _ => assert!(derived.is_empty(), "{id}: {rows:?}"),
            }
        }
        Ok(())
    }

    /// The breakends of a contig inserted behind 41215930, whose mates lie on the contig `<ctg1>`
    /// of an assembly, VCF 4.5 5.4.2, keep the rows of their own junctions, without `HGVS.g`, as
    /// the reference lacks the contig, and get no event.
    #[tokio::test]
    async fn breakends_of_inserted_contig() -> Result<(), anyhow::Error> {
        let records = "17\t41215930\tins1\tG\tG[<ctg1>:1[\t.\t.\tEVENT=ins\tGT\t0/1\n\
                       17\t41215931\tins2\tA\t]<ctg1>:329]A\t.\t.\tEVENT=ins\tGT\t0/1\n";
        let annotated = annotate_records(records, false).await?;
        assert_eq!(annotated.len(), 2);
        for (id, rows) in &annotated {
            assert!(
                rows.iter().any(|row| row.contains("feature_truncation")),
                "{id}: {rows:?}"
            );
            assert!(
                rows.iter().all(|row| !row.contains("|g.")),
                "{id}: {rows:?}"
            );
        }
        Ok(())
    }

    /// The output keeps the order of the input, and only the last record of a group has the rows
    /// of its event. In compound mode, the rows of the breakends, also those of their event, say
    /// that they were not evaluated with the variants on their haplotype, as the rows of the SNV
    /// do.
    #[rstest::rstest]
    #[tokio::test]
    async fn breakend_group_keeps_order(
        #[values(false, true)] compound: bool,
    ) -> Result<(), anyhow::Error> {
        let records = "\
            17\t41195000\tbnd1\tG\tG[17:41280000[\t.\t.\tMATEID=bnd2\tGT\t0|1\n\
            17\t41197701\tsnv\tG\tC\t.\t.\t.\tGT\t0|1\n\
            17\t41280000\tbnd2\tT\t]17:41195000]T\t.\t.\tMATEID=bnd1\tGT\t0|1\n";
        let annotated = annotate_records(records, compound).await?;
        let ids = annotated
            .iter()
            .map(|(id, _)| id.as_str())
            .collect::<Vec<_>>();
        assert_eq!(ids, ["bnd1", "snv", "bnd2"]);
        for (id, rows) in &annotated {
            assert_eq!(
                rows.iter().any(|row| row.contains("|transcript_ablation|")),
                id == "bnd2",
                "{id}: {rows:?}"
            );
            assert!(!rows.is_empty(), "{id}");
            assert!(
                rows.iter().all(|row| row.contains(INCOMPLETE) == compound),
                "{id}: {rows:?}"
            );
        }
        Ok(())
    }

    /// In compound mode, a breakend reaches the transcripts up to its mate behind it on the same
    /// contig, so the rows of a variant on its haplotype say that it was not evaluated with the
    /// junction, also on a transcript that neither breakend reaches, and also when the mate never
    /// arrives, so that there is no event.
    #[rstest::rstest]
    #[case::pair(
        "MATEID=bnd2",
        "17\t41290000\tbnd2\tT\t]17:41180000]T\t.\t.\tMATEID=bnd1\tGT\t0|1\n"
    )]
    #[case::missing_mate("MATEID=gone", "")]
    #[tokio::test]
    async fn breakend_reaches_transcripts_up_to_its_mate(
        #[case] info: &str,
        #[case] mate: &str,
    ) -> Result<(), anyhow::Error> {
        let records = format!(
            "17\t41180000\tbnd1\tG\tG[17:41290000[\t.\t.\t{info}\tGT\t0|1\n\
             17\t41197701\tsnv\tG\tC\t.\t.\t.\tGT\t0|1\n\
             {mate}"
        );
        let annotated = annotate_records(&records, true).await?;
        let snv = &annotated[1].1;
        assert!(!snv.is_empty());
        assert!(snv.iter().all(|row| row.contains(INCOMPLETE)), "{snv:?}");
        Ok(())
    }

    /// The breakends of a deletion of all BRCA1 transcripts, downstream and upstream of them.
    const DELETION_BREAKENDS: &str = "\
        17\t41195000\tbnd1\tG\tG[17:41280000[\t.\t.\tMATEID=bnd2\tGT\t0/1\n\
        17\t41280000\tbnd2\tT\t]17:41195000]T\t.\t.\tMATEID=bnd1\tGT\t0/1\n";

    /// With `--pick`, a transcript keeps one row on a breakend record, and the last record of the
    /// group has its own picked rows and those of its event: the row of higher impact, here the
    /// `transcript_ablation` of the event's row of BRCA1, that of NM_007294.4, over the
    /// `downstream_gene_variant` or `upstream_gene_variant` of the breakend.
    #[rstest::rstest]
    #[tokio::test]
    async fn breakend_group_with_pick_keeps_one_row_per_transcript(
        #[values(PickUnit::Gene, PickUnit::Transcript)] unit: PickUnit,
    ) -> Result<(), anyhow::Error> {
        let predictor_settings = PredictorSettings {
            transcript_settings: TranscriptSettings::most_severe(unit),
            ..Default::default()
        };
        for (id, rows) in annotate_records_with(DELETION_BREAKENDS, predictor_settings).await? {
            let transcripts = (rows.iter())
                .map(|row| row.split('|').nth(7).unwrap_or_default())
                .collect::<HashSet<_>>();
            assert_eq!(transcripts.len(), rows.len(), "{id}: {rows:?}");
            let brca1 = rows.iter().filter(|row| row.contains("|NM_007294.4|"));
            assert!(
                brca1.clone().count() > 0
                    && brca1
                        .into_iter()
                        .all(|row| row.contains("transcript_ablation") == (id == "bnd2")),
                "{id}: {rows:?}"
            );
        }
        Ok(())
    }

    /// A plugin sees an event row of a breakend group as a row of the event: `variants()` gives
    /// the derived `<DEL>` at the `POS` of the event, `structural-variant()` its deletion and
    /// `breakend()` nothing. The plugin runs once for the group, and its values go to the event
    /// rows of the last record; the first record has none.
    #[cfg(feature = "wasm")]
    #[tokio::test]
    async fn plugin_sees_event_rows_as_event() -> Result<(), anyhow::Error> {
        let mut predictor_settings = PredictorSettings::default();
        predictor_settings.reporting_settings.plugins = vec![PluginSpec::Wasm {
            path: "tests/data/plugins/probe.wasm".into(),
            options: Vec::new(),
        }];
        let annotated = annotate_records_with(DELETION_BREAKENDS, predictor_settings).await?;
        // The event rows without their allele, which is that of the record.
        let event_rows = |rows: &[String]| {
            (rows.iter())
                .filter(|row| row.contains("|g.41195001_41279999del|"))
                .filter_map(|row| row.split_once('|').map(|(_, rest)| rest.to_string()))
                .collect::<Vec<_>>()
        };
        let (bnd1, bnd2) = (&annotated[0].1, &annotated[1].1);
        let event = event_rows(bnd2);
        assert!(!event.is_empty(), "{bnd2:?}");
        assert!(event_rows(bnd1).is_empty(), "{bnd1:?}");
        let view = "~variant-position^41195000~structural-variant^deletion:41279999::false\
                    ~breakend^~";
        assert!(event.iter().all(|row| row.contains(view)), "{event:?}");
        // The rows of the breakend itself see it.
        let own = (bnd2.iter()).filter(|row| row.contains("delins["));
        let view =
            "~variant-position^41280000~structural-variant^~breakend^false::17:41195000:false~";
        assert!(own.clone().count() > 0, "{bnd2:?}");
        assert!(own.into_iter().all(|row| row.contains(view)), "{bnd2:?}");
        Ok(())
    }

    /// An imprecise breakend has an uncertain g. description, and its rows say so.
    #[tokio::test]
    async fn imprecise_breakend_is_annotated() -> Result<(), anyhow::Error> {
        let records =
            "17\t41215930\tbnd\tG\tG]13:32900000]\t.\t.\tIMPRECISE;CIPOS=-500,500\tGT\t0/1\n";
        let records = annotate_records(records, false).await?;
        let bnd = &records[0].1;
        assert!(!bnd.is_empty());
        assert!(
            bnd.iter().all(|row| row
                .contains("|g.(41215431_41216431)_qterdelins[NC_000013.10:g.pter_32900000inv]|")
                && row.contains("INFO_IMPRECISE_BREAKPOINTS")),
            "{bnd:?}"
        );
        Ok(())
    }

    /// An event reaches the transcripts along its length, not only those at its position. The
    /// event starts about 96 kb in front of BRCA1 at chr17:41196312.
    #[rstest::rstest]
    #[case::svlen_into_gene("SVLEN=-100000", true)]
    #[case::svlen_short("SVLEN=-1000", false)]
    #[case::end_into_gene("END=41200000", true)]
    #[tokio::test]
    async fn event_reaches_transcripts_by_length(
        #[case] info: &str,
        #[case] marked: bool,
    ) -> Result<(), anyhow::Error> {
        let records = format!(
            "17\t41100000\tdel\tN\t<DEL>\t.\t.\t{info}\tGT\t1/1\n\
             17\t41197701\tsnv\tG\tC\t.\t.\t.\tGT\t1/1\n"
        );
        let records = annotate_records(&records, true).await?;
        let rows = &records[1].1;
        assert!(!rows.is_empty());
        assert!(
            rows.iter().all(|row| row.contains(INCOMPLETE) == marked),
            "{rows:?}"
        );
        Ok(())
    }

    /// `*`, `<*>` and `<NON_REF>` are no events, so the group is complete without them, nor are
    /// a `<CNV>` with the copy number of the reference, a total of 2 before VCF 4.4, and `<CN1>`.
    #[tokio::test]
    async fn no_variant_alleles_leave_groups_complete() -> Result<(), anyhow::Error> {
        let records = "\
            17\t41197701\tsnv1\tG\tC\t.\t.\t.\tGT\t1/1\n\
            17\t41197710\tstar\tG\t*\t.\t.\t.\tGT\t1/1\n\
            17\t41197720\tany\tG\t<*>\t.\t.\t.\tGT\t1/1\n\
            17\t41197730\tnon_ref\tG\t<NON_REF>\t.\t.\t.\tGT\t1/1\n\
            17\t41197740\tcn2\tG\t<CNV>\t.\t.\tEND=41197800;CN=2\tGT\t1/1\n\
            17\t41197750\tcn1\tG\t<CN1>\t.\t.\tEND=41197800\tGT\t1/1\n\
            17\t41215920\tsnv2\tG\tT\t.\t.\t.\tGT\t1/1\n";
        let records = annotate_records(records, true).await?;
        let rows = records
            .iter()
            .flat_map(|(_, rows)| rows)
            .collect::<Vec<_>>();
        assert!(rows.iter().any(|row| row.contains("|comp_1|")), "{rows:?}");
        assert!(rows.iter().all(|row| !row.contains(INCOMPLETE)), "{rows:?}");
        Ok(())
    }

    /// Only an event on the same haplotype as a variant marks its rows. The SNV has `0|1`.
    #[rstest::rstest]
    #[case::unphased("0/1", false)]
    #[case::other_haplotype("1|0", false)]
    #[case::same_haplotype("0|1", true)]
    #[case::homozygous("1/1", true)]
    #[tokio::test]
    async fn event_marks_rows_on_its_haplotype(
        #[case] genotype: &str,
        #[case] marked: bool,
    ) -> Result<(), anyhow::Error> {
        let records = format!(
            "17\t41197701\tsnv\tG\tC\t.\t.\t.\tGT\t0|1\n\
             17\t41197710\tdel\tN\t<DEL>\t.\t.\tSVLEN=-50\tGT\t{genotype}\n"
        );
        let records = annotate_records(&records, true).await?;
        let rows = &records[0].1;
        assert!(!rows.is_empty());
        assert!(
            rows.iter().all(|row| row.contains(INCOMPLETE) == marked),
            "{rows:?}"
        );
        Ok(())
    }

    /// A group that fails validation is not evaluated together, which the rows of its variants
    /// say, not only the log. The deletion covers the SNV.
    #[tokio::test]
    async fn invalid_group_marks_rows() -> Result<(), anyhow::Error> {
        let records = "\
            17\t41197700\tdel\tAG\tA\t.\t.\t.\tGT\t1/1\n\
            17\t41197701\tsnv\tG\tC\t.\t.\t.\tGT\t1/1\n";
        for (_, rows) in annotate_records(records, true).await? {
            assert!(!rows.is_empty());
            assert!(rows.iter().all(|row| !row.contains("|comp_")), "{rows:?}");
            assert!(rows.iter().all(|row| row.contains(INCOMPLETE)), "{rows:?}");
        }
        Ok(())
    }

    fn transcripts() -> Sources {
        Sources {
            transcripts: Some(vec!["tests/data/annotate/db/grch37/txs.bin.zst".into()]),
            ..Default::default()
        }
    }

    fn records(lines: &[String]) -> Vec<&String> {
        lines.iter().filter(|line| !line.starts_with('#')).collect()
    }

    /// The header then describes mehari's columns, so no input annotation may remain.
    async fn assert_input_field_is_replaced(key: &str) -> Result<(), anyhow::Error> {
        let lines = annotate_input_with_field(transcripts(), key, key, OutputFormat::Vcf).await?;
        let records = records(&lines);

        assert_eq!(records.len(), 3);
        assert!(!records[0].contains(&format!("{key}=")), "{}", records[0]);
        assert!(records[1].contains(&format!("{key}=T|")), "{}", records[1]);
        assert!(!records[1].contains("from_input"), "{}", records[1]);
        assert!(!records[2].contains(&format!("{key}=")), "{}", records[2]);
        Ok(())
    }

    #[tokio::test]
    #[tracing_test::traced_test]
    async fn input_ann_is_replaced_when_annotating_transcripts() -> Result<(), anyhow::Error> {
        assert_input_field_is_replaced("ANN").await?;
        assert!(logs_contain("The input already has the INFO field ANN"));
        Ok(())
    }

    #[tokio::test]
    #[tracing_test::traced_test]
    async fn input_field_under_chosen_key_is_replaced() -> Result<(), anyhow::Error> {
        assert_input_field_is_replaced("MEHARI").await?;
        assert!(logs_contain("The input already has the INFO field MEHARI"));
        Ok(())
    }

    /// With another key, an `ANN` from another tool stays next to mehari's annotations.
    #[tokio::test]
    #[tracing_test::traced_test]
    async fn input_ann_is_kept_with_another_key() -> Result<(), anyhow::Error> {
        let lines =
            annotate_input_with_field(transcripts(), "ANN", "MEHARI", OutputFormat::Vcf).await?;
        let records = records(&lines);

        assert!(!logs_contain("The input already has the INFO field"));
        let header_line = |prefix: &str| lines.iter().find(|line| line.starts_with(prefix));
        let ann = header_line("##INFO=<ID=ANN,").unwrap();
        assert!(!ann.contains("Source="), "{ann}");
        let mehari = header_line("##INFO=<ID=MEHARI,").unwrap();
        assert!(mehari.contains(",Source=\"mehari\""), "{mehari}");

        assert_eq!(records.len(), 3);
        for record in &records {
            assert!(record.contains("from_input"), "{record}");
        }
        assert!(!records[0].contains("MEHARI="), "{}", records[0]);
        assert!(records[1].contains("MEHARI=T|"), "{}", records[1]);
        Ok(())
    }

    /// Without transcripts, mehari leaves the `ANN` header and the annotations alone.
    #[tokio::test]
    #[tracing_test::traced_test]
    async fn input_ann_is_kept_without_transcripts() -> Result<(), anyhow::Error> {
        let sources = Sources {
            frequencies: Some(vec!["tests/data/annotate/db/grch37/seqvars/freqs".into()]),
            ..Default::default()
        };
        let lines = annotate_input_with_field(sources, "ANN", "ANN", OutputFormat::Vcf).await?;
        let records = records(&lines);

        assert!(!logs_contain("The input already has the INFO field"));
        assert_eq!(records.len(), 3);
        for record in &records {
            assert!(record.contains("from_input"), "{record}");
        }
        Ok(())
    }

    /// `vcfInfo` of JSONL leaves out the field of the annotations, and only that one.
    #[rstest::rstest]
    #[case::default_key("ANN", false)]
    #[case::other_key("MEHARI", true)]
    #[tokio::test]
    async fn jsonl_info_keeps_fields_under_other_keys(
        #[case] ann_key: &str,
        #[case] ann_in_info: bool,
    ) -> Result<(), anyhow::Error> {
        let lines =
            annotate_input_with_field(transcripts(), "ANN", ann_key, OutputFormat::Jsonl).await?;

        assert_eq!(lines.len(), 3);
        for line in &lines {
            let record: serde_json::Value = serde_json::from_str(line)?;
            let info = record["vcfInfo"].as_object().unwrap();
            assert_eq!(info.contains_key("ANN"), ann_in_info, "{line}");
            assert!(!info.contains_key("MEHARI"), "{line}");
        }
        Ok(())
    }

    /// The default keeps `ANN`, and other keys must match the VCF specification's pattern.
    #[rstest::rstest]
    #[case("ANN", true)]
    #[case("_ann.2", true)]
    #[case("1000G", true)]
    #[case("", false)]
    #[case("1ANN", false)]
    #[case(".ANN", false)]
    #[case("A-B", false)]
    #[case("A B", false)]
    #[case("A=B", false)]
    #[case("A;B", false)]
    fn ann_key_follows_vcf_specification(#[case] key: &str, #[case] valid: bool) {
        let args = Args::try_parse_from([
            "seqvars",
            "--input",
            "in.vcf",
            "--transcripts",
            "txs",
            "--ann-key",
            key,
        ]);
        assert_eq!(args.is_ok(), valid, "{key}");
    }

    #[test]
    fn ann_key_defaults_to_ann() -> Result<(), anyhow::Error> {
        assert_eq!(
            Args::try_parse_from(["seqvars", "--input", "in.vcf", "--transcripts", "txs"])?.ann_key,
            "ANN"
        );
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
            ann_key: "ANN".into(),
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
            ann_key: "ANN".into(),
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
            ann_key: "ANN".into(),
            sources: Sources {
                transcripts: Some(vec!["tests/data/annotate/db/grch37/txs.bin.zst".into()]),
                ..Default::default()
            },
        };

        run(&args_common, &args).await?;

        let header = open_variant_reader(&args.output)
            .await?
            .read_header()
            .await?;
        let names = AnnField::ann_field_names_from_header(&header, "ANN").unwrap_or_default();
        let compound_columns = read_vcf(&args.output)
            .await?
            .iter()
            .flat_map(|record| match record.info().get("ANN").flatten() {
                Some(Value::Array(Array::String(anns))) => anns.clone(),
                _ => vec![],
            })
            .flatten()
            .map(|ann| -> Result<String, anyhow::Error> {
                let ann = AnnField::parse_named(&ann, &names)?;
                let column = |name: &str| ann.custom_fields[name].clone().unwrap_or_default();
                Ok(format!(
                    "{} {} {}",
                    ann.feature_id,
                    column(ANN_COMPOUND_IDS),
                    column(ANN_COMPOUND_VARIANTS)
                ))
            })
            .collect::<Result<Vec<_>, anyhow::Error>>()?;
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
            ann_key: "ANN".into(),
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
            ann_key: "ANN".into(),
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
            ann_key: "ANN".into(),
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
            ann_key: "ANN".into(),
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

    /// Annotate a chrM SNV in a VCF with the header line `contig` (or none), with `reference`.
    async fn annotate_chrm(contig: &str, reference: Option<PathBuf>) -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_in = temp.join("input.vcf");
        std::fs::write(
            &path_in,
            format!(
                "##fileformat=VCFv4.2\n{contig}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 chrM\t3\t.\tT\tC\t.\t.\t.\n"
            ),
        )?;
        let args = Args {
            threads: 1,
            reference,
            in_memory_reference: true,
            assembly: Some("grch38".into()),
            predictor_settings: Default::default(),
            input: path_in.to_string_lossy().into(),
            output: temp.join("output.vcf").to_string_lossy().into(),
            output_format: OutputFormat::Vcf,
            max_var_count: None,
            ann_key: "ANN".into(),
            sources: Sources {
                transcripts: Some(vec![
                    "tests/data/annotate/db/grch38/mt-origin.txs.bin.zst".into(),
                ]),
                ..Default::default()
            },
        };
        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        run(&args_common, &args).await
    }

    /// A chrM of 16,571 bp in the VCF header is the one of hg19, NC_001807.4.
    #[tokio::test]
    #[tracing_test::traced_test]
    async fn hg19_chrm_in_the_vcf_header_is_logged() -> Result<(), anyhow::Error> {
        annotate_chrm("##contig=<ID=chrM,length=16571>\n", None).await?;

        assert!(logs_contain("chrM has 16571 bp in the VCF header"));
        Ok(())
    }

    /// A chrM of 16,571 bp in the reference genome is the one of hg19, NC_001807.4.
    #[tokio::test]
    #[tracing_test::traced_test]
    async fn hg19_chrm_in_the_reference_is_logged() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let reference = temp.join("chrM.fasta");
        let rcrs = std::fs::read_to_string("tests/data/annotate/db/grch38/chrM.fasta")?;
        std::fs::write(&reference, format!("{}AC\n", rcrs.trim_end()))?;

        annotate_chrm("", Some(reference)).await?;

        assert!(logs_contain("chrM has 16571 bp in the reference genome"));
        Ok(())
    }

    /// The rCRS, NC_012920.1, has 16,569 bp.
    #[tokio::test]
    #[tracing_test::traced_test]
    async fn rcrs_chrm_is_not_logged() -> Result<(), anyhow::Error> {
        annotate_chrm(
            "##contig=<ID=chrM,length=16569>\n",
            Some("tests/data/annotate/db/grch38/chrM.fasta".into()),
        )
        .await?;

        assert!(!logs_contain("NC_001807.4"));
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
            ann_key: "ANN".into(),
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
            ann_key: "ANN".into(),
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
