//! Annotation of sequence variants.
use std::collections::{HashMap, HashSet};
use std::env;
use std::fmt::Display;
use std::fs::File;
use std::io::{Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;
use std::time::Instant;

use self::ann::{AnnField, FeatureBiotype};
use crate::annotate::cli::{Sources, TranscriptSettings};
use crate::annotate::genotype_string;
use crate::annotate::seqvars::csq::{
    ConfigBuilder as ConsequencePredictorConfigBuilder, ConsequencePredictor, VcfVariant,
};
use crate::annotate::seqvars::provider::{
    ConfigBuilder as MehariProviderConfigBuilder, Provider as MehariProvider,
};
use crate::common::noodles::{open_variant_reader, open_variant_writer, NoodlesVariantReader};
use crate::common::{guess_assembly, GenomeRelease};
use crate::db::merge::merge_transcript_databases;
use crate::pbs::txs::TxSeqDatabase;
use crate::ped::{PedigreeByName, Sex};
use annonars::common::cli::is_canonical;
use annonars::common::keys;
use annonars::freqs::serialized::{auto, mt, xy};
use anyhow::{anyhow, Error};
use biocommons_bioutils::assemblies::Assembly;
use clap::{Args as ClapArgs, Parser};
use flate2::write::GzEncoder;
use flate2::Compression;
use itertools::Itertools;
use noodles::vcf::header::record::value::map::format::Number as FormatNumber;
use noodles::vcf::header::record::value::map::format::Type as FormatType;
use noodles::vcf::header::record::value::map::info::Number;
use noodles::vcf::header::record::value::map::{info::Type as InfoType, Info};
use noodles::vcf::header::FileFormat;
use noodles::vcf::io::writer::Writer as VcfWriter;
use noodles::vcf::variant::io::Write as VcfWrite;
use noodles::vcf::variant::record::samples::keys::key::{
    CONDITIONAL_GENOTYPE_QUALITY, GENOTYPE, READ_DEPTH,
};
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record_buf::info::field;
use noodles::vcf::variant::RecordBuf as VcfRecord;
use noodles::vcf::{header::record::value::map::Map, Header as VcfHeader};
use once_cell::sync::Lazy;
use prost::Message;
use rocksdb::{DBWithThreadMode, MultiThreaded};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use strum::Display;
use thousands::Separable;
use tokio::io::AsyncWriteExt;

pub mod ann;
pub mod binning;
pub mod csq;
pub mod provider;
pub(crate) mod reference;

/// Parsing of HGNC xlink records.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct HgncRecord {
    /// HGNC ID.
    pub hgnc_id: String,

    /// Ensembl gene ID.
    pub ensembl_gene_id: String,

    /// Entrez/NCBI gene ID.
    pub entrez_id: String,

    /// HGNC approved gene symbol.
    #[serde(alias = "symbol")]
    pub gene_symbol: String,
}

/// Command line arguments for `annotate seqvars` sub command.
#[derive(Parser, Debug)]
#[command(about = "Annotate sequence variant VCF files", long_about = None)]
pub struct Args {
    /// Genome release to use, default is to auto-detect.
    #[arg(long, value_enum)]
    pub genome_release: Option<GenomeRelease>,

    /// Reference genome FASTA file (with accompanying index).
    #[arg(long)]
    pub reference: Option<PathBuf>,

    /// Read the reference genome into memory.
    #[arg(long, requires = "reference")]
    pub in_memory_reference: bool,

    /// Path to the input PED file.
    #[arg(long)]
    pub path_input_ped: Option<String>,

    /// Path to the input VCF file.
    #[arg(long)]
    pub path_input_vcf: String,

    #[command(flatten)]
    pub output: PathOutput,

    /// For debug purposes, maximal number of variants to annotate.
    #[arg(long)]
    pub max_var_count: Option<usize>,

    /// Path to HGNC TSV file.
    #[arg(long, required_unless_present = "path_output_vcf")]
    pub hgnc: Option<String>,

    /// What to annotate and which source to use.
    #[command(flatten)]
    pub sources: Sources,

    /// Transcript annotation related settings
    #[command(flatten)]
    pub transcript_settings: TranscriptSettings,
}

#[derive(Debug, Display, Copy, Clone, clap::ValueEnum, PartialEq, Eq, parse_display::FromStr)]
pub enum AnnotationOption {
    Transcripts,
    Frequencies,
    ClinVar,
}

/// Command line arguments to enforce either `--path-output-vcf` or `--path-output-tsv`.
#[derive(Debug, ClapArgs)]
#[group(required = true, multiple = false)]
pub struct PathOutput {
    /// Path to the output VCF file.
    #[arg(long)]
    pub path_output_vcf: Option<String>,

    /// Path to the output TSV file (for import into VarFish).
    #[arg(long, requires = "hgnc")]
    pub path_output_tsv: Option<String>,
}

fn build_header(
    header_in: &VcfHeader,
    with_annotations: bool,
    with_frequencies: bool,
    with_clinvar: bool,
    additional_records: &[(String, String)],
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
        let fields = AnnField::ann_field_names().join(" | ");
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

pub static CHROM_MT: Lazy<HashSet<&'static str>> =
    Lazy::new(|| HashSet::from_iter(["M", "MT", "chrM", "chrMT"]));
pub static CHROM_XY: Lazy<HashSet<&'static str>> =
    Lazy::new(|| HashSet::from_iter(["X", "Y", "chrX", "chrY"]));

pub static CHROM_AUTO: Lazy<HashSet<&'static str>> = Lazy::new(|| {
    HashSet::from_iter([
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
        "17", "18", "19", "20", "21", "22", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
        "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
        "chr18", "chr19", "chr20", "chr21", "chr22",
    ])
});

/// Mapping from chromosome name to chromsome number.
pub static CHROM_TO_CHROM_NO: Lazy<HashMap<String, u32>> = Lazy::new(|| {
    let mut m = HashMap::new();

    for i in 1..=22 {
        m.insert(format!("chr{}", i), i);
        m.insert(format!("{}", i), i);
    }
    m.insert(String::from("X"), 23);
    m.insert(String::from("chrX"), 23);
    m.insert(String::from("Y"), 24);
    m.insert(String::from("chrY"), 24);
    m.insert(String::from("M"), 25);
    m.insert(String::from("chrM"), 25);
    m.insert(String::from("MT"), 25);
    m.insert(String::from("chrMT"), 25);

    m
});

/// Mapping from chromosome number to canonical chromosome name.
///
/// We use the names without `"chr"` prefix for the canonical name.  In the case of GRCh38,
/// the the prefix must be prepended.
pub static CHROM_NO_TO_NAME: Lazy<Vec<String>> = Lazy::new(|| {
    let mut v = Vec::new();
    v.push(String::from("")); // 0
    for i in 1..=22 {
        v.push(format!("{}", i));
    }
    v.push(String::from("X"));
    v.push(String::from("Y"));
    v.push(String::from("MT"));
    v
});

/// Return path component for the assembly.
pub fn path_component(assembly: Assembly) -> &'static str {
    match assembly {
        Assembly::Grch37 | Assembly::Grch37p10 => "grch37",
        Assembly::Grch38 => "grch38",
    }
}

/// Load protobuf transcripts.
pub fn load_tx_db(tx_path: impl AsRef<Path> + Display) -> anyhow::Result<TxSeqDatabase> {
    // Open file and if necessary, wrap in a decompressor.
    let file = File::open(tx_path.as_ref())
        .map_err(|e| anyhow!("failed to open file {}: {}", tx_path, e))?;
    let mut reader: Box<dyn Read> = match tx_path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gz") => Box::new(flate2::read::MultiGzDecoder::new(file)),
        Some("zst") => Box::new(
            zstd::Decoder::new(file)
                .map_err(|e| anyhow!("failed to open zstd decoder for {}: {}", tx_path, e))?,
        ),
        _ => Box::new(file),
    };

    // Now read the whole file into a byte buffer.
    let mut buffer = Vec::new();
    reader
        .read_to_end(&mut buffer)
        .map_err(|e| anyhow!("failed to read file {}: {}", tx_path, e))?;

    // Deserialize the buffer with prost.
    TxSeqDatabase::decode(&mut Cursor::new(buffer))
        .map_err(|e| anyhow!("failed to decode protobuf file {}: {}", tx_path, e))
}

/// Trait for writing out annotated VCF records as VCF or VarFish TSV.
///
/// TODO: use async_trait crate
pub trait AsyncAnnotatedVariantWriter {
    #[allow(async_fn_in_trait)]
    async fn write_noodles_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error>;
    #[allow(async_fn_in_trait)]
    async fn write_noodles_record(
        &mut self,
        header: &VcfHeader,
        record: &VcfRecord,
    ) -> Result<(), anyhow::Error>;

    #[allow(async_fn_in_trait)]
    async fn shutdown(&mut self) -> Result<(), anyhow::Error>;

    fn set_hgnc_map(&mut self, _hgnc_map: FxHashMap<String, HgncRecord>) {
        // nop
    }
    fn set_assembly(&mut self, _assembly: Assembly) {
        // nop
    }
    fn set_pedigree(&mut self, _pedigree: &PedigreeByName) {
        // nop
    }
}

/// Implement `AnnotatedVcfWriter` for `VcfWriter`.
impl<Inner: Write> AsyncAnnotatedVariantWriter for VcfWriter<Inner> {
    async fn write_noodles_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error> {
        self.write_header(header)
            .map_err(|e| anyhow::anyhow!("Error writing VCF header: {}", e))
    }

    async fn write_noodles_record(
        &mut self,
        header: &VcfHeader,
        record: &VcfRecord,
    ) -> Result<(), anyhow::Error> {
        self.write_variant_record(header, record)
            .map_err(|e| anyhow::anyhow!("Error writing VCF record: {}", e))
    }

    async fn shutdown(&mut self) -> Result<(), Error> {
        Ok(<VcfWriter<Inner>>::get_mut(self).flush()?)
    }
}

/// Implement `AsyncAnnotatedVariantWriter` for `AsyncWriter`.
impl<Inner: tokio::io::AsyncWrite + Unpin> AsyncAnnotatedVariantWriter
    for noodles::vcf::AsyncWriter<Inner>
{
    async fn write_noodles_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error> {
        self.write_header(header)
            .await
            .map_err(|e| anyhow::anyhow!("Error writing VCF header: {}", e))
    }

    async fn write_noodles_record(
        &mut self,
        header: &VcfHeader,
        record: &VcfRecord,
    ) -> Result<(), anyhow::Error> {
        noodles::vcf::AsyncWriter::write_variant_record(self, header, record)
            .await
            .map_err(|e| anyhow::anyhow!("Error writing VCF record: {}", e))
    }

    async fn shutdown(&mut self) -> Result<(), Error> {
        Ok(<noodles::vcf::AsyncWriter<Inner>>::get_mut(self)
            .flush()
            .await?)
    }
}

/// Implement `AsyncAnnotatedVariantWriter` for `AsyncWriter`.
impl<Inner: tokio::io::AsyncWrite + Unpin> AsyncAnnotatedVariantWriter
    for noodles::bcf::AsyncWriter<Inner>
{
    async fn write_noodles_header(&mut self, header: &VcfHeader) -> Result<(), Error> {
        self.write_header(header)
            .await
            .map_err(|e| anyhow::anyhow!("Error writing VCF header: {}", e))
    }

    async fn write_noodles_record(
        &mut self,
        header: &VcfHeader,
        record: &VcfRecord,
    ) -> Result<(), anyhow::Error> {
        noodles::bcf::AsyncWriter::write_variant_record(self, header, record)
            .await
            .map_err(|e| anyhow::anyhow!("Error writing VCF record: {}", e))
    }

    async fn shutdown(&mut self) -> Result<(), Error> {
        Ok(<noodles::bcf::AsyncWriter<Inner>>::get_mut(self)
            .shutdown()
            .await?)
    }
}

/// Writing of sequence variants to VarFish TSV files.
struct VarFishSeqvarTsvWriter {
    inner: Box<dyn Write>,
    assembly: Option<Assembly>,
    pedigree: Option<PedigreeByName>,
    header: Option<VcfHeader>,
    hgnc_map: Option<FxHashMap<String, HgncRecord>>,
}

/// Entry with genotype (`gt`), coverage (`dp`), allele depth (`ad`), somatic quality (`sq`) and
/// genotype quality (`gq`).
#[derive(Debug, Default)]
struct GenotypeInfo {
    pub name: String,
    pub gt: Option<String>,
    pub dp: Option<i32>,
    pub ad: Option<i32>,
    pub gq: Option<i32>,
    pub sq: Option<f32>,
}

#[derive(Debug, Default)]
struct GenotypeCalls {
    pub entries: Vec<GenotypeInfo>,
}

impl GenotypeCalls {
    /// Generate and return the dict in Postgres JSON syntax.
    ///
    /// The returned string is suitable for a direct TSV import into Postgres.
    pub fn for_tsv(&self) -> String {
        let mut result = String::new();
        result.push('{');

        let mut first = true;
        for entry in &self.entries {
            if first {
                first = false;
            } else {
                result.push(',');
            }
            result.push_str(&format!("\"\"\"{}\"\"\":{{", entry.name));

            let mut prev = false;
            if let Some(gt) = &entry.gt {
                prev = true;
                result.push_str(&format!("\"\"\"gt\"\"\":\"\"\"{}\"\"\"", gt));
            }

            if let Some(ad) = &entry.ad {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"ad\"\"\":{}", ad));
            }

            if let Some(dp) = &entry.dp {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"dp\"\"\":{}", dp));
            }

            // The DRAGEN variant caller writes out FORMAT/SQ for chrMT ("somatic quality") as
            // it uses the somatic genotyping model for mitochondrial callers.  We just write
            // this out as GQ for now.
            if let Some(gq) = &entry.gq {
                if prev {
                    result.push(',');
                }
                // prev = true;
                result.push_str(&format!("\"\"\"gq\"\"\":{}", gq));
            } else if let Some(sq) = &entry.sq {
                if prev {
                    result.push(',');
                }
                // prev = true;
                result.push_str(&format!("\"\"\"gq\"\"\":{}", sq.round() as i32));
            }

            result.push('}');
        }

        result.push('}');
        result
    }
}

impl VarFishSeqvarTsvWriter {
    // Create new TSV writer from path.
    pub fn with_path<P>(p: P) -> Self
    where
        P: AsRef<Path>,
    {
        Self {
            inner: if p.as_ref().extension().unwrap_or_default() == "gz" {
                Box::new(GzEncoder::new(
                    File::create(p).unwrap(),
                    Compression::default(),
                ))
            } else {
                Box::new(File::create(p).unwrap())
            },
            hgnc_map: None,
            assembly: None,
            pedigree: None,
            header: None,
        }
    }

    /// Flush buffers.
    pub fn flush(&mut self) -> Result<(), anyhow::Error> {
        self.inner.flush()?;
        Ok(())
    }

    /// Fill `record` coordinate fields.
    ///
    /// # Returns
    ///
    /// Returns `true` if the record is to be written to the TSV file and `false` if
    /// it is to be skipped because it was not on a canonical chromosome (chr1..chr22,
    /// chrX, chrY, chrM/chrMT).
    fn fill_coords(
        &self,
        record: &VcfRecord,
        tsv_record: &mut VarFishSeqvarTsvRecord,
    ) -> Result<bool, anyhow::Error> {
        let assembly = self.assembly.expect("assembly must have been set");
        tsv_record.release = match assembly {
            Assembly::Grch37 | Assembly::Grch37p10 => String::from("GRCh37"),
            Assembly::Grch38 => String::from("GRCh38"),
        };
        let name = record.reference_sequence_name();
        // strip "chr" prefix from `name`
        let name = if let Some(stripped) = name.strip_prefix("chr") {
            stripped
        } else {
            name
        };
        // add back "chr" prefix if necessary (for GRCh38)
        tsv_record.chromosome = if assembly == Assembly::Grch38 {
            format!("chr{}", name)
        } else {
            name.to_string()
        };

        if let Some(chromosome_no) = CHROM_TO_CHROM_NO.get(&tsv_record.chromosome) {
            tsv_record.chromosome_no = *chromosome_no;
        } else {
            return Ok(false);
        }

        tsv_record.reference = record.reference_bases().to_string();
        tsv_record.alternative = record.alternate_bases().as_ref()[0].to_string();

        tsv_record.start = record
            .variant_start()
            .expect("Telomeres unsupported")
            .into();
        tsv_record.end = tsv_record.start + tsv_record.reference.len() - 1;
        tsv_record.bin =
            binning::bin_from_range(tsv_record.start as i32 - 1, tsv_record.end as i32)? as u32;

        tsv_record.var_type = if tsv_record.reference.len() == tsv_record.alternative.len() {
            if tsv_record.reference.len() == 1 {
                String::from("snv")
            } else {
                String::from("mnv")
            }
        } else {
            String::from("indel")
        };

        Ok(true)
    }

    /// Fill `record` genotype and family-local allele frequencies.
    fn fill_genotype_and_freqs(
        &self,
        record: &VcfRecord,
        tsv_record: &mut VarFishSeqvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        use noodles::vcf::variant::record_buf::samples::sample::value::Array;
        use noodles::vcf::variant::record_buf::samples::sample::Value;

        // Extract genotype information.
        let hdr = self
            .header
            .as_ref()
            .expect("VCF header must be set/written");
        let file_format_gt = FileFormat::new(4, 3);
        let mut gt_calls = GenotypeCalls::default();
        let samples = record.samples();
        let sample_names = hdr.sample_names().iter();

        // The following couple of lines are quite ugly,
        // there's probably a more elegant way to do this
        let genotypes = samples.select(GENOTYPE);
        let get_genotype = |sample_idx| {
            genotypes.as_ref().and_then(|gt| {
                gt.get(sample_idx).map(|value| match value {
                    Some(Value::String(s)) => s.to_owned(),
                    Some(Value::Genotype(gt)) => genotype_string(gt, file_format_gt),
                    _ => ".".into(),
                })
            })
        };

        let read_depths = samples.select(READ_DEPTH);
        let get_read_depth = |sample_idx| {
            read_depths.as_ref().and_then(|rd| {
                rd.get(sample_idx)
                    .map(|value| match value {
                        Some(Value::Integer(i)) => Ok(*i),
                        None => Ok(-1),
                        _ => anyhow::bail!(format!("invalid DP value {:?}", value)),
                    })
                    .transpose()
                    .unwrap()
            })
        };

        let allele_depths = samples.select("AD");
        let get_allele_depths = |sample_idx| {
            allele_depths.as_ref().and_then(|ad| {
                ad.get(sample_idx)
                    .map(|value| match value {
                        Some(Value::Array(Array::Integer(arr))) => Ok(arr[1].unwrap_or(0)),
                        None => Ok(-1),
                        _ => anyhow::bail!(format!("invalid AD value {:?}", value)),
                    })
                    .transpose()
                    .unwrap()
            })
        };

        let conditional_gt_quality = samples.select(CONDITIONAL_GENOTYPE_QUALITY);
        let get_conditional_gt_quality = |sample_idx| {
            conditional_gt_quality.as_ref().and_then(|cgq| {
                cgq.get(sample_idx)
                    .map(|value| match value {
                        Some(Value::Integer(i)) => Ok(*i),
                        None => Ok(-1),
                        _ => anyhow::bail!(format!(
                            "invalid GQ value {:?} at {:?}:{:?}",
                            value,
                            record.reference_sequence_name(),
                            record.variant_start()
                        )),
                    })
                    .transpose()
                    .unwrap()
            })
        };

        let sq = samples.select("SQ");
        let get_sq = |sample_idx| {
            sq.as_ref().and_then(|sq| {
                sq.get(sample_idx)
                    .map(|value| match value {
                        Some(Value::Float(f)) => Ok(*f),
                        Some(Value::Array(Array::Float(f))) => {
                            Ok(f[0].expect("SQ should be a single float value"))
                        }
                        None => Ok(-1.0),
                        _ => {
                            anyhow::bail!(format!(
                                "invalid SQ value {:?} at {:?}:{:?}",
                                value,
                                record.reference_sequence_name(),
                                record.variant_start()
                            ))
                        }
                    })
                    .transpose()
                    .unwrap()
            })
        };

        for (i, name) in sample_names.enumerate() {
            let mut gt_info = GenotypeInfo {
                name: name.clone(),
                ..Default::default()
            };

            if let Some(gt) = get_genotype(i) {
                let individual = self
                    .pedigree
                    .as_ref()
                    .expect("pedigree must be set")
                    .individuals
                    .get(name)
                    .unwrap_or_else(|| panic!("individual {} not found in pedigree", name));
                // Update per-family counts.
                if ["X", "chrX"]
                    .iter()
                    .any(|c| *c == tsv_record.chromosome.as_str())
                {
                    match individual.sex {
                        Sex::Male => {
                            if gt.contains('1') {
                                tsv_record.num_hemi_alt += 1;
                            } else {
                                tsv_record.num_hemi_ref += 1;
                            }
                        }
                        Sex::Female | Sex::Unknown => {
                            // assume diploid/female if unknown
                            let matches_1 = gt.matches('1').count();
                            let matches_0 = gt.matches('0').count();
                            if matches_0 == 2 {
                                tsv_record.num_hom_ref += 1;
                            } else if matches_1 == 1 {
                                tsv_record.num_het += 1;
                            } else if matches_1 == 2 {
                                tsv_record.num_hom_alt += 1;
                            }
                        }
                    }
                } else if ["Y", "chrY"]
                    .iter()
                    .any(|c| *c == tsv_record.chromosome.as_str())
                {
                    if individual.sex == Sex::Male {
                        if gt.contains('1') {
                            tsv_record.num_hemi_alt += 1;
                        } else if gt.contains('0') {
                            tsv_record.num_hemi_ref += 1;
                        }
                    }
                } else {
                    let matches_1 = gt.matches('1').count();
                    let matches_0 = gt.matches('0').count();
                    if matches_0 == 2 {
                        tsv_record.num_hom_ref += 1;
                    } else if matches_1 == 1 {
                        tsv_record.num_het += 1;
                    } else if matches_1 == 2 {
                        tsv_record.num_hom_alt += 1;
                    }
                }

                // Store genotype value.
                gt_info.gt = Some(gt);
            }

            if let Some(dp) = get_read_depth(i) {
                gt_info.dp = Some(dp);
            }

            if let Some(ad) = get_allele_depths(i) {
                gt_info.ad = Some(ad);
            }

            if let Some(gq) = get_conditional_gt_quality(i) {
                gt_info.gq = Some(gq);
            }
            if let Some(sq) = get_sq(i) {
                gt_info.sq = Some(sq);
            }

            gt_calls.entries.push(gt_info);
        }
        tsv_record.genotype = gt_calls.for_tsv();

        Ok(())
    }

    /// Fill `record` background frequencies.
    fn fill_bg_freqs(
        &self,
        record: &VcfRecord,
        tsv_record: &mut VarFishSeqvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        // Extract gnomAD frequencies.
        let gnomad_exomes_an = record
            .info()
            .get("gnomad_exomes_an")
            .unwrap_or_default()
            .map(|v| match v {
                field::Value::Integer(value) => *value,
                _ => panic!("Unexpected value type for GNOMAD_EXOMES_AN"),
            })
            .unwrap_or_default();
        if gnomad_exomes_an > 0 {
            tsv_record.gnomad_exomes_homozygous = record
                .info()
                .get("gnomad_exomes_hom")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_EXOMES_HOM"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_exomes_heterozygous = record
                .info()
                .get("gnomad_exomes_het")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_EXOMES_HET"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_exomes_hemizygous = record
                .info()
                .get("gnomad_exomes_hemi")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_EXOMES_HEMI"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_exomes_frequency = (tsv_record.gnomad_exomes_hemizygous
                + tsv_record.gnomad_exomes_heterozygous
                + tsv_record.gnomad_exomes_homozygous * 2)
                as f32
                / gnomad_exomes_an as f32;
        }

        let gnomad_genomes_an = record
            .info()
            .get("gnomad_genomes_an")
            .unwrap_or_default()
            .map(|v| match v {
                field::Value::Integer(value) => *value,
                _ => panic!("Unexpected value type for GNOMAD_GENOMES_AN"),
            })
            .unwrap_or_default();
        if gnomad_genomes_an > 0 {
            tsv_record.gnomad_genomes_homozygous = record
                .info()
                .get("gnomad_genomes_hom")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_GENOMES_HOM"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_genomes_heterozygous = record
                .info()
                .get("gnomad_genomes_het")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_GENOMES_HET"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_genomes_hemizygous = record
                .info()
                .get("gnomad_genomes_hemi")
                .unwrap_or(Some(&field::Value::Integer(0)))
                .map(|v| match v {
                    field::Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_GENOMES_HEMI"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_genomes_frequency = (tsv_record.gnomad_genomes_hemizygous
                + tsv_record.gnomad_genomes_heterozygous
                + tsv_record.gnomad_genomes_homozygous * 2)
                as f32
                / gnomad_genomes_an as f32;
        }

        Ok(())
    }

    /// Fill `record` RefSEq and ENSEMBL and fields and write to `self.inner`.
    ///
    /// First, the values in the `ANN` field are parsed and all predictions are extracted.
    /// For each gene, we find the single worst prediction for a RefSeq transcript and the
    /// single worst prediction for an ENSEMBL transcript. We then combine these to create
    /// one output row per gene.
    fn expand_refseq_ensembl_and_write(
        &mut self,
        record: &VcfRecord,
        tsv_record: &mut VarFishSeqvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        if let Some(anns) = record
            .info()
            .get("ANN")
            .unwrap_or_default()
            .map(|v| match v {
                field::Value::Array(field::value::Array::String(values)) => values
                    .iter()
                    .filter(|v| v.is_some())
                    .map(|v| AnnField::from_str(v.as_ref().unwrap())),
                _ => panic!("Unexpected value type for INFO/ANN"),
            })
        {
            // Extract `AnnField` records, letting the errors bubble up.
            let anns: Result<Vec<_>, _> = anns.collect();
            let anns = anns?;

            // Collect the `AnnField` records by `gene_id`.
            let mut anns_by_gene: FxHashMap<String, Vec<AnnField>> = FxHashMap::default();
            for ann in anns {
                let gene_id = ann.gene_id.clone();
                anns_by_gene.entry(gene_id).or_default().push(ann);
            }

            // Within each gene's list of annotations, sort by consequence severity (worst first).
            for anns in anns_by_gene.values_mut() {
                anns.sort_by_key(|ann| ann.consequences[0]);
            }

            let empty_hgnc_record = HgncRecord {
                hgnc_id: "".to_string(),
                ensembl_gene_id: "".to_string(),
                entrez_id: "".to_string(),
                gene_symbol: "".to_string(),
            };

            fn is_ensembl_id(feature_id: &str) -> bool {
                feature_id.starts_with("ENST")
            }

            fn is_refseq_id(feature_id: &str) -> bool {
                feature_id.starts_with("N") || feature_id.starts_with("X")
            }

            // For each gene, find the worst RefSeq and ENSEMBL annotation and write one record.
            for (hgnc_id, anns) in anns_by_gene.iter() {
                // Get `HgncRecord` for the `hgnc_id` or skip this gene.
                let hgnc_record = match self.hgnc_map.as_ref().unwrap().get(hgnc_id) {
                    Some(hgnc_record) => hgnc_record,
                    None if hgnc_id.is_empty() => {
                        tracing::warn!(
                            "Empty HGNC id for {}:{}-{}, writing empty record.",
                            tsv_record.chromosome,
                            tsv_record.start,
                            tsv_record.end
                        );
                        &empty_hgnc_record
                    }
                    None => {
                        tracing::warn!(
                            "HGNC record for {} not found in HGNC map, writing empty record",
                            hgnc_id
                        );
                        &HgncRecord {
                            hgnc_id: hgnc_id.to_string(),
                            ensembl_gene_id: "".to_string(),
                            entrez_id: "".to_string(),
                            gene_symbol: "".to_string(),
                        }
                    }
                };

                tsv_record.clear_refseq_ensembl();

                // Since `anns` is sorted by severity, the first match we find is the worst one.
                let worst_refseq_ann = anns.iter().find(|ann| is_refseq_id(&ann.feature_id));
                let worst_ensembl_ann = anns.iter().find(|ann| is_ensembl_id(&ann.feature_id));

                let (refseq_source_ann, ensembl_source_ann) =
                    match (worst_refseq_ann, worst_ensembl_ann) {
                        (Some(refseq), Some(ensembl)) => (Some(refseq), Some(ensembl)),
                        (Some(refseq), None) => {
                            tracing::debug!(
                                "No Ensembl annotation for {}:{}-{}, using RefSeq annotation",
                                tsv_record.chromosome,
                                tsv_record.start,
                                tsv_record.end
                            );
                            (Some(refseq), Some(refseq))
                        }
                        (None, Some(ensembl)) => {
                            tracing::debug!(
                                "No RefSeq annotation for {}:{}-{}, using Ensembl annotation",
                                tsv_record.chromosome,
                                tsv_record.start,
                                tsv_record.end
                            );
                            (Some(ensembl), Some(ensembl))
                        }
                        (None, None) => {
                            tracing::debug!(
                                "No RefSeq or Ensembl annotation for {}:{}-{}.",
                                tsv_record.chromosome,
                                tsv_record.start,
                                tsv_record.end
                            );
                            (None, None)
                        }
                    };

                if let Some(ann) = ensembl_source_ann {
                    tsv_record.ensembl_gene_id = Some(hgnc_record.ensembl_gene_id.clone());
                    tsv_record.ensembl_transcript_id = Some(ann.feature_id.clone());
                    tsv_record.ensembl_transcript_coding =
                        Some(ann.feature_biotype.contains(&FeatureBiotype::Coding));
                    tsv_record.ensembl_hgvs_c.clone_from(&ann.hgvs_c);
                    tsv_record.ensembl_hgvs_p.clone_from(&ann.hgvs_p);
                    if !ann.consequences.is_empty() {
                        tsv_record.ensembl_effect = Some(
                            ann.consequences
                                .iter()
                                .map(|c| format!("\"{}\"", &c))
                                .collect::<Vec<_>>(),
                        );
                    }
                    tsv_record.ensembl_exon_dist = ann.distance;
                }

                if let Some(ann) = refseq_source_ann {
                    tsv_record.refseq_gene_id = Some(hgnc_record.entrez_id.clone());
                    tsv_record.refseq_transcript_id = Some(ann.feature_id.clone());
                    tsv_record.refseq_transcript_coding =
                        Some(ann.feature_biotype.contains(&FeatureBiotype::Coding));
                    tsv_record.refseq_hgvs_c.clone_from(&ann.hgvs_c);
                    tsv_record.refseq_hgvs_p.clone_from(&ann.hgvs_p);
                    if !ann.consequences.is_empty() {
                        tsv_record.refseq_effect = Some(
                            ann.consequences
                                .iter()
                                .map(|c| format!("\"{}\"", &c))
                                .collect::<Vec<_>>(),
                        );
                    }
                    tsv_record.refseq_exon_dist = ann.distance;
                }

                writeln!(self.inner, "{}", tsv_record.to_tsv().join("\t"))
                    .map_err(|e| anyhow::anyhow!("Error writing VarFish TSV record: {}", e))?;
            }

            Ok(())
        } else {
            // No annotations, write out without expanding `INFO/ANN` fields.
            writeln!(self.inner, "{}", tsv_record.to_tsv().join("\t"))
                .map_err(|e| anyhow::anyhow!("Error writing VarFish TSV record: {}", e))
        }
    }

    /// Fill `record` ClinVar fields.
    fn fill_clinvar(
        &self,
        record: &VcfRecord,
        tsv_record: &mut VarFishSeqvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        tsv_record.in_clinvar = record
            .info()
            .get("clinvar_germline_classification")
            .is_some();

        Ok(())
    }
}

/// A record, as written out to a VarFish TSV file.
#[derive(Debug, Default)]
pub struct VarFishSeqvarTsvRecord {
    pub release: String,
    pub chromosome: String,
    pub chromosome_no: u32,
    pub start: usize,
    pub end: usize,
    pub bin: u32,
    pub reference: String,
    pub alternative: String,
    pub var_type: String,

    // Writing out case_id and set_info is not used anyway.
    // pub case_id: String,
    // pub set_id: String,

    // The info field is not populated anyway.
    // pub info: String,
    pub genotype: String,

    pub num_hom_alt: u32,
    pub num_hom_ref: u32,
    pub num_het: u32,
    pub num_hemi_alt: u32,
    pub num_hemi_ref: u32,

    pub in_clinvar: bool,

    // NB: ExAc and 1000 Genomes are not written out anymore.
    // pub exac_frequency: String,
    // pub exac_homozygous: String,
    // pub exac_heterozygous: String,
    // pub exac_hemizygous: String,
    // pub thousand_genomes_frequency: String,
    // pub thousand_genomes_homozygous: String,
    // pub thousand_genomes_heterozygous: String,
    // pub thousand_genomes_hemizygous: String,
    pub gnomad_exomes_frequency: f32,
    pub gnomad_exomes_homozygous: i32,
    pub gnomad_exomes_heterozygous: i32,
    pub gnomad_exomes_hemizygous: i32,
    pub gnomad_genomes_frequency: f32,
    pub gnomad_genomes_homozygous: i32,
    pub gnomad_genomes_heterozygous: i32,
    pub gnomad_genomes_hemizygous: i32,

    pub refseq_gene_id: Option<String>,
    pub refseq_transcript_id: Option<String>,
    pub refseq_transcript_coding: Option<bool>,
    pub refseq_hgvs_c: Option<String>,
    pub refseq_hgvs_p: Option<String>,
    pub refseq_effect: Option<Vec<String>>,
    pub refseq_exon_dist: Option<i32>,

    pub ensembl_gene_id: Option<String>,
    pub ensembl_transcript_id: Option<String>,
    pub ensembl_transcript_coding: Option<bool>,
    pub ensembl_hgvs_c: Option<String>,
    pub ensembl_hgvs_p: Option<String>,
    pub ensembl_effect: Option<Vec<String>>,
    pub ensembl_exon_dist: Option<i32>,
}

impl VarFishSeqvarTsvRecord {
    /// Clear the `refseq_*` and `ensembl_*` fields.
    pub fn clear_refseq_ensembl(&mut self) {
        self.refseq_gene_id = None;
        self.refseq_transcript_id = None;
        self.refseq_transcript_coding = None;
        self.refseq_hgvs_c = None;
        self.refseq_hgvs_p = None;
        self.refseq_effect = None;
        self.refseq_exon_dist = None;

        self.ensembl_gene_id = None;
        self.ensembl_transcript_id = None;
        self.ensembl_transcript_coding = None;
        self.ensembl_hgvs_c = None;
        self.ensembl_hgvs_p = None;
        self.ensembl_effect = None;
        self.ensembl_exon_dist = None;
    }

    /// Convert to a `Vec<String>` suitable for writing to a VarFish TSV file.
    pub fn to_tsv(&self) -> Vec<String> {
        vec![
            self.release.clone(),
            self.chromosome.clone(),
            format!("{}", self.chromosome_no),
            format!("{}", self.start),
            format!("{}", self.end),
            format!("{}", self.bin),
            self.reference.clone(),
            self.alternative.clone(),
            self.var_type.clone(),
            String::from("."),
            String::from("."),
            String::from("{}"),
            self.genotype.clone(),
            format!("{}", self.num_hom_alt),
            format!("{}", self.num_hom_ref),
            format!("{}", self.num_het),
            format!("{}", self.num_hemi_alt),
            format!("{}", self.num_hemi_ref),
            if self.in_clinvar { "TRUE" } else { "FALSE" }.to_string(),
            // exac
            String::from("0"),
            String::from("0"),
            String::from("0"),
            String::from("0"),
            // thousand genomes
            String::from("0"),
            String::from("0"),
            String::from("0"),
            String::from("0"),
            format!("{}", self.gnomad_exomes_frequency),
            format!("{}", self.gnomad_exomes_homozygous),
            format!("{}", self.gnomad_exomes_heterozygous),
            format!("{}", self.gnomad_exomes_hemizygous),
            format!("{}", self.gnomad_genomes_frequency),
            format!("{}", self.gnomad_genomes_homozygous),
            format!("{}", self.gnomad_genomes_heterozygous),
            format!("{}", self.gnomad_genomes_hemizygous),
            self.refseq_gene_id.clone().unwrap_or(String::from(".")),
            self.refseq_transcript_id
                .clone()
                .unwrap_or(String::from(".")),
            self.refseq_transcript_coding
                .map(|refseq_transcript_coding| {
                    if refseq_transcript_coding {
                        String::from("TRUE")
                    } else {
                        String::from("FALSE")
                    }
                })
                .unwrap_or(String::from(".")),
            self.refseq_hgvs_c.clone().unwrap_or(String::from(".")),
            self.refseq_hgvs_p.clone().unwrap_or(String::from(".")),
            format!(
                "{{{}}}",
                self.refseq_effect
                    .as_ref()
                    .map(|refseq_effect| refseq_effect.join(","))
                    .unwrap_or_default()
            ),
            self.refseq_exon_dist
                .map(|refseq_exon_dist| format!("{}", refseq_exon_dist))
                .unwrap_or(String::from(".")),
            self.ensembl_gene_id.clone().unwrap_or(String::from(".")),
            self.ensembl_transcript_id
                .clone()
                .unwrap_or(String::from(".")),
            self.ensembl_transcript_coding
                .map(|ensembl_transcript_coding| {
                    if ensembl_transcript_coding {
                        String::from("TRUE")
                    } else {
                        String::from("FALSE")
                    }
                })
                .unwrap_or(String::from(".")),
            self.ensembl_hgvs_c.clone().unwrap_or(String::from(".")),
            self.ensembl_hgvs_p.clone().unwrap_or(String::from(".")),
            format!(
                "{{{}}}",
                self.ensembl_effect
                    .as_ref()
                    .map(|ensembl_effect| ensembl_effect.join(","))
                    .unwrap_or_default()
            ),
            self.ensembl_exon_dist
                .map(|ensembl_exon_dist| format!("{}", ensembl_exon_dist))
                .unwrap_or(String::from(".")),
        ]
    }
}

/// Implement `AnnotatedVcfWriter` for `VarFishTsvWriter`.
impl AsyncAnnotatedVariantWriter for VarFishSeqvarTsvWriter {
    async fn write_noodles_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error> {
        self.header = Some(header.clone());
        let header = &[
            "release",
            "chromosome",
            "chromosome_no",
            "start",
            "end",
            "bin",
            "reference",
            "alternative",
            "var_type",
            "case_id",
            "set_id",
            "info",
            "genotype",
            "num_hom_alt",
            "num_hom_ref",
            "num_het",
            "num_hemi_alt",
            "num_hemi_ref",
            "in_clinvar",
            "exac_frequency",
            "exac_homozygous",
            "exac_heterozygous",
            "exac_hemizygous",
            "thousand_genomes_frequency",
            "thousand_genomes_homozygous",
            "thousand_genomes_heterozygous",
            "thousand_genomes_hemizygous",
            "gnomad_exomes_frequency",
            "gnomad_exomes_homozygous",
            "gnomad_exomes_heterozygous",
            "gnomad_exomes_hemizygous",
            "gnomad_genomes_frequency",
            "gnomad_genomes_homozygous",
            "gnomad_genomes_heterozygous",
            "gnomad_genomes_hemizygous",
            "refseq_gene_id",
            "refseq_transcript_id",
            "refseq_transcript_coding",
            "refseq_hgvs_c",
            "refseq_hgvs_p",
            "refseq_effect",
            "refseq_exon_dist",
            "ensembl_gene_id",
            "ensembl_transcript_id",
            "ensembl_transcript_coding",
            "ensembl_hgvs_c",
            "ensembl_hgvs_p",
            "ensembl_effect",
            "ensembl_exon_dist",
        ];
        writeln!(self.inner, "{}", header.join("\t"))
            .map_err(|e| anyhow::anyhow!("Error writing VarFish TSV header: {}", e))
    }

    async fn write_noodles_record(
        &mut self,
        _header: &VcfHeader,
        record: &VcfRecord,
    ) -> Result<(), anyhow::Error> {
        let mut tsv_record = VarFishSeqvarTsvRecord::default();

        if !self.fill_coords(record, &mut tsv_record)? {
            // Record was not on canonical chromosome and should not be written out.
            return Ok(());
        }
        self.fill_genotype_and_freqs(record, &mut tsv_record)?;
        self.fill_bg_freqs(record, &mut tsv_record)?;
        self.fill_clinvar(record, &mut tsv_record)?;
        self.expand_refseq_ensembl_and_write(record, &mut tsv_record)
    }

    async fn shutdown(&mut self) -> Result<(), Error> {
        Ok(self.inner.flush()?)
    }

    fn set_hgnc_map(&mut self, hgnc_map: FxHashMap<String, HgncRecord>) {
        self.hgnc_map = Some(hgnc_map)
    }

    fn set_assembly(&mut self, assembly: Assembly) {
        self.assembly = Some(assembly)
    }

    fn set_pedigree(&mut self, pedigree: &PedigreeByName) {
        self.pedigree = Some(pedigree.clone())
    }
}

#[allow(clippy::large_enum_variant)]
pub(crate) enum AnnotatorEnum {
    Frequency(FrequencyAnnotator),
    Clinvar(ClinvarAnnotator),
    Consequence(ConsequenceAnnotator),
}

impl AnnotatorEnum {
    fn annotate(&self, var: &keys::Var, record: &mut VcfRecord) -> anyhow::Result<()> {
        match self {
            AnnotatorEnum::Frequency(a) => a.annotate(var, record),
            AnnotatorEnum::Clinvar(a) => a.annotate(var, record),
            AnnotatorEnum::Consequence(a) => a.annotate(var, record),
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

#[derive(Debug)]
pub struct FrequencyAnnotator {
    db: DBWithThreadMode<MultiThreaded>,
}
impl FrequencyAnnotator {
    pub fn new(db: DBWithThreadMode<MultiThreaded>) -> Self {
        Self { db }
    }

    pub(crate) fn from_path(path: impl AsRef<Path> + Display) -> anyhow::Result<Self> {
        // Open the frequency RocksDB database in read only mode.
        tracing::info!("Opening frequency database");
        tracing::debug!("RocksDB path = {}", &path);
        let options = rocksdb::Options::default();
        let db_freq = rocksdb::DB::open_cf_for_read_only(
            &options,
            &path,
            ["meta", "autosomal", "gonosomal", "mitochondrial"],
            false,
        )?;
        Ok(Self::new(db_freq))
    }

    /// Annotate record on autosomal chromosome with gnomAD exomes/genomes.
    pub fn annotate_record_auto(
        &self,
        key: &[u8],
        vcf_record: &mut noodles::vcf::variant::RecordBuf,
    ) -> Result<(), Error> {
        if let Some(freq) = self
            .db
            .get_cf(self.db.cf_handle("autosomal").as_ref().unwrap(), key)?
        {
            let auto_record = auto::Record::from_buf(&freq);

            vcf_record.info_mut().insert(
                "gnomad_exomes_an".into(),
                Some(field::Value::Integer(auto_record.gnomad_exomes.an as i32)),
            );
            vcf_record.info_mut().insert(
                "gnomad_exomes_hom".into(),
                Some(field::Value::Integer(
                    auto_record.gnomad_exomes.ac_hom as i32,
                )),
            );
            vcf_record.info_mut().insert(
                "gnomad_exomes_het".into(),
                Some(field::Value::Integer(
                    auto_record.gnomad_exomes.ac_het as i32,
                )),
            );

            vcf_record.info_mut().insert(
                "gnomad_genomes_an".into(),
                Some(field::Value::Integer(auto_record.gnomad_genomes.an as i32)),
            );
            vcf_record.info_mut().insert(
                "gnomad_genomes_hom".into(),
                Some(field::Value::Integer(
                    auto_record.gnomad_genomes.ac_hom as i32,
                )),
            );
            vcf_record.info_mut().insert(
                "gnomad_genomes_het".into(),
                Some(field::Value::Integer(
                    auto_record.gnomad_genomes.ac_het as i32,
                )),
            );
        };
        Ok(())
    }

    /// Annotate record on gonosomal chromosome with gnomAD exomes/genomes.
    pub fn annotate_record_xy(
        &self,
        key: &[u8],
        vcf_record: &mut noodles::vcf::variant::RecordBuf,
    ) -> Result<(), Error> {
        if let Some(freq) = self
            .db
            .get_cf(self.db.cf_handle("gonosomal").as_ref().unwrap(), key)?
        {
            let xy_record = xy::Record::from_buf(&freq);

            vcf_record.info_mut().insert(
                "gnomad_exomes_an".into(),
                Some(field::Value::Integer(xy_record.gnomad_exomes.an as i32)),
            );
            vcf_record.info_mut().insert(
                "gnomad_exomes_hom".into(),
                Some(field::Value::Integer(xy_record.gnomad_exomes.ac_hom as i32)),
            );
            vcf_record.info_mut().insert(
                "gnomad_exomes_het".into(),
                Some(field::Value::Integer(xy_record.gnomad_exomes.ac_het as i32)),
            );
            vcf_record.info_mut().insert(
                "gnomad_exomes_hemi".into(),
                Some(field::Value::Integer(
                    xy_record.gnomad_exomes.ac_hemi as i32,
                )),
            );

            vcf_record.info_mut().insert(
                "gnomad_genomes_an".into(),
                Some(field::Value::Integer(xy_record.gnomad_genomes.an as i32)),
            );
            vcf_record.info_mut().insert(
                "gnomad_genomes_hom".into(),
                Some(field::Value::Integer(
                    xy_record.gnomad_genomes.ac_hom as i32,
                )),
            );
            vcf_record.info_mut().insert(
                "gnomad_genomes_het".into(),
                Some(field::Value::Integer(
                    xy_record.gnomad_genomes.ac_het as i32,
                )),
            );
            vcf_record.info_mut().insert(
                "gnomad_genomes_hemi".into(),
                Some(field::Value::Integer(
                    xy_record.gnomad_genomes.ac_hemi as i32,
                )),
            );
        };
        Ok(())
    }

    /// Annotate record on mitochondrial genome with gnomAD mtDNA and HelixMtDb.
    pub fn annotate_record_mt(
        &self,
        key: &[u8],
        vcf_record: &mut noodles::vcf::variant::RecordBuf,
    ) -> Result<(), Error> {
        if let Some(freq) = self
            .db
            .get_cf(self.db.cf_handle("mitochondrial").as_ref().unwrap(), key)?
        {
            let mt_record = mt::Record::from_buf(&freq);

            vcf_record.info_mut().insert(
                "helix_an".into(),
                Some(field::Value::Integer(mt_record.helixmtdb.an as i32)),
            );
            vcf_record.info_mut().insert(
                "helix_hom".into(),
                Some(field::Value::Integer(mt_record.helixmtdb.ac_hom as i32)),
            );
            vcf_record.info_mut().insert(
                "helix_het".into(),
                Some(field::Value::Integer(mt_record.helixmtdb.ac_het as i32)),
            );

            vcf_record.info_mut().insert(
                "gnomad_genomes_an".into(),
                Some(field::Value::Integer(mt_record.gnomad_mtdna.an as i32)),
            );
            vcf_record.info_mut().insert(
                "gnomad_genomes_hom".into(),
                Some(field::Value::Integer(mt_record.gnomad_mtdna.ac_hom as i32)),
            );
            vcf_record.info_mut().insert(
                "gnomad_genomes_het".into(),
                Some(field::Value::Integer(mt_record.gnomad_mtdna.ac_het as i32)),
            );
        };
        Ok(())
    }

    fn annotate(&self, vcf_var: &keys::Var, record: &mut VcfRecord) -> anyhow::Result<()> {
        // Only attempt lookups into RocksDB for canonical contigs.
        if is_canonical(vcf_var.chrom.as_str()) {
            // Build key for RocksDB database from `vcf_var`.
            let key: Vec<u8> = vcf_var.clone().into();

            // Annotate with frequency.
            if CHROM_AUTO.contains(vcf_var.chrom.as_str()) {
                self.annotate_record_auto(&key, record)?;
            } else if CHROM_XY.contains(vcf_var.chrom.as_str()) {
                self.annotate_record_xy(&key, record)?;
            } else if CHROM_MT.contains(vcf_var.chrom.as_str()) {
                self.annotate_record_mt(&key, record)?;
            } else {
                tracing::trace!(
                    "Record @{:?} on non-canonical chromosome, skipping.",
                    &vcf_var
                );
            }
        }
        Ok(())
    }

    pub(crate) fn annotate_variant(
        &self,
        vcf_var: &VcfVariant,
    ) -> anyhow::Result<
        Option<crate::server::run::actix_server::seqvars_frequencies::FrequencyResultEntry>,
    > {
        // Only attempt lookups into RocksDB for canonical contigs.
        if !is_canonical(&vcf_var.chromosome) {
            return Ok(None);
        }

        // Build key for RocksDB database
        let vcf_var = keys::Var::from(
            &vcf_var.chromosome,
            vcf_var.position,
            &vcf_var.reference,
            &vcf_var.alternative,
        );
        let key: Vec<u8> = vcf_var.clone().into();
        use crate::server::run::actix_server::seqvars_frequencies::*;
        // Annotate with frequency.
        if CHROM_AUTO.contains(vcf_var.chrom.as_str()) {
            match self
                .db
                .get_cf(self.db.cf_handle("autosomal").as_ref().unwrap(), key)?
            {
                Some(freq) => {
                    let val = auto::Record::from_buf(&freq);
                    Ok(Some(FrequencyResultEntry::Autosomal(
                        AutosomalResultEntry {
                            gnomad_exomes_an: val.gnomad_exomes.an,
                            gnomad_exomes_hom: val.gnomad_exomes.ac_hom,
                            gnomad_exomes_het: val.gnomad_exomes.ac_het,
                            gnomad_genomes_an: val.gnomad_genomes.an,
                            gnomad_genomes_hom: val.gnomad_genomes.ac_hom,
                            gnomad_genomes_het: val.gnomad_genomes.ac_het,
                        },
                    )))
                }
                _ => Err(anyhow!("No frequency data found for variant {:?}", vcf_var)),
            }
        } else if CHROM_XY.contains(vcf_var.chrom.as_str()) {
            match self
                .db
                .get_cf(self.db.cf_handle("gonosomal").as_ref().unwrap(), key)?
            {
                Some(freq) => {
                    let val = xy::Record::from_buf(&freq);
                    Ok(Some(FrequencyResultEntry::Gonosomal(
                        GonosomalResultEntry {
                            gnomad_exomes_an: val.gnomad_exomes.an,
                            gnomad_exomes_hom: val.gnomad_exomes.ac_hom,
                            gnomad_exomes_het: val.gnomad_exomes.ac_het,
                            gnomad_exomes_hemi: val.gnomad_exomes.ac_hemi,
                            gnomad_genomes_an: val.gnomad_genomes.an,
                            gnomad_genomes_hom: val.gnomad_genomes.ac_hom,
                            gnomad_genomes_het: val.gnomad_genomes.ac_het,
                            gnomad_genomes_hemi: val.gnomad_genomes.ac_hemi,
                        },
                    )))
                }
                _ => Err(anyhow!("No frequency data found for variant {:?}", vcf_var)),
            }
        } else if CHROM_MT.contains(vcf_var.chrom.as_str()) {
            match self
                .db
                .get_cf(self.db.cf_handle("mitochondrial").as_ref().unwrap(), key)?
            {
                Some(freq) => {
                    let val = mt::Record::from_buf(&freq);
                    Ok(Some(FrequencyResultEntry::Mitochondrial(
                        MitochondrialResultEntry {
                            helix_an: val.helixmtdb.an,
                            helix_hom: val.helixmtdb.ac_hom,
                            helix_het: val.helixmtdb.ac_het,
                            gnomad_genomes_an: val.gnomad_mtdna.an,
                            gnomad_genomes_hom: val.gnomad_mtdna.ac_hom,
                            gnomad_genomes_het: val.gnomad_mtdna.ac_het,
                        },
                    )))
                }
                _ => Err(anyhow!("No frequency data found for variant {:?}", vcf_var)),
            }
        } else {
            tracing::trace!(
                "Record @{:?} on non-canonical chromosome, skipping.",
                &vcf_var
            );
            Ok(None)
        }
    }
}

#[derive(Debug)]
pub struct ClinvarAnnotator {
    db: DBWithThreadMode<rocksdb::MultiThreaded>,
}

impl ClinvarAnnotator {
    pub fn new(db: DBWithThreadMode<rocksdb::MultiThreaded>) -> Self {
        Self { db }
    }

    pub(crate) fn from_path(path: impl AsRef<Path> + Display) -> anyhow::Result<Self> {
        tracing::info!("Opening ClinVar database");
        tracing::debug!("RocksDB path = {}", &path);
        let options = rocksdb::Options::default();
        let db_clinvar =
            rocksdb::DB::open_cf_for_read_only(&options, &path, ["meta", "clinvar"], false)?;
        Ok(Self::new(db_clinvar))
    }

    /// Annotate record with ClinVar information
    pub fn annotate_record_clinvar(
        &self,
        key: &[u8],
        vcf_record: &mut noodles::vcf::variant::RecordBuf,
    ) -> Result<(), Error> {
        if let Some(raw_value) = self
            .db
            .get_cf(self.db.cf_handle("clinvar").as_ref().unwrap(), key)?
        {
            let record_list = annonars::pbs::clinvar::minimal::ExtractedVcvRecordList::decode(
                &mut Cursor::new(&raw_value),
            )?;

            let mut clinvar_vcvs = Vec::new();
            let mut clinvar_germline_classifications = Vec::new();
            for clinvar_record in record_list.records.iter() {
                let accession = clinvar_record.accession.as_ref().expect("must have VCV");
                let vcv = format!("{}.{}", accession.accession, accession.version);
                let classifications = clinvar_record
                    .classifications
                    .as_ref()
                    .expect("must have classifications");
                if let Some(germline_classification) = &classifications.germline_classification {
                    let description = germline_classification
                        .description
                        .as_ref()
                        .expect("description missing")
                        .to_string();
                    clinvar_vcvs.push(vcv);
                    clinvar_germline_classifications.push(description);
                }
            }

            vcf_record.info_mut().insert(
                "clinvar_vcv".into(),
                Some(field::Value::Array(field::value::Array::String(
                    clinvar_vcvs.into_iter().map(Some).collect::<Vec<_>>(),
                ))),
            );
            vcf_record.info_mut().insert(
                "clinvar_germline_classification".into(),
                Some(field::Value::Array(field::value::Array::String(
                    clinvar_germline_classifications
                        .into_iter()
                        .map(Some)
                        .collect::<Vec<_>>(),
                ))),
            );
        }

        Ok(())
    }

    fn annotate(&self, vcf_var: &keys::Var, record: &mut VcfRecord) -> anyhow::Result<()> {
        // Only attempt lookups into RocksDB for canonical contigs.
        if is_canonical(vcf_var.chrom.as_str()) {
            // Build key for RocksDB database from `vcf_var`.
            let key: Vec<u8> = vcf_var.clone().into();

            // Annotate with ClinVar information.
            self.annotate_record_clinvar(&key, record)?;
        }
        Ok(())
    }

    pub(crate) fn annotate_variant(
        &self,
        vcf_var: &VcfVariant,
    ) -> anyhow::Result<Option<crate::server::run::actix_server::seqvars_clinvar::ClinvarResultEntry>>
    {
        // Only attempt lookups into RocksDB for canonical contigs.
        if !is_canonical(&vcf_var.chromosome) {
            return Ok(None);
        }

        // Build key for RocksDB database
        let vcf_var = keys::Var::from(
            &vcf_var.chromosome,
            vcf_var.position,
            &vcf_var.reference,
            &vcf_var.alternative,
        );
        let key: Vec<u8> = vcf_var.clone().into();

        match self
            .db
            .get_cf(self.db.cf_handle("clinvar").as_ref().unwrap(), key)?
        {
            Some(raw_value) => {
                let record_list = annonars::pbs::clinvar::minimal::ExtractedVcvRecordList::decode(
                    &mut Cursor::new(&raw_value),
                )?;

                let mut clinvar_vcvs = Vec::new();
                let mut clinvar_germline_classifications = Vec::new();
                for clinvar_record in record_list.records.iter() {
                    let accession = clinvar_record.accession.as_ref().expect("must have VCV");
                    let vcv = format!("{}.{}", accession.accession, accession.version);
                    let classifications = clinvar_record
                        .classifications
                        .as_ref()
                        .expect("must have classifications");
                    if let Some(germline_classification) = &classifications.germline_classification
                    {
                        let description = germline_classification
                            .description
                            .as_ref()
                            .expect("description missing")
                            .to_string();
                        clinvar_vcvs.push(vcv);
                        clinvar_germline_classifications.push(description);
                    }
                }

                Ok(Some(
                    crate::server::run::actix_server::seqvars_clinvar::ClinvarResultEntry {
                        clinvar_vcv: clinvar_vcvs,
                        clinvar_germline_classification: clinvar_germline_classifications,
                    },
                ))
            }
            _ => Ok(None),
        }
    }
}

pub(crate) struct ConsequenceAnnotator {
    pub(crate) predictor: ConsequencePredictor,
}

impl ConsequenceAnnotator {
    fn new(predictor: ConsequencePredictor) -> Self {
        Self { predictor }
    }

    pub(crate) fn from_db_and_settings(
        tx_db: TxSeqDatabase,
        reference: Option<impl AsRef<Path>>,
        in_memory_reference: bool,
        transcript_settings: &TranscriptSettings,
    ) -> anyhow::Result<Self> {
        let args = transcript_settings;

        let provider = Arc::new(MehariProvider::new(
            tx_db,
            reference,
            in_memory_reference,
            MehariProviderConfigBuilder::default()
                .pick_transcript(args.pick_transcript.clone())
                .pick_transcript_mode(args.pick_transcript_mode)
                .build()?,
        ));
        let predictor = ConsequencePredictor::new(
            provider,
            ConsequencePredictorConfigBuilder::default()
                .report_most_severe_consequence_by(args.report_most_severe_consequence_by)
                .transcript_source(args.transcript_source)
                .keep_intergenic(args.keep_intergenic)
                .discard_utr_splice_variants(args.discard_utr_splice_variants)
                .vep_hgvs_shift(args.vep_hgvs_shift)
                .build()?,
        );
        Ok(Self::new(predictor))
    }

    fn annotate(&self, vcf_var: &keys::Var, record: &mut VcfRecord) -> anyhow::Result<()> {
        let keys::Var {
            chrom,
            pos,
            reference,
            alternative,
        } = vcf_var.clone();

        // Annotate with variant effect.
        if let Some(ann_fields) = self.predictor.predict(&VcfVariant {
            chromosome: chrom,
            position: pos,
            reference,
            alternative,
        })? {
            if !ann_fields.is_empty() {
                record.info_mut().insert(
                    "ANN".into(),
                    Some(field::Value::Array(field::value::Array::String(
                        ann_fields.iter().map(|ann| Some(ann.to_string())).collect(),
                    ))),
                );
            }
        }

        Ok(())
    }
}

impl Annotator {
    pub(crate) fn new(annotators: Vec<AnnotatorEnum>) -> Self {
        Self { annotators }
    }

    fn annotate(&self, vcf_record: &mut VcfRecord) -> anyhow::Result<()> {
        // Get first alternate allele record.
        let vcf_var = from_vcf_allele(vcf_record, 0);

        // Skip records with a deletion as alternative allele.
        if vcf_var.alternative == "*" {
            return Ok(());
        }

        for annotator in &self.annotators {
            annotator.annotate(&vcf_var, vcf_record)?;
        }

        Ok(())
    }
}

/// Main entry point for `annotate seqvars` sub command.
pub async fn run(_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("config = {:#?}", &args);
    if let Some(path_output_vcf) = &args.output.path_output_vcf {
        let mut writer = open_variant_writer(path_output_vcf).await?;
        run_with_writer(&mut writer, args).await?;
        writer.shutdown().await?;
    } else {
        // Load the HGNC xlink map.
        let hgnc_map = {
            tracing::info!("Loading HGNC map ...");
            let mut result = FxHashMap::default();

            let tsv_file = File::open(args.hgnc.as_ref().unwrap())?;
            let mut tsv_reader = csv::ReaderBuilder::new()
                .comment(Some(b'#'))
                .delimiter(b'\t')
                .from_reader(tsv_file);
            for record in tsv_reader.deserialize() {
                let record: HgncRecord = record?;
                result.insert(record.hgnc_id.clone(), record);
            }
            tracing::info!("... done loading HGNC map");

            result
        };

        let path_output_tsv = args
            .output
            .path_output_tsv
            .as_ref()
            .expect("tsv path must be set; vcf and tsv are mutually exclusive, vcf unset");
        let mut writer = VarFishSeqvarTsvWriter::with_path(path_output_tsv);

        // Load the pedigree.
        tracing::info!("Loading pedigree...");
        let pedigree = match &args.path_input_ped {
            Some(p) => {
                tracing::info!("Loading pedigree from file {}", p);
                PedigreeByName::from_path(p)?
            }
            None => {
                std::panic!("No pedigree file provided. This is required for tsv annotation.")
            }
        };
        writer.set_pedigree(&pedigree);
        tracing::info!("... done loading pedigree");

        writer.set_hgnc_map(hgnc_map);
        run_with_writer(&mut writer, args).await?;

        writer
            .flush()
            .map_err(|e| anyhow::anyhow!("problem flushing file: {}", e))?;
    }

    Ok(())
}

/// Run the annotation with the given `Write` within the `VcfWriter`.
async fn run_with_writer(
    writer: &mut impl AsyncAnnotatedVariantWriter,
    args: &Args,
) -> Result<(), Error> {
    tracing::info!("Open VCF and read header");
    let mut reader = open_variant_reader(&args.path_input_vcf).await?;

    let mut header_in = reader.read_header().await?;

    // Work around glnexus issue with RNC.
    if let Some(format) = header_in.formats_mut().get_mut("RNC") {
        *format.number_mut() = FormatNumber::Count(1);
        *format.type_mut() = FormatType::String;
    }

    // Guess genome release from contigs in VCF header.
    let genome_release = args.genome_release.map(|gr| match gr {
        GenomeRelease::Grch37 => Assembly::Grch37p10, // has chrMT!
        GenomeRelease::Grch38 => Assembly::Grch38,
    });
    let assembly = guess_assembly(&header_in, false, genome_release)?;
    writer.set_assembly(assembly);
    tracing::info!("Determined input assembly to be {:?}", &assembly);

    let annotator = setup_seqvars_annotator(
        &args.sources,
        args.reference.as_ref(),
        args.in_memory_reference,
        &args.transcript_settings,
        Some(assembly),
    )?;

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

    let header_out = build_header(
        &header_in,
        with_annotations,
        with_frequencies,
        with_clinvar,
        &additional_header_info,
    );

    // Perform the VCF annotation.
    tracing::info!("Annotating VCF ...");
    let start = Instant::now();
    let mut prev = Instant::now();
    let mut total_written = 0usize;

    writer.write_noodles_header(&header_out).await?;

    use futures::TryStreamExt;
    let mut records = reader.records(&header_in).await;
    loop {
        match records.try_next().await? {
            Some(mut vcf_record) => {
                // We currently can only process records with one alternate allele.
                if vcf_record.alternate_bases().len() != 1 {
                    tracing::error!(
                    "Found record with more than one alternate allele.  This is currently not supported. \
                    Please use `bcftools norm` to split multi-allelic records.  Record: {:?}",
                    &vcf_record
                );
                    anyhow::bail!("multi-allelic records not supported");
                }

                annotator.annotate(&mut vcf_record)?;

                if prev.elapsed().as_secs() >= 60 {
                    tracing::info!("at {:?}", from_vcf_allele(&vcf_record, 0));
                    prev = Instant::now();
                }

                // Write out the record.
                writer
                    .write_noodles_record(&header_out, &vcf_record)
                    .await?;
            }
            _ => {
                break; // all done
            }
        }

        total_written += 1;
        if let Some(max_var_count) = args.max_var_count {
            if total_written >= max_var_count {
                tracing::warn!(
                    "Stopping after {} records as requested by --max-var-count",
                    total_written
                );
                break;
            }
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

pub(crate) fn proto_assembly_from(assembly: &Assembly) -> Option<crate::pbs::txs::Assembly> {
    crate::pbs::txs::Assembly::from_str_name(&format!(
        "ASSEMBLY_{}",
        match assembly {
            Assembly::Grch38 => "GRCH38",
            _ => "GRCH37",
        }
    ))
}

pub(crate) fn setup_seqvars_annotator(
    sources: &Sources,
    reference: Option<impl AsRef<Path>>,
    in_memory_reference: bool,
    transcript_settings: &TranscriptSettings,
    assembly: Option<Assembly>,
) -> Result<Annotator, Error> {
    let mut annotators = vec![];

    // Add the frequency annotator if requested.
    if let Some(rocksdb_paths) = &sources.frequencies {
        let freq_dbs = initialize_frequency_annotators_for_assembly(rocksdb_paths, assembly)?;
        for freq_db in freq_dbs {
            annotators.push(AnnotatorEnum::Frequency(freq_db))
        }
    }

    // Add the ClinVar annotator if requested.
    if let Some(rocksdb_paths) = &sources.clinvar {
        let clinvar_dbs = initialize_clinvar_annotators_for_assembly(rocksdb_paths, assembly)?;
        for clinvar_db in clinvar_dbs {
            annotators.push(AnnotatorEnum::Clinvar(clinvar_db))
        }
    }

    // Add the consequence annotator if requested.
    if let Some(tx_sources) = &sources.transcripts {
        tracing::info!("Opening transcript database(s)");

        let databases = load_transcript_dbs_for_assembly(tx_sources, assembly)?;

        if databases.is_empty() {
            tracing::warn!("No suitable transcript databases found for requested assembly {:?}, therefore no consequence prediction will occur.", &assembly);
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
                    transcript_settings,
                )
                .map(AnnotatorEnum::Consequence)?,
            );
        }
    }

    let annotator = Annotator::new(annotators);
    Ok(annotator)
}

pub(crate) fn initialize_clinvar_annotators_for_assembly(
    rocksdb_paths: &[String],
    assembly: Option<Assembly>,
) -> Result<Vec<ClinvarAnnotator>, Error> {
    rocksdb_paths
        .iter()
        .filter_map(|rocksdb_path| {
            let skip = assembly.is_some_and(|a| !rocksdb_path.contains(path_component(a)));
            if !skip {
                tracing::info!(
                    "Loading ClinVar database for assembly {:?} from {}",
                    &assembly,
                    &rocksdb_path
                );
                Some(ClinvarAnnotator::from_path(rocksdb_path))
            } else {
                tracing::warn!(
                "Skipping ClinVar database as its assembly does not match the requested one ({:?})",
                &assembly
            );
                None
            }
        })
        .collect()
}

pub(crate) fn initialize_frequency_annotators_for_assembly(
    rocksdb_paths: &[String],
    assembly: Option<Assembly>,
) -> Result<Vec<FrequencyAnnotator>, Error> {
    rocksdb_paths.iter().filter_map(|rocksdb_path| {
        let skip = assembly.is_some_and(|a| !rocksdb_path.contains(path_component(a)));
        if !skip {
            tracing::info!(
                    "Loading frequency database for assembly {:?} from {}",
                    &assembly,
                    &rocksdb_path
                );
            Some(FrequencyAnnotator::from_path(rocksdb_path))
        } else {
            tracing::warn!("Skipping frequency database as its assembly does not match the requested one ({:?})", &assembly);
            None
        }
    }).collect()
}

pub(crate) fn load_transcript_dbs_for_assembly(
    tx_sources: &[String],
    assembly: Option<Assembly>,
) -> Result<Vec<TxSeqDatabase>, Error> {
    let pb_assembly = assembly.as_ref().and_then(proto_assembly_from);

    // Filter out any transcript databases that do not match the requested assembly.
    let check_assembly = |db: &TxSeqDatabase, assembly: crate::pbs::txs::Assembly| {
        db.source_version
            .iter()
            .map(|s| s.assembly)
            .any(|a| a == i32::from(assembly))
    };
    let databases = tx_sources
        .iter()
        .enumerate()
        .map(|(i, path)| (i, load_tx_db(path)))
        .filter_map(|(i, txdb)| match txdb {
            Ok(db) => match pb_assembly {
                Some(assembly) if check_assembly(&db, assembly) => Some(Ok(db)),
                Some(_) => {
                    tracing::info!("Skipping transcript database {} as its version {:?} does not support the requested assembly ({:?})", &tx_sources[i], &db.source_version, &assembly);
                    None
                }
                None => Some(Ok(db)),
            },
            Err(_) => Some(txdb),
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    Ok(databases)
}

/// Create for all alternate alleles from the given VCF record.
pub fn from_vcf_allele(value: &noodles::vcf::variant::RecordBuf, allele_no: usize) -> keys::Var {
    let chrom = value.reference_sequence_name().to_string();
    let pos: usize = value
        .variant_start()
        .expect("Telomeric breakends not supported")
        .get();
    let pos = i32::try_from(pos).unwrap();
    let reference = value.reference_bases().to_string();
    keys::Var {
        chrom,
        pos,
        reference,
        alternative: value.alternate_bases().as_ref()[allele_no].to_string(),
    }
}

#[cfg(test)]
mod test {
    use super::binning::bin_from_range;
    use super::{run, Args, PathOutput};
    use crate::annotate::cli::ConsequenceBy;
    use crate::annotate::cli::{Sources, TranscriptSettings};
    use clap_verbosity_flag::Verbosity;
    use pretty_assertions::assert_eq;
    use temp_testdir::TempDir;

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
            reference: None,
            in_memory_reference: true,
            genome_release: None,
            transcript_settings: TranscriptSettings {
                report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                ..Default::default()
            },
            path_input_vcf: String::from("tests/data/annotate/seqvars/brca1.examples.vcf"),
            output: PathOutput {
                path_output_vcf: Some(path_out.into_os_string().into_string().unwrap()),
                path_output_tsv: None,
            },
            max_var_count: None,
            path_input_ped: Some(String::from(
                "tests/data/annotate/seqvars/brca1.examples.ped",
            )),
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
            },
            hgnc: None,
        };

        run(&args_common, &args).await?;

        let actual = std::fs::read_to_string(args.output.path_output_vcf.unwrap())?;
        // remove vcf header lines starting with ##mehari
        let actual = actual
            .lines()
            .filter(|line| !line.starts_with("##mehari"))
            .collect::<Vec<_>>()
            .join("\n");
        insta::assert_snapshot!(actual);

        Ok(())
    }

    #[tokio::test]
    async fn smoke_test_output_tsv() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.tsv");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db";
        let assembly = "grch37";
        let args = Args {
            reference: None,
            in_memory_reference: true,
            genome_release: None,
            transcript_settings: TranscriptSettings {
                report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                ..Default::default()
            },
            path_input_vcf: String::from("tests/data/annotate/seqvars/brca1.examples.vcf"),
            output: PathOutput {
                path_output_vcf: None,
                path_output_tsv: Some(path_out.into_os_string().into_string().unwrap()),
            },
            max_var_count: None,
            path_input_ped: Some(String::from(
                "tests/data/annotate/seqvars/brca1.examples.ped",
            )),
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
            },
            hgnc: Some(format!("{prefix}/hgnc.tsv")),
        };

        run(&args_common, &args).await?;

        let actual = std::fs::read_to_string(args.output.path_output_tsv.unwrap())?;
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
        let path_out = temp.join("output.tsv");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db";
        let assembly = "grch37";
        let args = Args {
            reference: None,
            in_memory_reference: true,
            genome_release: None,
            transcript_settings: TranscriptSettings {
                report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                ..Default::default()
            },
            path_input_vcf: String::from("tests/data/annotate/seqvars/badly_formed_vcf_entry.vcf"),
            output: PathOutput {
                path_output_vcf: None,
                path_output_tsv: Some(path_out.into_os_string().into_string().unwrap()),
            },
            max_var_count: None,
            path_input_ped: Some(String::from(
                "tests/data/annotate/seqvars/badly_formed_vcf_entry.ped",
            )),
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
            },
            hgnc: Some(format!("{prefix}/hgnc.tsv")),
        };

        run(&args_common, &args).await?;

        let actual = std::fs::read_to_string(args.output.path_output_tsv.unwrap())?;
        let expected =
            std::fs::read_to_string("tests/data/annotate/seqvars/badly_formed_vcf_entry.tsv")?;
        assert_eq!(&expected, &actual);

        Ok(())
    }

    /// Mitochondnrial variants called by the DRAGEN v310 germline caller have special format
    /// considerations.
    ///
    /// See: https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/DRAGEN/MitochondrialCalling.htm
    #[tokio::test]
    async fn test_dragen_mitochondrial_variant() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.tsv");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db";
        let assembly = "grch37";
        let args = Args {
            reference: None,
            in_memory_reference: true,
            genome_release: None,
            transcript_settings: TranscriptSettings {
                report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                ..Default::default()
            },
            path_input_vcf: String::from("tests/data/annotate/seqvars/mitochondrial_variants.vcf"),
            output: PathOutput {
                path_output_vcf: None,
                path_output_tsv: Some(path_out.into_os_string().into_string().unwrap()),
            },
            max_var_count: None,
            path_input_ped: Some(String::from(
                "tests/data/annotate/seqvars/mitochondrial_variants.ped",
            )),
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
            },
            hgnc: Some(format!("{prefix}/hgnc.tsv")),
        };

        run(&args_common, &args).await?;

        let actual = std::fs::read_to_string(args.output.path_output_tsv.unwrap())?;
        let expected =
            std::fs::read_to_string("tests/data/annotate/seqvars/mitochondrial_variants.tsv")?;
        assert_eq!(&expected, &actual);

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
        let path_out = temp.join("output.tsv");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db";
        let assembly = "grch37";
        let args = Args {
            reference: None,
            in_memory_reference: true,
            genome_release: None,
            transcript_settings: TranscriptSettings {
                report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                ..Default::default()
            },
            path_input_vcf: String::from("tests/data/annotate/seqvars/clair3-glnexus-min.vcf"),
            output: PathOutput {
                path_output_vcf: None,
                path_output_tsv: Some(path_out.into_os_string().into_string().unwrap()),
            },
            max_var_count: None,
            path_input_ped: Some(String::from(
                "tests/data/annotate/seqvars/clair3-glnexus-min.ped",
            )),
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
            },
            hgnc: Some(format!("{prefix}/hgnc.tsv")),
        };

        run(&args_common, &args).await?;

        let actual = std::fs::read_to_string(args.output.path_output_tsv.unwrap())?;
        let expected =
            std::fs::read_to_string("tests/data/annotate/seqvars/clair3-glnexus-min.tsv")?;
        assert_eq!(&expected, &actual);

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
        let path_out = temp.join("output.tsv");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let prefix = "tests/data/annotate/db";
        let assembly = "grch38";
        let args = Args {
            reference: None,
            in_memory_reference: true,
            genome_release: None,
            transcript_settings: TranscriptSettings {
                report_most_severe_consequence_by: Some(ConsequenceBy::Gene),
                ..Default::default()
            },
            path_input_vcf: String::from("tests/data/annotate/seqvars/brca2_zar1l/brca2_zar1l.vcf"),
            output: PathOutput {
                path_output_vcf: None,
                path_output_tsv: Some(path_out.into_os_string().into_string().unwrap()),
            },
            max_var_count: None,
            path_input_ped: Some(String::from(
                "tests/data/annotate/seqvars/brca2_zar1l/brca2_zar1l.ped",
            )),
            sources: Sources {
                frequencies: Some(vec![format!("{prefix}/{assembly}/seqvars/freqs")]),
                clinvar: Some(vec![format!("{prefix}/{assembly}/seqvars/clinvar")]),
                transcripts: Some(vec![format!("{prefix}/{assembly}/txs.bin.zst")]),
            },
            hgnc: Some(format!("{prefix}/hgnc.tsv")),
        };

        run(&args_common, &args).await?;

        let actual = std::fs::read_to_string(args.output.path_output_tsv.unwrap())?;
        let header_actual = actual
            .lines()
            .filter(|l| l.starts_with('#'))
            .collect::<Vec<_>>();
        let records_actual = actual
            .lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<std::collections::HashSet<_>>();

        let expected =
            std::fs::read_to_string("tests/data/annotate/seqvars/brca2_zar1l/brca2_zar1l.tsv")?;
        let header_expected = expected
            .lines()
            .filter(|l| l.starts_with('#'))
            .collect::<Vec<_>>();
        let records_expected = expected
            .lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<std::collections::HashSet<_>>();

        assert_eq!(&header_actual, &header_expected);
        assert_eq!(&records_actual, &records_expected);

        Ok(())
    }
}
