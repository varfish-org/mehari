//! Annotation of structural variant VCF files.

use std::collections::HashSet;
use std::fs::OpenOptions;
use std::io::{BufReader, Write};
use std::path::Path;
use std::str::FromStr;
use std::{fs::File, io::BufWriter};

use crate::common::GenomeRelease;
use crate::db::create::seqvar_freqs::reading::{guess_assembly, CANONICAL};
use crate::ped::PedigreeByName;
use bio::data_structures::interval_tree::IntervalTree;
use chrono::Utc;
use clap::{Args as ClapArgs, Parser};
use flate2::write::GzEncoder;
use flate2::Compression;
use hgvs::static_data::Assembly;
use noodles::bgzf::Writer as BgzfWriter;
use noodles::vcf::record::alternate_bases::Allele;
use noodles::vcf::record::genotypes::Keys;
use noodles::vcf::record::{Genotypes, Position};
use noodles::vcf::Record as VcfRecord;
use noodles::vcf::{self, Header as VcfHeader};
use noodles::vcf::{
    header::info::key::Key as InfoKey, header::info::key::Other as InfoKeyOther,
    header::info::key::Standard as InfoKeyStandard,
    header::record::value::Other as HeaderValueOther,
    record::genotypes::genotype::field::value::Value as GenotypeValue,
    record::info::field::value::Value as InfoValue, Writer as VcfWriter,
};
use noodles_util::variant::reader::Builder as VariantReaderBuilder;
use rand::rngs::StdRng;
use rand::RngCore;
use rand_core::SeedableRng;
use serde::{Deserialize, Serialize};
use std::ops::Deref;
use strum::{Display, EnumIter, IntoEnumIterator};
use tempdir::TempDir;
use uuid::Uuid;

use self::bnd::Breakend;

use super::seqvars::binning::bin_from_range;
use super::seqvars::{binning, AnnotatedVcfWriter, CHROM_TO_CHROM_NO};

/// Command line arguments for `annotate strucvars` sub command.
#[derive(Parser, Debug, Clone)]
#[command(about = "Annotate structural variant VCF files", long_about = None)]
pub struct Args {
    /// Path to the mehari database folder.
    #[arg(long)]
    pub path_db: String,

    /// Genome release to use, default is to auto-detect.
    #[arg(long, value_enum)]
    pub genome_release: Option<GenomeRelease>,
    /// Path to the input PED file.
    #[arg(long)]
    pub path_input_ped: String,
    /// Path to the input VCF files.
    #[arg(long, required = true)]
    pub path_input_vcf: Vec<String>,
    #[command(flatten)]
    pub output: PathOutput,

    /// For debug purposes, maximal number of variants to annotate.
    #[arg(long)]
    pub max_var_count: Option<usize>,
    /// Maximal number of flatbuffers tables, should not need tweaking.  However, if you see a
    /// "too many tables" error output, increase this value.
    #[arg(long, default_value_t = 5_000_000)]
    pub max_fb_tables: usize,

    /// Seed for random number generator (UUIDs), if any.
    #[arg(long)]
    pub rng_seed: Option<u64>,
    /// Minimal reciprocal overlap to require.
    #[arg(long, default_value_t = 0.8)]
    pub min_overlap: f32,
    /// Slack to use around break-ends.
    #[arg(long, default_value_t = 50)]
    pub slack_bnd: i32,
    /// Slack to use around insertions.
    #[arg(long, default_value_t = 50)]
    pub slack_ins: i32,
}

/// Command line arguments to enforce either `--path-output-vcf` or `--path-output-tsv`.
#[derive(Debug, ClapArgs, Clone)]
#[group(required = true, multiple = false)]
pub struct PathOutput {
    /// Path to the output VCF file.
    #[arg(long)]
    pub path_output_vcf: Option<String>,

    /// Path to the output TSV file (for import into VarFish).
    #[arg(long)]
    pub path_output_tsv: Option<String>,
}

pub mod keys {
    use std::str::FromStr;

    use noodles::vcf::{
        header::format::key::Key as FormatKey, header::format::key::Other as FormatKeyOther,
        header::info::key::Key as InfoKey, header::info::key::Other as InfoKeyOther,
    };

    lazy_static::lazy_static! {
        /// File format.
        pub static ref FILE_FORMAT: String = String::from("VCFv4.3");

        /// Key for `PASS`.
        pub static ref PASS: String = String::from("PASS");

        /// Key for `INFO/PE_ORIENTATION`.
        pub static ref PE_ORIENTATION: InfoKey = InfoKey::Other(InfoKeyOther::from_str("PE_ORIENTATION").unwrap());
        /// Key for `INFO/CIPOS`.
        pub static ref CIPOS: InfoKey = InfoKey::Other(InfoKeyOther::from_str("CIPOS").unwrap());
        /// Key for `INFO/CIEND`.
        pub static ref CIEND: InfoKey = InfoKey::Other(InfoKeyOther::from_str("CIEND").unwrap());
        /// Key for `INFO/CALLERS`.
        pub static ref CALLERS: InfoKey = InfoKey::Other(InfoKeyOther::from_str("CALLERS").unwrap());
        /// Key for `INFO/SVTYPE`.
        pub static ref SV_TYPE: InfoKey = InfoKey::Other(InfoKeyOther::from_str("SVTYPE").unwrap());
        /// Key for `INFO/SVSUBTYPE`.
        pub static ref SV_SUB_TYPE: InfoKey = InfoKey::Other(InfoKeyOther::from_str("SVSUBTYPE").unwrap());

        /// Key for `FORMAT/GT`.
        pub static ref GT: FormatKey = FormatKey::Other(FormatKeyOther::from_str("GT").unwrap());
        /// Key for `FORMAT/FT`.
        pub static ref FT: FormatKey = FormatKey::Other(FormatKeyOther::from_str("FT").unwrap());
        /// Key for `FORMAT/GQ`.
        pub static ref GQ: FormatKey = FormatKey::Other(FormatKeyOther::from_str("GQ").unwrap());
        /// Key for `FORMAT/PEC`.
        pub static ref PEC: FormatKey = FormatKey::Other(FormatKeyOther::from_str("PEC").unwrap());
        /// Key for `FORMAT/PEV`.
        pub static ref PEV: FormatKey = FormatKey::Other(FormatKeyOther::from_str("PEV").unwrap());
        /// Key for `FORMAT/SRC`.
        pub static ref SRC: FormatKey = FormatKey::Other(FormatKeyOther::from_str("SRC").unwrap());
        /// Key for `FORMAT/SRV`.
        pub static ref SRV: FormatKey = FormatKey::Other(FormatKeyOther::from_str("SRV").unwrap());
        /// Key for `FORMAT/AMQ`.
        pub static ref AMQ: FormatKey = FormatKey::Other(FormatKeyOther::from_str("AMQ").unwrap());
        /// Key for `FORMAT/CN`.
        pub static ref CN: FormatKey = FormatKey::Other(FormatKeyOther::from_str("CN").unwrap());
        /// Key for `FORMAT/ANC`.
        pub static ref ANC: FormatKey = FormatKey::Other(FormatKeyOther::from_str("ANC").unwrap());
        /// Key for `FORMAT/PC`.
        pub static ref PC: FormatKey = FormatKey::Other(FormatKeyOther::from_str("PC").unwrap());
    }
}

/// Code for building the VCF header to be written out.
pub mod vcf_header {
    use std::str::FromStr;

    use chrono::NaiveDate;
    use hgvs::static_data::{Assembly, ASSEMBLY_INFOS};
    use noodles::vcf::header::record::key::Key as HeaderKey;
    use noodles::vcf::header::record::value::map::{Contig, Filter, Format, Info, Meta, Other};
    use noodles::vcf::header::record::value::Map;
    use noodles::vcf::header::{self, record, Number};
    use noodles::vcf::{
        header::{Builder, FileFormat},
        Header,
    };

    use crate::db::create::seqvar_freqs::reading::CANONICAL;
    use crate::ped::{Disease, PedigreeByName, Sex};

    /// Major VCF version to use.
    static FILE_FORMAT_MAJOR: u32 = 4;
    /// Minor VCF version to use.
    static FILE_FORMAT_MINOR: u32 = 3;
    /// The string to write out as the source.
    static SOURCE: &str = "mehari";

    /// Construct VCF header.
    ///
    /// # Arguments
    ///
    /// * `assembly` - Genome assembly to use.  The canonical contigs will be taken from here.
    /// * `pedigree` - Pedigree to use.  Will write out appropriate `META`, `SAMPLE`, and
    ///    `PEDIGREE` header lines.
    /// * `date` - Date to use for the `fileDate` header line.
    ///
    /// # Returns
    ///
    /// The constructed VCF header.
    pub fn build(
        assembly: Assembly,
        pedigree: &PedigreeByName,
        date: &NaiveDate,
    ) -> Result<Header, anyhow::Error> {
        let builder = add_meta_leading(Header::builder(), date)?;
        let builder = add_meta_contigs(builder, assembly)?;
        let builder = add_meta_alt(builder)?;
        let builder = add_meta_info(builder)?;
        let builder = add_meta_filter(builder)?;
        let builder = add_meta_format(builder)?;
        let builder = add_meta_pedigree(builder, pedigree)?;

        Ok(builder.build())
    }

    /// Add the leading header lines, such as `fileformat`, `fileDate`, `source`, etc.
    fn add_meta_leading(builder: Builder, date: &NaiveDate) -> Result<Builder, anyhow::Error> {
        Ok(builder
            .set_file_format(FileFormat::new(FILE_FORMAT_MAJOR, FILE_FORMAT_MINOR))
            .insert(
                HeaderKey::other("fileDate").unwrap(),
                record::value::Other::from(date.format("%Y%m%d").to_string()),
            )
            .insert(
                HeaderKey::other("source").unwrap(),
                record::value::Other::from(SOURCE),
            ))
    }

    /// Add the `contig` header lines.
    ///
    /// NB: `Assembly::Grch37` does not contain chrMT, but `Assembly::Grch37p10` does.
    fn add_meta_contigs(builder: Builder, assembly: Assembly) -> Result<Builder, anyhow::Error> {
        let mut builder = builder;
        let assembly_info = &ASSEMBLY_INFOS[assembly];
        let assembly_name = match assembly {
            Assembly::Grch37 | Assembly::Grch37p10 => String::from("GRCh37"),
            Assembly::Grch38 => String::from("GRCh38"),
        };

        for sequence in &assembly_info.sequences {
            if CANONICAL.contains(&sequence.name.as_ref()) {
                let mut contig = Map::<Contig>::try_from(vec![
                    (String::from("length"), format!("{}", sequence.length)),
                    (String::from("assembly"), assembly_name.clone()),
                    (String::from("accession"), sequence.refseq_ac.clone()),
                ])?;

                Map::<Contig>::new();
                *contig.length_mut() = Some(sequence.length);
                builder = builder.add_contig(sequence.name.parse()?, contig);
            }
        }

        Ok(builder)
    }

    /// Add the `ALT` header lines.
    fn add_meta_alt(builder: Builder) -> Result<Builder, anyhow::Error> {
        use noodles::vcf::{
            header::record::value::map::AlternativeAllele,
            record::alternate_bases::allele::{
                symbol::{structural_variant::Type, StructuralVariant},
                Symbol,
            },
        };

        let del_id = Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion));
        let del_alt = Map::<AlternativeAllele>::new("Deletion");
        let del_me_id = Symbol::StructuralVariant(StructuralVariant::new(
            Type::Deletion,
            vec![String::from("ME")],
        ));
        let del_me_alt = Map::<AlternativeAllele>::new("Deletion of mobile element");
        let ins_id = Symbol::StructuralVariant(StructuralVariant::from(Type::Insertion));
        let ins_alt = Map::<AlternativeAllele>::new("Insertion");
        let ins_me_id = Symbol::StructuralVariant(StructuralVariant::new(
            Type::Insertion,
            vec![String::from("ME")],
        ));
        let ins_me_alt = Map::<AlternativeAllele>::new("Insertion of mobile element");
        let dup_id = Symbol::StructuralVariant(StructuralVariant::from(Type::Duplication));
        let dup_alt = Map::<AlternativeAllele>::new("Duplication");
        let dup_tnd_id = Symbol::StructuralVariant(StructuralVariant::new(
            Type::Duplication,
            vec![String::from("TANDEM")],
        ));
        let dup_tnd_alt = Map::<AlternativeAllele>::new("Tandem Duplication");
        let dup_dsp_id = Symbol::StructuralVariant(StructuralVariant::new(
            Type::Duplication,
            vec![String::from("DISPERSED")],
        ));
        let dup_dsp_alt = Map::<AlternativeAllele>::new("Dispersed Duplication");
        let cnv_id = Symbol::StructuralVariant(StructuralVariant::from(Type::CopyNumberVariation));
        let cnv_alt = Map::<AlternativeAllele>::new("Copy number variation");
        let bnd_id = Symbol::StructuralVariant(StructuralVariant::from(Type::Breakend));
        let bnd_alt = Map::<AlternativeAllele>::new("Breakend");

        Ok(builder
            .add_alternative_allele(del_id, del_alt)
            .add_alternative_allele(del_me_id, del_me_alt)
            .add_alternative_allele(ins_id, ins_alt)
            .add_alternative_allele(ins_me_id, ins_me_alt)
            .add_alternative_allele(dup_id, dup_alt)
            .add_alternative_allele(dup_tnd_id, dup_tnd_alt)
            .add_alternative_allele(dup_dsp_id, dup_dsp_alt)
            .add_alternative_allele(cnv_id, cnv_alt)
            .add_alternative_allele(bnd_id, bnd_alt))
    }

    /// Add the `INFO` header lines.
    fn add_meta_info(builder: Builder) -> Result<Builder, anyhow::Error> {
        use header::info::key::*;
        use header::record::value::map::info::Type;

        Ok(builder
            .add_info(END_POSITION, Map::<Info>::from(&END_POSITION))
            .add_info(
                POSITION_CONFIDENCE_INTERVALS,
                Map::<Info>::from(&POSITION_CONFIDENCE_INTERVALS),
            )
            .add_info(
                END_CONFIDENCE_INTERVALS,
                Map::<Info>::from(&END_CONFIDENCE_INTERVALS),
            )
            .add_info(
                Key::from_str("callers")?,
                Map::<Info>::new(
                    Number::Unknown,
                    Type::String,
                    "Callers that detected the variant",
                ),
            )
            // The SV UUID will only be written out temporarily until we don't need TSV anymore.
            // TODO: remove this once we don't need TSV anymore.
            .add_info(
                Key::from_str("sv_uuid")?,
                Map::<Info>::new(
                    Number::Unknown,
                    Type::String,
                    "Temporary UUID; needed for TSV output",
                ),
            )
            // Note that we will write out the sub type here, actually.
            .add_info(SV_TYPE, Map::<Info>::from(&SV_TYPE)))
    }

    /// Add the `FILTER` header lines.
    fn add_meta_filter(builder: Builder) -> Result<Builder, anyhow::Error> {
        Ok(builder.add_filter("PASS", Map::<Filter>::new("All filters passed")))
    }

    /// Add the `FORMAT` header lines.
    fn add_meta_format(builder: Builder) -> Result<Builder, anyhow::Error> {
        use header::format::key::*;
        use header::record::value::map::format::Type;

        Ok(builder
            .add_format(GENOTYPE, Map::<Format>::from(&GENOTYPE))
            .add_format(FILTER, Map::<Format>::from(&FILTER))
            .add_format(
                CONDITIONAL_GENOTYPE_QUALITY,
                Map::<Format>::from(&CONDITIONAL_GENOTYPE_QUALITY),
            )
            .add_format(
                Key::from_str("pec")?,
                Map::<Format>::new(Number::Count(1), Type::Integer, "Paired-end coverage"),
            )
            .add_format(
                Key::from_str("pev")?,
                Map::<Format>::new(
                    Number::Count(1),
                    Type::Integer,
                    "Paired-end variant support",
                ),
            )
            .add_format(
                Key::from_str("src")?,
                Map::<Format>::new(Number::Count(1), Type::Integer, "Split-end coverage"),
            )
            .add_format(
                Key::from_str("src")?,
                Map::<Format>::new(Number::Count(1), Type::Integer, "Split-end variant support"),
            )
            .add_format(
                Key::from_str("amq")?,
                Map::<Format>::new(Number::Count(1), Type::Integer, "Average mapping quality"),
            )
            .add_format(
                GENOTYPE_COPY_NUMBER,
                Map::<Format>::from(&GENOTYPE_COPY_NUMBER),
            )
            .add_format(
                Key::from_str("anc")?,
                Map::<Format>::new(
                    Number::Count(1),
                    Type::Integer,
                    "Average normalied coverage",
                ),
            )
            .add_format(
                Key::from_str("pc")?,
                Map::<Format>::new(
                    Number::Count(1),
                    Type::Integer,
                    "Point count (windows/targets/probes)",
                ),
            ))
    }

    /// Helper that returns header string value for `sex`.
    fn sex_str(sex: Sex) -> String {
        match sex {
            Sex::Male => String::from("Male"),
            Sex::Female => String::from("Female"),
            Sex::Unknown => String::from("Unknown"),
        }
    }

    // Helper that returns header string value for `disease`.
    fn disease_str(disease: Disease) -> String {
        match disease {
            Disease::Affected => String::from("Affected"),
            Disease::Unaffected => String::from("Unaffected"),
            Disease::Unknown => String::from("Unknown"),
        }
    }

    /// Add the `PEDIGREE` and supporting `SAMPLE` and `META` lines; set sample names.
    fn add_meta_pedigree(
        builder: Builder,
        pedigree: &PedigreeByName,
    ) -> Result<Builder, anyhow::Error> {
        let pedigree_key =
            header::record::key::Key::other("PEDIGREE").expect("invalid other meta key");
        let sample_key = header::record::key::Key::other("SAMPLE").expect("invalid other meta key");

        let mut builder = add_meta_fields(builder);

        // Wait for https://github.com/zaeleus/noodles/issues/162#issuecomment-1514444101
        // let mut b: record::value::map::Builder<record::value::map::Other> = Map::<noodles::vcf::header::record::value::map::Other>::builder();

        for i in pedigree.individuals.values() {
            builder = builder.add_sample_name(i.name.clone());

            // Add SAMPLE entry.
            {
                let mut entries = vec![(String::from("ID"), i.name.clone())];
                entries.push((String::from("Sex"), sex_str(i.sex)));
                entries.push((String::from("Disease"), disease_str(i.disease)));

                builder = builder.insert(
                    sample_key.clone(),
                    header::record::value::Other::Map(
                        i.name.clone(),
                        Map::<Other>::try_from(entries)?,
                    ),
                );
            }

            // Add PEDIGREE entry.
            {
                let mut entries = vec![(String::from("ID"), i.name.clone())];
                if let Some(father) = i.father.as_ref() {
                    entries.push((String::from("Father"), father.clone()));
                }
                if let Some(mother) = i.mother.as_ref() {
                    entries.push((String::from("Mother"), mother.clone()));
                }

                builder = builder.insert(
                    pedigree_key.clone(),
                    header::record::value::Other::Map(
                        i.name.clone(),
                        Map::<Other>::try_from(entries)?,
                    ),
                );
            }
        }

        Ok(builder)
    }

    // Define fields for the gonosomal karyotype, (sex for canonical), affected status
    fn add_meta_fields(builder: Builder) -> Builder {
        builder
            .add_meta(
                "Sex",
                Map::<Meta>::new(vec![
                    String::from("Male"),
                    String::from("Female"),
                    String::from("Other"),
                    String::from("Unknown"),
                ]),
            )
            .add_meta(
                "GonosomalKaryotype",
                Map::<Meta>::new(vec![
                    // "canonical" female
                    String::from("XX"),
                    // "canonical" male
                    String::from("XY"),
                    // Turner syndrome
                    String::from("XO"),
                    // Klinefelter syndrome
                    String::from("XXY"),
                    // Triple X syndrome
                    String::from("XXX"),
                    // Jacobs syndrome
                    String::from("XYY"),
                    // Other
                    String::from("Other"),
                    // Unknown
                    String::from("Unknown"),
                ]),
            )
            .add_meta(
                "Affected",
                Map::<Meta>::new(vec![
                    String::from("Yes"),
                    String::from("No"),
                    String::from("Unknown"),
                ]),
            )
    }
}

/// Writing of structural variants to VarFish TSV files.
struct VarFishStrucvarTsvWriter {
    /// The actual (compressed) text output writer.
    inner: Box<dyn Write>,
    /// Assembly information.
    assembly: Option<Assembly>,
    /// Pedigree information.
    pedigree: Option<PedigreeByName>,
    /// VCF Header for equivalent output file.
    header: Option<VcfHeader>,
}

/// Per-genotype call information.
#[derive(Debug, Default, Serialize, Deserialize, PartialEq, Clone)]
pub struct GenotypeInfo {
    /// Sample name.
    pub name: String,
    /// Genotype value.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gt: Option<String>,
    /// Per-genotype filter values.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ft: Option<Vec<String>>,
    /// Genotype quality.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gq: Option<i32>,
    /// Paired-end coverage.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pec: Option<i32>,
    /// Paired-end variant support.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pev: Option<i32>,
    /// Split-read coverage.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub src: Option<i32>,
    /// Split-read variant support.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub srv: Option<i32>,
    /// Average mapping quality.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub amq: Option<i32>,
    /// Copy number.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cn: Option<i32>,
    /// Average normalized coverage.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anc: Option<f32>,
    /// Point count (windows/targets/probes).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pc: Option<i32>,
}

#[derive(Debug, Default, PartialEq, Serialize, Deserialize, Clone)]
pub struct GenotypeCalls {
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

            if let Some(ft) = &entry.ft {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str("\"\"\"ft\"\"\":[");
                let mut first_ft = true;
                for ft_entry in ft {
                    if first_ft {
                        first_ft = false;
                    } else {
                        result.push(',');
                    }
                    result.push_str(&format!("\"\"\"{}\"\"\"", ft_entry));
                }
                result.push(']');
            }

            if let Some(gq) = &entry.gq {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"ad\"\"\":{}", gq));
            }

            if let Some(pec) = &entry.pec {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"pec\"\"\":{}", pec));
            }

            if let Some(pev) = &entry.pev {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"pev\"\"\":{}", pev));
            }

            if let Some(src) = &entry.src {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"src\"\"\":{}", src));
            }

            if let Some(srv) = &entry.srv {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"srv\"\"\":{}", srv));
            }

            if let Some(amq) = &entry.amq {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"amq\"\"\":{}", amq));
            }

            if let Some(cn) = &entry.cn {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"cn\"\"\":{}", cn));
            }

            if let Some(anc) = &entry.anc {
                if prev {
                    result.push(',');
                }
                prev = true;
                result.push_str(&format!("\"\"\"anc\"\"\":{}", anc));
            }

            if let Some(pc) = &entry.pc {
                if prev {
                    result.push(',');
                }
                // prev = true;
                result.push_str(&format!("\"\"\"pc\"\"\":{}", pc));
            }

            result.push('}');
        }

        result.push('}');
        result
    }
}

/// Implement `AnnotatedVcfWriter` for `VarFishTsvWriter`.
impl AnnotatedVcfWriter for VarFishStrucvarTsvWriter {
    fn write_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error> {
        self.header = Some(header.clone());
        let header = &[
            "release",
            "chromosome",
            "chromosome_no",
            "bin",
            "chromosome2",
            "chromosome_no2",
            "bin2",
            "pe_orientation",
            "start",
            "end",
            "start_ci_left",
            "start_ci_right",
            "end_ci_left",
            "end_ci_right",
            "case_id",
            "set_id",
            "sv_uuid",
            "callers",
            "sv_type",
            "sv_sub_type",
            "info",
            "genotype",
        ];
        writeln!(self.inner, "{}", header.join("\t"))
            .map_err(|e| anyhow::anyhow!("Error writing VarFish TSV header: {}", e))
    }

    fn write_record(&mut self, record: &VcfRecord) -> Result<(), anyhow::Error> {
        let mut tsv_record = VarFishStrucvarTsvRecord::default();

        match self.assembly {
            Some(Assembly::Grch37) | Some(Assembly::Grch37p10) => {
                tsv_record.release = String::from("GRCh37");
            }
            Some(Assembly::Grch38) => {
                tsv_record.release = String::from("GRCh38");
            }
            _ => panic!("assembly must have been set"),
        }

        tsv_record.chromosome = record.chromosome().to_string();
        tsv_record.chromosome_no = *CHROM_TO_CHROM_NO
            .get(&tsv_record.chromosome)
            .expect("chromosome not canonical");

        tsv_record.chromosome2 = record.chromosome().to_string();
        tsv_record.chromosome_no2 = *CHROM_TO_CHROM_NO
            .get(&tsv_record.chromosome2)
            .expect("chromosome not canonical");

        tsv_record.start = {
            let start: usize = record.position().into();
            start as i32
        };
        tsv_record.end = {
            let pos_end = record
                .info()
                .get(&noodles::vcf::header::info::key::END_POSITION);
            if let Some(Some(noodles::vcf::record::info::field::Value::Integer(pos_end))) = pos_end
            {
                *pos_end
            } else {
                // E.g., if INS
                tsv_record.start
            }
        };

        let sv_uuid = record
            .info()
            .get(&noodles::vcf::header::info::key::Key::from_str("sv_uuid")?);
        if let Some(Some(noodles::vcf::record::info::field::Value::String(sv_uuid))) = sv_uuid {
            tsv_record.sv_uuid = Uuid::from_str(sv_uuid)?;
        }
        let callers = record
            .info()
            .get(&noodles::vcf::header::info::key::Key::from_str("callers")?);
        if let Some(Some(noodles::vcf::record::info::field::Value::String(callers))) = callers {
            tsv_record.callers = callers.split(',').map(|x| x.to_string()).collect();
        }
        let sv_sub_type = record.info().get(&noodles::vcf::header::info::key::SV_TYPE);
        if let Some(Some(noodles::vcf::record::info::field::Value::String(sv_sub_type))) =
            sv_sub_type
        {
            tsv_record.sv_type =
                SvType::from_str(sv_sub_type.split(':').next().expect("invalid INFO/SVTYPE"))?;
            tsv_record.sv_sub_type = SvSubType::from_str(sv_sub_type)?;
        }

        tsv_record.bin = {
            if tsv_record.sv_type != SvType::Bnd {
                bin_from_range(tsv_record.start - 1, tsv_record.end)? as u32
            } else {
                bin_from_range(tsv_record.start - 1, tsv_record.start)? as u32
            }
        };
        tsv_record.bin2 = {
            if tsv_record.sv_type == SvType::Bnd {
                bin_from_range(tsv_record.end - 1, tsv_record.end)? as u32
            } else {
                tsv_record.bin
            }
        };

        tsv_record.pe_orientation = if tsv_record.sv_type == SvType::Bnd {
            let alt = record.alternate_bases().deref()[0].to_string();
            bnd::Breakend::from_ref_alt_str("N", alt.as_ref())?.pe_orientation
        } else {
            tsv_record.sv_type.into()
        };

        // Fill `tsv_record.genotypes`.
        let individuals = &self
            .pedigree
            .as_ref()
            .expect("pedigree must have been set")
            .individuals;
        // First, create genotype info records.
        let mut gt_it = record.genotypes().deref().iter();
        for (_, indiv) in individuals {
            tsv_record.genotype.entries.push(GenotypeInfo {
                name: indiv.name.clone(),
                ..Default::default()
            });

            let mut entry = tsv_record.genotype.entries.last_mut().expect("just pushed");
            let gt = gt_it.next().expect("genotype iterator exhausted");

            for (key, value) in gt.deref().iter() {
                match (key.as_ref(), value) {
                    ("GT", Some(GenotypeValue::String(gt))) => {
                        entry.gt = Some(gt.clone());
                    }
                    ("FT", Some(GenotypeValue::String(ft))) => {
                        entry.ft = Some(ft.split(';').map(|s| s.to_string()).collect());
                    }
                    ("GQ", Some(GenotypeValue::Integer(gq))) => {
                        entry.gq = Some(*gq);
                    }
                    // pec
                    ("pec", Some(GenotypeValue::Integer(pec))) => {
                        entry.pec = Some(*pec);
                    }
                    // pev
                    ("pev", Some(GenotypeValue::Integer(pev))) => {
                        entry.pev = Some(*pev);
                    }
                    // src
                    ("src", Some(GenotypeValue::Integer(src))) => {
                        entry.src = Some(*src);
                    }
                    // amq
                    ("CN", Some(GenotypeValue::Integer(cn))) => {
                        entry.cn = Some(*cn);
                    }
                    // anc
                    ("anc", Some(GenotypeValue::Float(anc))) => {
                        entry.anc = Some(*anc);
                    }
                    // pc
                    ("pc", Some(GenotypeValue::Integer(pc))) => {
                        entry.pc = Some(*pc);
                    }
                    // Ignore all other keys.
                    _ => (),
                }
            }
        }

        writeln!(
            self.inner,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{{}}\t{}",
            &tsv_record.release,
            &tsv_record.chromosome,
            &tsv_record.chromosome_no,
            &tsv_record.bin,
            &tsv_record.chromosome2,
            &tsv_record.chromosome_no2,
            &tsv_record.bin2,
            &tsv_record.pe_orientation,
            &tsv_record.start,
            &tsv_record.end,
            &tsv_record.start_ci_left,
            &tsv_record.start_ci_right,
            &tsv_record.end_ci_left,
            &tsv_record.end_ci_right,
            &tsv_record.case_id,
            &tsv_record.set_id,
            &tsv_record.sv_uuid,
            tsv_record.callers.join(";"),
            &tsv_record.sv_type,
            &tsv_record.sv_sub_type,
            &tsv_record.genotype.for_tsv(),
        ).map_err(|e| anyhow::anyhow!("Error writing VarFish TSV record: {}", e))
    }

    fn set_assembly(&mut self, assembly: Assembly) {
        self.assembly = Some(assembly)
    }

    fn set_pedigree(&mut self, pedigree: &PedigreeByName) {
        self.pedigree = Some(pedigree.clone())
    }
}

impl VarFishStrucvarTsvWriter {
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
            assembly: None,
            pedigree: None,
            header: None,
        }
    }
}

/// Enumeration for describing the orientation of a paired-end read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, Default, Serialize, Deserialize)]
pub enum PeOrientation {
    #[strum(serialize = "3to3")]
    ThreeToThree,
    #[strum(serialize = "5to5")]
    FiveToFive,
    #[strum(serialize = "3to5")]
    ThreeToFive,
    #[strum(serialize = "5to3")]
    FiveToThree,
    #[strum(serialize = "NtoN")]
    #[default]
    Other,
}

impl From<SvType> for PeOrientation {
    fn from(sv_type: SvType) -> Self {
        match sv_type {
            SvType::Del => Self::ThreeToFive,
            SvType::Dup => Self::ThreeToThree,
            SvType::Inv => Self::FiveToFive,
            SvType::Ins | SvType::Cnv | SvType::Bnd => Self::Other,
        }
    }
}

impl FromStr for PeOrientation {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "3to3" => Ok(PeOrientation::ThreeToThree),
            "5to5" => Ok(PeOrientation::FiveToFive),
            "3to5" => Ok(PeOrientation::ThreeToFive),
            "5to3" => Ok(PeOrientation::FiveToThree),
            "NtoN" => Ok(PeOrientation::Other),
            _ => Err(anyhow::anyhow!("Invalid PE orientation: {}", s)),
        }
    }
}

/// Encode the type of an SV
#[derive(
    EnumIter,
    PartialEq,
    Eq,
    Ord,
    PartialOrd,
    Hash,
    Debug,
    Clone,
    Copy,
    Default,
    Display,
    Serialize,
    Deserialize,
)]
pub enum SvType {
    /// Deletion
    #[strum(serialize = "DEL")]
    #[default]
    Del,
    /// Duplication
    #[strum(serialize = "DUP")]
    Dup,
    /// Inversion
    #[strum(serialize = "INV")]
    Inv,
    /// Insertion
    #[strum(serialize = "INS")]
    Ins,
    /// Break-end
    #[strum(serialize = "BND")]
    Bnd,
    /// Copy number variable region
    #[strum(serialize = "CNV")]
    Cnv,
}

impl FromStr for SvType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "DEL" => Ok(SvType::Del),
            "DUP" => Ok(SvType::Dup),
            "INV" => Ok(SvType::Inv),
            "INS" => Ok(SvType::Ins),
            "BND" => Ok(SvType::Bnd),
            "CNV" => Ok(SvType::Cnv),
            _ => Err(anyhow::anyhow!("Invalid SV type: {}", s)),
        }
    }
}

/// Enumeration for describing the SV sub type.
#[derive(
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Debug,
    Clone,
    Copy,
    Default,
    Display,
    Serialize,
    Deserialize,
)]
pub enum SvSubType {
    /// Deletion
    #[strum(serialize = "DEL")]
    #[default]
    Del,
    /// Mobile element deletion
    #[strum(serialize = "DEL:ME")]
    DelMe,
    /// Mobile element deletion (SVA)
    #[strum(serialize = "DEL:ME:SVA")]
    DelMeSva,
    /// Mobile element deletion (L1)
    #[strum(serialize = "DEL:ME:L1")]
    DelMeL1,
    /// Mobile element deletion (ALU)
    #[strum(serialize = "DEL:ME:ALU")]
    DelMeAlu,
    /// Duplication
    #[strum(serialize = "DUP")]
    Dup,
    /// Tandem duplication
    #[strum(serialize = "DUP:TANDEM")]
    DupTandem,
    /// Inversion
    #[strum(serialize = "INV")]
    Inv,
    /// Insertion
    #[strum(serialize = "INS")]
    Ins,
    /// Mobile element insertion
    #[strum(serialize = "INS:ME")]
    InsMe,
    /// Mobile element insertion (SVA)
    #[strum(serialize = "INS:ME:SVA")]
    InsMeSva,
    /// Mobile element insertion (L1)
    #[strum(serialize = "INS:ME:L1")]
    InsMeL1,
    /// Mobile element insertion (ALU)
    #[strum(serialize = "INS:ME:ALU")]
    InsMeAlu,
    /// Break-end
    #[strum(serialize = "BND")]
    Bnd,
    /// Copy number variable region
    #[strum(serialize = "CNV")]
    Cnv,
}

impl From<SvSubType> for SvType {
    fn from(other: SvSubType) -> SvType {
        match other {
            SvSubType::Del
            | SvSubType::DelMe
            | SvSubType::DelMeSva
            | SvSubType::DelMeL1
            | SvSubType::DelMeAlu => SvType::Del,
            SvSubType::Dup | SvSubType::DupTandem => SvType::Dup,
            SvSubType::Inv => SvType::Inv,
            SvSubType::Ins
            | SvSubType::InsMe
            | SvSubType::InsMeSva
            | SvSubType::InsMeL1
            | SvSubType::InsMeAlu => SvType::Ins,
            SvSubType::Bnd => SvType::Bnd,
            SvSubType::Cnv => SvType::Cnv,
        }
    }
}

impl FromStr for SvSubType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "DEL" => Ok(SvSubType::Del),
            "DEL:ME" => Ok(SvSubType::DelMe),
            "DEL:ME:SVA" => Ok(SvSubType::DelMeSva),
            "DEL:ME:L1" => Ok(SvSubType::DelMeL1),
            "DEL:ME:ALU" => Ok(SvSubType::DelMeAlu),
            "DUP" => Ok(SvSubType::Dup),
            "DUP:TANDEM" => Ok(SvSubType::DupTandem),
            "INV" => Ok(SvSubType::Inv),
            "INS" => Ok(SvSubType::Ins),
            "INS:ME" => Ok(SvSubType::InsMe),
            "INS:ME:SVA" => Ok(SvSubType::InsMeSva),
            "INS:ME:L1" => Ok(SvSubType::InsMeL1),
            "INS:ME:ALU" => Ok(SvSubType::InsMeAlu),
            "BND" => Ok(SvSubType::Bnd),
            "CNV" => Ok(SvSubType::Cnv),
            _ => Err(anyhow::anyhow!("Invalid SV sub type: {}", s)),
        }
    }
}

/// Additional infomation for `VarFishStrucvarTsvRecord`.
///
/// Mostly used to store the original `ALT` allele for break-ends.
#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct InfoRecord {
    /// The original ALT allele.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub alt: Option<String>,
}

/// A record, as written out to a VarFish TSV file.
#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct VarFishStrucvarTsvRecord {
    /// The genome assembly/build.
    pub release: String,
    /// Chromosome name of start position.
    pub chromosome: String,
    /// Chromosome number of start position.
    pub chromosome_no: u32,
    /// UCSC bin of start position..
    pub bin: u32,
    /// Chromosome name of end position.
    pub chromosome2: String,
    /// Chromosome number of end position.
    pub chromosome_no2: u32,
    /// UCSC bin of end position.
    pub bin2: u32,

    /// The paired-end orientation.
    pub pe_orientation: PeOrientation,

    /// The start position.
    pub start: i32,
    /// The end position.
    pub end: i32,
    /// Left boundary of confidence around start position.
    pub start_ci_left: i32,
    /// Right boundary of confidence around start position.
    pub start_ci_right: i32,
    /// Left boundary of confidence around end position.
    pub end_ci_left: i32,
    /// Right boundary of confidence around end position.
    pub end_ci_right: i32,

    /// Case ID.
    pub case_id: usize,
    /// Structural variant set ID.
    pub set_id: usize,
    /// Structural variant UUID.
    pub sv_uuid: Uuid,
    /// The callers that this variant was called with.
    pub callers: Vec<String>,
    /// The SV type.
    pub sv_type: SvType,
    /// The SV sub type.
    pub sv_sub_type: SvSubType,

    // /// Additional information (currently not used).
    pub info: InfoRecord,
    /// Genotype call information.
    pub genotype: GenotypeCalls,
}

impl VarFishStrucvarTsvRecord {
    /// Compute reciprocal overlap between `self` and `other`.
    pub fn overlap(&self, other: &Self) -> f32 {
        let s1 = if self.start > 0 { self.start - 1 } else { 0 };
        let e1 = self.end + 1;
        let s2 = if other.start > 0 { other.start - 1 } else { 0 };
        let e2 = other.end;

        let ovl_s = std::cmp::max(s1, s2);
        let ovl_e = std::cmp::min(e1, e2);
        if ovl_e <= ovl_s {
            0.0
        } else {
            let len1 = (e1 - s1) as f32;
            let len2 = (e2 - s2) as f32;
            let ovl_len = (ovl_e - ovl_s) as f32;
            (ovl_len / len1).min(ovl_len / len2)
        }
    }

    /// Merge the other record into this one.
    ///
    /// Currently, the first record's properties are kept except for the `genotype` field.
    /// Here, we take the values with the hights variant count per category (only `gq`,
    /// `pe{c,v}`, `sr{c,v}`).
    ///
    /// Of course, the list of callers is merged as well.
    ///
    /// # Args
    ///
    /// * `other` - The other record to merge into this one.
    ///
    /// # Preconditions
    ///
    /// Both records must have the same SV type and have an appropriate reciprocal overlap.
    fn merge(&mut self, other: &VarFishStrucvarTsvRecord) {
        assert_eq!(self.sv_type, other.sv_type);
        assert_eq!(self.genotype.entries.len(), other.genotype.entries.len());

        self.callers.extend_from_slice(&other.callers);
        self.callers.sort();
        self.callers.dedup();

        for i in 0..self.genotype.entries.len() {
            let mut lhs = self
                .genotype
                .entries
                .get_mut(i)
                .expect("vecs have same size");
            let rhs = other.genotype.entries.get(i).expect("vecs have same size");

            if let Some(lhs_gq) = lhs.gq {
                if let Some(rhs_gq) = rhs.gq {
                    if lhs_gq < rhs_gq {
                        lhs.gq = rhs.gq;
                    }
                }
            } else {
                lhs.gq = rhs.gq;
            }

            if let (Some(_), Some(lhs_pev)) = (lhs.pec, lhs.pec) {
                if let (Some(_), Some(rhs_pev)) = (rhs.pec, rhs.pec) {
                    if lhs_pev < rhs_pev {
                        lhs.pec = rhs.pec;
                        lhs.pev = rhs.pev;
                    }
                }
            } else {
                lhs.pec = rhs.pec;
                lhs.pev = rhs.pev;
            }

            if let (Some(_), Some(lhs_srv)) = (lhs.src, lhs.src) {
                if let (Some(_), Some(rhs_srv)) = (rhs.src, rhs.src) {
                    if lhs_srv < rhs_srv {
                        lhs.src = rhs.src;
                        lhs.srv = rhs.srv;
                    }
                }
            } else {
                lhs.src = rhs.src;
                lhs.srv = rhs.srv;
            }
        }
    }
}

/// Conversion to VCF record.
impl TryInto<VcfRecord> for VarFishStrucvarTsvRecord {
    type Error = anyhow::Error;

    fn try_into(self) -> Result<VcfRecord, Self::Error> {
        use vcf::header::format::key::*;

        let mut genotypes = Vec::new();
        for genotype in &self.genotype.entries {
            genotypes.push(
                [
                    (
                        GENOTYPE,
                        genotype
                            .gt
                            .as_ref()
                            .map(|gt| GenotypeValue::String(gt.clone())),
                    ),
                    (
                        FILTER,
                        genotype.ft.as_ref().map(|ft| {
                            GenotypeValue::StringArray(ft.iter().map(|s| Some(s.clone())).collect())
                        }),
                    ),
                    (
                        CONDITIONAL_GENOTYPE_QUALITY,
                        genotype.gq.as_ref().map(|gq| GenotypeValue::Integer(*gq)),
                    ),
                    (
                        Key::from_str("pec")?,
                        genotype
                            .pec
                            .as_ref()
                            .map(|pec| GenotypeValue::Integer(*pec)),
                    ),
                    (
                        Key::from_str("pev")?,
                        genotype
                            .pev
                            .as_ref()
                            .map(|pev| GenotypeValue::Integer(*pev)),
                    ),
                    (
                        Key::from_str("src")?,
                        genotype
                            .src
                            .as_ref()
                            .map(|src| GenotypeValue::Integer(*src)),
                    ),
                    (
                        Key::from_str("srv")?,
                        genotype
                            .srv
                            .as_ref()
                            .map(|srv| GenotypeValue::Integer(*srv)),
                    ),
                    (
                        Key::from_str("amq")?,
                        genotype
                            .amq
                            .as_ref()
                            .map(|amq| GenotypeValue::Integer(*amq)),
                    ),
                    (
                        GENOTYPE_COPY_NUMBER,
                        genotype.cn.as_ref().map(|cn| GenotypeValue::Integer(*cn)),
                    ),
                    (
                        Key::from_str("anc")?,
                        genotype.anc.as_ref().map(|anc| GenotypeValue::Float(*anc)),
                    ),
                    (
                        Key::from_str("pc")?,
                        genotype.pc.as_ref().map(|pc| GenotypeValue::Integer(*pc)),
                    ),
                ]
                .into_iter()
                .collect(),
            )
        }
        let genotypes = Genotypes::new(
            Keys::try_from(vec![
                GENOTYPE,
                FILTER,
                CONDITIONAL_GENOTYPE_QUALITY,
                Key::from_str("pec")?,
                Key::from_str("pev")?,
                Key::from_str("src")?,
                Key::from_str("srv")?,
                Key::from_str("amq")?,
                GENOTYPE_COPY_NUMBER,
                Key::from_str("anc")?,
                Key::from_str("pc")?,
            ])?,
            genotypes,
        );

        let info = format!(
            "END={};sv_uuid={};callers={};SVTYPE={}",
            self.end,
            self.sv_uuid,
            self.callers.join(","),
            self.sv_sub_type
        );

        let builder = VcfRecord::builder()
            .set_chromosome(self.chromosome.parse()?)
            .set_position(Position::from(self.start as usize))
            .set_reference_bases("N".parse()?)
            .set_info(info.parse()?)
            .set_genotypes(genotypes);

        let builder = if self.sv_sub_type == SvSubType::Bnd {
            builder.set_alternate_bases(self.info.alt.unwrap().parse()?)
        } else {
            builder.set_alternate_bases(format!("<{}>", self.sv_sub_type).parse()?)
        };

        builder.build().map_err(|e| anyhow::anyhow!(e))
    }
}

/// Enumeration for the supported variant callers.
#[derive(Debug, Clone, PartialEq, EnumIter)]
pub enum SvCaller {
    Delly { version: String },
    DragenSv { version: String },
    DragenCnv { version: String },
    Gcnv { version: String },
    Manta { version: String },
    Popdel { version: String },
}

impl SvCaller {
    /// Consider the VCF header and return whether the caller is compatible with `self`.
    fn caller_compatible(&self, header: &VcfHeader) -> bool {
        match self {
            SvCaller::Delly { .. } => {
                self.all_format_defined(header, &["RC", "RCL", "RCR"])
                    && self.all_info_defined(header, &["CHR2", "INSLEN", "HOMLEN"])
            }
            SvCaller::DragenSv { .. } => {
                self.all_format_defined(header, &["PL", "PR", "SR"])
                    && self.all_info_defined(header, &["BND_DEPTH", "MATE_BND_DEPTH"])
                    && self.source_starts_with(header, "DRAGEN")
            }
            SvCaller::Manta { .. } => {
                self.all_format_defined(header, &["PL", "PR", "SR"])
                    && self.all_info_defined(header, &["BND_DEPTH", "MATE_BND_DEPTH"])
                    && self.source_starts_with(header, "GenerateSVCandidates")
            }
            SvCaller::DragenCnv { .. } => {
                self.all_format_defined(header, &["SM", "CN", "BC", "PE"])
                    && self.all_info_defined(header, &["REFLEN", "CIPOS", "CIEND"])
            }
            SvCaller::Gcnv { .. } => {
                self.all_format_defined(header, &["QA", "QS", "QSE", "QSS"])
                    && self.all_info_defined(header, &["END", "AC_Orig", "AN_Orig"])
            }
            SvCaller::Popdel { .. } => {
                self.all_format_defined(header, &["LAD", "DAD", "FL"])
                    && self.all_info_defined(header, &["AF", "IMPRECISE", "SVLEN"])
            }
        }
    }

    /// Extract the caller version from compatible VCF header or first record.
    fn extract_version(
        &self,
        header: &VcfHeader,
        record: &VcfRecord,
    ) -> Result<Self, anyhow::Error> {
        Ok(match self {
            SvCaller::Delly { .. } => SvCaller::Delly {
                version: self.version_from_info_svmethod(record)?,
            },
            SvCaller::DragenSv { .. } => SvCaller::DragenSv {
                version: self.version_from_source_trailing(header)?,
            },
            SvCaller::Gcnv { .. } => SvCaller::Gcnv {
                version: self.version_from_mapping_key(header, "GATKCommandLine")?,
            },
            SvCaller::DragenCnv { .. } => {
                let raw_version = self.version_from_mapping_key(header, "DRAGENVersion")?;
                let mut version = raw_version
                    .split(' ')
                    .nth(1)
                    .expect("Problem extracting version from DRAGENVersion")
                    .chars();
                version.next_back();
                let version = version.as_str().to_owned();
                SvCaller::DragenCnv { version }
            }
            SvCaller::Manta { .. } => SvCaller::Manta {
                version: self.version_from_source_trailing(header)?,
            },
            SvCaller::Popdel { .. } => SvCaller::Popdel {
                version: self.version_from_info_svmethod(record)?,
            },
        })
    }

    /// Parse out version from `INFO/SVMETHOD` in `record`.
    fn version_from_info_svmethod(&self, record: &VcfRecord) -> Result<String, anyhow::Error> {
        let value = record
            .info()
            .get(&InfoKey::Other(InfoKeyOther::from_str("SVMETHOD")?))
            .ok_or(anyhow::anyhow!("Problem with INFO/SVMETHOD field"))?
            .ok_or(anyhow::anyhow!("Problem with INFO/SVMETHOD INFO field"))?;
        if let InfoValue::String(value) = value {
            Ok(value.split('v').last().unwrap().to_string())
        } else {
            anyhow::bail!("Problem with INFO/SVMETHOD INFO field")
        }
    }

    /// Parse out version from `##source=<x> <version>` header.
    fn version_from_source_trailing(&self, header: &VcfHeader) -> Result<String, anyhow::Error> {
        for (key, values) in header.other_records() {
            if key.as_ref() == "source" {
                for value in values {
                    if let HeaderValueOther::String(value) = value {
                        if let Some(version) = value.split(' ').last() {
                            return Ok(version.to_string());
                        }
                    }
                }
            }
        }

        anyhow::bail!("Could not extract ##source header")
    }

    /// Parse out version from `$row_key=<ID=...,CommandLine="...",\
    /// Version="<VERSION>",Date="...">`
    fn version_from_mapping_key(
        &self,
        header: &VcfHeader,
        row_key: &str,
    ) -> Result<String, anyhow::Error> {
        for (key, values) in header.other_records() {
            if key.as_ref() == row_key {
                for value in values {
                    if let HeaderValueOther::Map(_, m) = value {
                        for (k, v) in m.other_fields() {
                            if k.as_ref() == "Version" {
                                return Ok(v.clone());
                            }
                        }
                    }
                }
            }
        }

        anyhow::bail!("Could not extract ##{} header", row_key)
    }

    /// Whether all FORMAT fields are defined.
    fn all_format_defined(&self, header: &VcfHeader, names: &[&str]) -> bool {
        let mut missing = names.iter().map(|s| s.to_string()).collect::<HashSet<_>>();

        for (key, _) in header.formats() {
            missing.remove(key.as_ref());
        }

        missing.is_empty()
    }

    /// Whether all INFO fields are defined.
    fn all_info_defined(&self, header: &VcfHeader, names: &[&str]) -> bool {
        let mut missing = names.iter().map(|s| s.to_string()).collect::<HashSet<_>>();

        for (key, _) in header.infos() {
            missing.remove(key.as_ref());
        }

        missing.is_empty()
    }

    /// Returns whether a `##source=` line's value starts with `prefix`.
    fn source_starts_with(&self, header: &VcfHeader, prefix: &str) -> bool {
        for (key, values) in header.other_records() {
            if key.as_ref() == "source" {
                for value in values {
                    if let HeaderValueOther::String(value) = value {
                        if value.starts_with(prefix) {
                            return true;
                        }
                    }
                }
            }
        }

        false
    }
}

/// Guess the `SvCaller` from the VCF file at the given path.
pub fn guess_sv_caller<P>(p: P) -> Result<SvCaller, anyhow::Error>
where
    P: AsRef<Path>,
{
    let mut reader = vcf::reader::Builder::default().build_from_path(p)?;
    let header: VcfHeader = reader.read_header()?.parse()?;
    let mut records = reader.records(&header);
    let record = records
        .next()
        .ok_or(anyhow::anyhow!("No records found"))??;

    for caller in SvCaller::iter() {
        if caller.caller_compatible(&header) {
            return caller.extract_version(&header, &record);
        }
    }

    anyhow::bail!("Could not guess SV caller from VCF header and first record.");
}

/// Trait that allows the conversion from a VCF record into a `VarFishStrucvarTsvRecord`.
pub trait VcfRecordConverter {
    /// Convert the VCF record into a `VarFishStrucvarTsvRecord`.
    ///
    /// The UUID is passed separately to allow for deterministic UUIDs when necessary.
    ///
    /// # Arguments
    ///
    /// * `record` - The VCF record to convert.
    /// * `uuid` - The UUID to use for the `sv_uuid` field.
    /// * `genome_release` - The genome release to use for the `release` field.
    ///
    /// # Returns
    ///
    /// * `Ok(tsv_record)` if the VCF record was successfully converted.
    ///
    /// # Errors
    ///
    /// * If the VCF record could not be converted.
    fn convert(
        &self,
        vcf_record: &VcfRecord,
        uuid: Uuid,
        genome_release: GenomeRelease,
    ) -> Result<VarFishStrucvarTsvRecord, anyhow::Error> {
        let mut tsv_record = VarFishStrucvarTsvRecord::default();

        self.fill_sv_type(vcf_record, &mut tsv_record)?;
        self.fill_coords(vcf_record, genome_release, &mut tsv_record)?;
        self.fill_callers_uuid(uuid, &mut tsv_record)?;
        self.fill_genotypes(vcf_record, &mut tsv_record)?;
        self.fill_cis(vcf_record, &mut tsv_record)?;
        self.fill_info(vcf_record, &mut tsv_record)?;

        Ok(tsv_record)
    }

    /// Fill the info field of `vcf_record` from `tsv_record`.
    fn fill_info(
        &self,
        vcf_record: &VcfRecord,
        tsv_record: &mut VarFishStrucvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        if let Some(alt) = vcf_record.alternate_bases().deref().first() {
            let ref_allele = vcf_record.reference_bases().to_string();
            let alt_allele = alt.to_string();
            if Breakend::from_ref_alt_str(&ref_allele, &alt_allele).is_ok() {
                tsv_record.info.alt = Some(alt_allele);
            }
        }

        Ok(())
    }

    /// Fill the genotyeps field of `vcf_record` from `tsv_record`.
    ///
    /// # Arguments
    ///
    /// * `vcf_record` - The VCF record to derive the genotypes from.
    /// * `tsv_record` - The TSV record to write the genotypes field to.
    ///
    /// # Returns
    ///
    /// * `Ok(())` if the genotypes were successfully filled.
    ///
    /// # Errors
    ///
    /// * If the genotypes could not be derived from the VCF record.
    fn fill_genotypes(
        &self,
        vcf_record: &VcfRecord,
        tsv_record: &mut VarFishStrucvarTsvRecord,
    ) -> Result<(), anyhow::Error>;

    /// Fill the SV type and sub type of `tsv_record` from `vcf_record`.
    ///
    /// Also fills paired-end orientation (will be overridden later in `fill_coords` in
    /// case of BND).
    ///
    /// # Arguments
    ///
    /// * `vcf_record` - The VCF record to read the SV type and sub type from.
    /// * `tsv_record` - The TSV record to write the SV type and sub type to.
    ///
    /// # Returns
    ///
    /// * `Ok(())` if the SV type and sub type were successfully filled.
    ///
    /// # Errors
    ///
    /// * If the SV type or sub type could not be read from the VCF record.
    fn fill_sv_type(
        &self,
        vcf_record: &VcfRecord,
        tsv_record: &mut VarFishStrucvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        let sv_type = vcf_record
            .info()
            .get(&InfoKey::Standard(InfoKeyStandard::SvType))
            .ok_or_else(|| anyhow::anyhow!("SVTYPE not found"))?
            .ok_or_else(|| anyhow::anyhow!("SVTYPE empty"))?;
        if let InfoValue::String(value) = sv_type {
            tsv_record.sv_sub_type = SvSubType::from_str(value)?;
        } else {
            anyhow::bail!("SVTYPE is not a string");
        }
        tsv_record.sv_type = tsv_record.sv_sub_type.into();

        // Fill paired-end orientation.
        tsv_record.pe_orientation = tsv_record.sv_type.into();

        Ok(())
    }

    /// Fill the coordinates of the given record.
    ///
    /// In the case of breakends, the second chromosome and end position will be parsed from
    /// the (only by assumption) alternate allele in `vcf_record`.  Will also fill the paired-
    /// end orientation in the case of breakend.
    ///
    /// # Preconditions
    ///
    /// - `tsv_record.sv_type` must have been set
    /// - `tsv_record.sv_sub_type` must have been set
    ///
    /// # Arguments
    ///
    /// * `vcf_record` - The VCF record to read the coordinates from.
    /// * `genome_release` - The genome release to use for writing to.
    /// * `tsv_record` - The TSV record to write the coordinates to.
    ///
    /// # Returns
    ///
    /// * `Ok(())` if the coordinates were successfully filled.
    ///
    /// # Errors
    ///
    /// * If the coordinates could not be read from the VCF record.
    fn fill_coords(
        &self,
        vcf_record: &VcfRecord,
        genome_release: GenomeRelease,
        tsv_record: &mut VarFishStrucvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        // Genome release.
        tsv_record.release = match genome_release {
            GenomeRelease::Grch37 => String::from("GRCh37"),
            GenomeRelease::Grch38 => String::from("GRCh38"),
        };

        // Chromosome (of start position if BND).
        let chrom = match vcf_record.chromosome() {
            vcf::record::Chromosome::Name(name) => name.clone(),
            vcf::record::Chromosome::Symbol(_) => panic!("Chromosome symbols are not supported"),
        };
        tsv_record.chromosome = if let Some(chrom) = chrom.strip_prefix("chr") {
            chrom.to_string()
        } else {
            chrom
        };
        // Compute chromosome number.
        tsv_record.chromosome_no = CHROM_TO_CHROM_NO
            .get(&tsv_record.chromosome)
            .copied()
            .unwrap_or(0);

        // Start position.
        let start: usize = vcf_record.position().into();
        tsv_record.start = start as i32;

        // Extract chromosome 2, from alternative allele for BND.  In the case of BND, also extract
        // end position.  For non-BND, set chromosome 2 to chromosome 1.
        let mut end: Option<i32> = None;
        let alleles = &**vcf_record.alternate_bases();
        if alleles.len() != 1 {
            panic!("Only one alternative allele is supported for SVs");
        }
        if let Allele::Breakend(bnd_string) = &alleles[0] {
            let reference = vcf_record
                .reference_bases()
                .iter()
                .map(|c| char::from(*c))
                .collect::<String>();
            let bnd = bnd::Breakend::from_ref_alt_str(&reference, bnd_string)?;

            tsv_record.chromosome2 = bnd.chrom.clone();
            end = Some(bnd.pos);

            // Obtain paired-end orientation.
            tsv_record.pe_orientation = bnd.pe_orientation;
        } else {
            tsv_record.chromosome2 = tsv_record.chromosome.clone();
        }

        // Compute chromosome number of second chromosome.
        tsv_record.chromosome_no2 = CHROM_TO_CHROM_NO
            .get(&tsv_record.chromosome2)
            .copied()
            .unwrap_or(0);

        // End position.  If derived from alternative allele in the case of breakends, use this.
        // Otherwise, try to derive from INFO/END.  If this is not present, use start position
        // as will be the case for insertions.
        tsv_record.end = if let Some(end) = end {
            end
        } else {
            let tmp_end = vcf_record
                .info()
                .get(&InfoKey::Standard(InfoKeyStandard::EndPosition))
                .map(|end| {
                    end.map(|end| match end {
                        vcf::record::info::field::Value::Integer(value) => Ok(Some(*value)),
                        _ => anyhow::bail!("END is not an integer"),
                    })
                })
                .unwrap_or_default()
                .transpose()?
                .unwrap_or_default();

            if let Some(end) = tmp_end {
                end
            } else {
                tsv_record.start
            }
        };

        // Fill bin/bin2.  In the case of BND, the first and second bin are computed indepdenently, even
        // if both positions are on the same chromosome.  Otherwise, the second bin is the same as the
        // first one.  In the latter case, we need to consider insertions as a special case as we simply
        // compute the bin around the single breakpoint position.
        if tsv_record.sv_type == SvType::Bnd {
            tsv_record.bin =
                binning::bin_from_range(tsv_record.start - 1, tsv_record.start)? as u32;
            tsv_record.bin2 = binning::bin_from_range(tsv_record.end - 1, tsv_record.end)? as u32;
        } else {
            if tsv_record.sv_type == SvType::Ins {
                tsv_record.bin =
                    binning::bin_from_range(tsv_record.start - 1, tsv_record.start)? as u32;
            } else {
                tsv_record.bin =
                    binning::bin_from_range(tsv_record.start - 1, tsv_record.end)? as u32;
            }
            tsv_record.bin2 = tsv_record.bin;
        }

        Ok(())
    }

    /// Fill the `callers` and `sv_uuid` field of `tsv_record`.
    fn fill_callers_uuid(
        &self,
        uuid: Uuid,
        tsv_record: &mut VarFishStrucvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        tsv_record.callers = vec![self.caller_version()];
        tsv_record.sv_uuid = uuid;

        Ok(())
    }

    /// Return the caller/version string for this.
    fn caller_version(&self) -> String;

    /// Fill the confidence interval fields.
    fn fill_cis(
        &self,
        vcf_record: &VcfRecord,
        tsv_record: &mut VarFishStrucvarTsvRecord,
    ) -> Result<(), anyhow::Error>;
}

/// Conversion from VCF records to `VarFishStrucvarTsvRecord`.
mod conv {
    use super::GenotypeInfo;
    use super::VarFishStrucvarTsvRecord;
    use super::VcfRecordConverter;

    use lazy_static::__Deref;
    use noodles::vcf::{
        header::info::key::Key as InfoKey, header::info::key::Standard as InfoKeyStandard,
        record::genotypes::genotype::field::value::Value as GenotypeValue,
        record::info::field::value::Value as InfoValue, Record as VcfRecord,
    };

    /// Helper function that extract the CIPOS and CIEND fields from `vcf_record` into `tsv_record`.
    pub fn extract_standard_cis(
        vcf_record: &VcfRecord,
        tsv_record: &mut VarFishStrucvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        let cipos = vcf_record.info().get(&InfoKey::Standard(
            InfoKeyStandard::PositionConfidenceIntervals,
        ));
        // Extract CIPOS; missing field is OK, but if present, must be integer array of length 2.
        if let Some(Some(InfoValue::IntegerArray(cipos))) = cipos {
            if cipos.len() == 2 {
                tsv_record.start_ci_left =
                    cipos[0].ok_or(anyhow::anyhow!("CIPOS[0] is missing"))?;
                tsv_record.start_ci_right =
                    cipos[1].ok_or(anyhow::anyhow!("CIPOS[1] is missing"))?;
            } else {
                anyhow::bail!("CIPOS has wrong number of elements");
            }
        }
        // Extract CIEND; missing field is OK, but if present, must be integer array of length 2.
        let ciend = vcf_record
            .info()
            .get(&InfoKey::Standard(InfoKeyStandard::EndConfidenceIntervals));
        if let Some(Some(InfoValue::IntegerArray(ciend))) = ciend {
            if ciend.len() == 2 {
                tsv_record.end_ci_left = ciend[0].ok_or(anyhow::anyhow!("CIEND[0] is missing"))?;
                tsv_record.end_ci_right = ciend[1].ok_or(anyhow::anyhow!("CIEND[1] is missing"))?;
            } else {
                anyhow::bail!("CIEND has wrong number of elements");
            }
        }

        Ok(())
    }

    pub struct DellyVcfRecordConverter {
        /// The samples from the VCF file.
        pub samples: Vec<String>,
        /// The Delly caller version.
        pub version: String,
    }

    impl DellyVcfRecordConverter {
        pub fn new<T: AsRef<str>>(version: &str, samples: &[T]) -> Self {
            Self {
                samples: samples.iter().map(|s| s.as_ref().to_string()).collect(),
                version: version.to_string(),
            }
        }
    }

    impl VcfRecordConverter for DellyVcfRecordConverter {
        fn caller_version(&self) -> String {
            format!("DELLYv{}", self.version)
        }

        fn fill_cis(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            extract_standard_cis(vcf_record, tsv_record)
        }

        fn fill_genotypes(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            let mut entries: Vec<GenotypeInfo> = vec![Default::default(); self.samples.len()];

            // Extract `FORMAT/*` values.
            for (sample_no, gt) in vcf_record.genotypes().deref().iter().enumerate() {
                entries[sample_no].name = self.samples[sample_no].clone();

                let mut pec = 0;
                let mut src = 0;

                for (key, value) in gt.deref().iter() {
                    match (key.as_ref(), value) {
                        // Obtain `GenotypeInfo::gt` from `FORMAT/GT`.
                        ("GT", Some(GenotypeValue::String(gt))) => {
                            entries[sample_no].gt = Some(gt.to_string());
                        }
                        // Obtain `GenotypeInfo::gq` from `FORMAT/GQ`.
                        ("GQ", Some(GenotypeValue::Integer(gq))) => {
                            entries[sample_no].gq = Some(*gq);
                        }
                        // Obtain `GenotypeInfo::ft` from `FORMAT/FT`.
                        ("FT", Some(GenotypeValue::String(ft))) => {
                            entries[sample_no].ft =
                                Some(ft.split(';').map(|s| s.to_string()).collect());
                        }
                        // Obtain `GenotypeInfo::pev` from `FORMAT/DV`, and accumulate pec.
                        ("DV", Some(GenotypeValue::Integer(dv))) => {
                            entries[sample_no].pev = Some(*dv);
                            pec += *dv;
                        }
                        // Accumulate `FORMAT/DR` into pec.
                        ("DR", Some(GenotypeValue::Integer(dr))) => {
                            pec += *dr;
                        }
                        // Obtain `GenotypeInfo::srv` from `FORMAT/DV`, and accumulate src.
                        ("RV", Some(GenotypeValue::Integer(rv))) => {
                            entries[sample_no].srv = Some(*rv);
                            src += *rv;
                        }
                        // Accumulate `FORMAT/RR` into src.
                        ("RR", Some(GenotypeValue::Integer(rr))) => {
                            src += *rr;
                        }
                        // Ignore all other keys.
                        _ => (),
                    }
                }

                entries[sample_no].pec = Some(pec);
                entries[sample_no].src = Some(src);
            }

            // TODO: get average mapping quality/amq from `maelstrom-core bam-collect-doc` output.

            tsv_record.genotype.entries = entries;

            Ok(())
        }
    }

    pub struct DragenCnvVcfRecordConverter {
        /// The samples from the VCF file.
        pub samples: Vec<String>,
        /// The Dragen SV caller version.
        pub version: String,
    }

    impl DragenCnvVcfRecordConverter {
        pub fn new<T: AsRef<str>>(version: &str, samples: &[T]) -> Self {
            Self {
                samples: samples.iter().map(|s| s.as_ref().to_string()).collect(),
                version: version.to_string(),
            }
        }
    }

    impl VcfRecordConverter for DragenCnvVcfRecordConverter {
        fn caller_version(&self) -> String {
            format!("DRAGEN_CNVv{}", self.version)
        }

        fn fill_cis(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            extract_standard_cis(vcf_record, tsv_record)
        }

        fn fill_genotypes(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            let mut entries: Vec<GenotypeInfo> = vec![Default::default(); self.samples.len()];

            // Extract `FORMAT/*` values.
            for (sample_no, gt) in vcf_record.genotypes().deref().iter().enumerate() {
                entries[sample_no].name = self.samples[sample_no].clone();

                for (key, value) in gt.deref().iter() {
                    match (key.as_ref(), value) {
                        // Obtain `GenotypeInfo::gt` from `FORMAT/GT`.
                        ("GT", Some(GenotypeValue::String(gt))) => {
                            entries[sample_no].gt = Some(gt.to_string());
                        }
                        // Obtain `GenotypeInfo::pev` from `FORMAT/PE`; no pec is computed.
                        ("PE", Some(GenotypeValue::Integer(pe))) => {
                            entries[sample_no].pev = Some(*pe);
                        }
                        // Obtain `GenotypeInfo::cn` from `FORMAT/CN`.
                        ("CN", Some(GenotypeValue::Integer(cn))) => {
                            entries[sample_no].cn = Some(*cn);
                        }
                        // Obtain `GenotypeInfo::pc` from `FORMAT/BC`.
                        ("BC", Some(GenotypeValue::Integer(bc))) => {
                            entries[sample_no].pc = Some(*bc);
                        }
                        // Ignore all other keys.
                        _ => (),
                    }
                }
            }

            // TODO: get average mapping quality/amq from `maelstrom-core bam-collect-doc` output.

            tsv_record.genotype.entries = entries;

            Ok(())
        }
    }

    pub struct DragenSvVcfRecordConverter {
        /// The samples from the VCF file.
        pub samples: Vec<String>,
        /// The Dragen SV caller version.
        pub version: String,
    }

    impl DragenSvVcfRecordConverter {
        pub fn new<T: AsRef<str>>(version: &str, samples: &[T]) -> Self {
            Self {
                samples: samples.iter().map(|s| s.as_ref().to_string()).collect(),
                version: version.to_string(),
            }
        }
    }

    impl VcfRecordConverter for DragenSvVcfRecordConverter {
        fn caller_version(&self) -> String {
            format!("DRAGEN_SVv{}", self.version)
        }

        fn fill_cis(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            extract_standard_cis(vcf_record, tsv_record)
        }

        fn fill_genotypes(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            fill_genotypes_illumina_sv(&self.samples, vcf_record, tsv_record)
        }
    }

    pub struct GcnvVcfRecordConverter {
        /// The samples from the VCF file.
        pub samples: Vec<String>,
        /// The Dragen SV caller version.
        pub version: String,
    }

    impl GcnvVcfRecordConverter {
        pub fn new<T: AsRef<str>>(version: &str, samples: &[T]) -> Self {
            Self {
                samples: samples.iter().map(|s| s.as_ref().to_string()).collect(),
                version: version.to_string(),
            }
        }
    }

    impl VcfRecordConverter for GcnvVcfRecordConverter {
        fn caller_version(&self) -> String {
            format!("GATK_GCNVv{}", self.version)
        }

        fn fill_cis(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            extract_standard_cis(vcf_record, tsv_record)
        }

        fn fill_genotypes(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            let mut entries: Vec<GenotypeInfo> = vec![Default::default(); self.samples.len()];

            // Extract `FORMAT/*` values.
            for (sample_no, gt) in vcf_record.genotypes().deref().iter().enumerate() {
                entries[sample_no].name = self.samples[sample_no].clone();

                for (key, value) in gt.deref().iter() {
                    match (key.as_ref(), value) {
                        // Obtain `GenotypeInfo::gt` from `FORMAT/GT`.
                        ("GT", Some(GenotypeValue::String(gt))) => {
                            entries[sample_no].gt = Some(gt.to_string());
                        }
                        // Obtain `GenotypeInfo::cn` from `FORMAT/CN`.
                        ("CN", Some(GenotypeValue::Integer(cn))) => {
                            entries[sample_no].cn = Some(*cn);
                        }
                        // Obtain `GenotypeInfo::pc` from `FORMAT/NP`.
                        ("NP", Some(GenotypeValue::Integer(np))) => {
                            entries[sample_no].pc = Some(*np);
                        }
                        // Ignore all other keys.
                        _ => (),
                    }
                }
            }

            // TODO: get average mapping quality/amq from `maelstrom-core bam-collect-doc` output.

            tsv_record.genotype.entries = entries;

            Ok(())
        }
    }

    pub struct MantaVcfRecordConverter {
        /// The samples from the VCF file.
        pub samples: Vec<String>,
        /// The Manta caller version.
        pub version: String,
    }

    impl MantaVcfRecordConverter {
        pub fn new<T: AsRef<str>>(version: &str, samples: &[T]) -> Self {
            Self {
                samples: samples.iter().map(|s| s.as_ref().to_string()).collect(),
                version: version.to_string(),
            }
        }
    }

    /// Implementation for filling genotypes of VarFishStrucvarTsvRecord from
    /// Illumina tooling VCF records (Manta/DragenSV).
    fn fill_genotypes_illumina_sv(
        samples: &Vec<String>,
        vcf_record: &VcfRecord,
        tsv_record: &mut VarFishStrucvarTsvRecord,
    ) -> Result<(), anyhow::Error> {
        let mut entries: Vec<GenotypeInfo> = vec![Default::default(); samples.len()];

        // Extract `FORMAT/*` values.
        for (sample_no, gt) in vcf_record.genotypes().deref().iter().enumerate() {
            entries[sample_no].name = samples[sample_no].clone();

            for (key, value) in gt.deref().iter() {
                match (key.as_ref(), value) {
                    // Obtain `GenotypeInfo::gt` from `FORMAT/GT`.
                    ("GT", Some(GenotypeValue::String(gt))) => {
                        entries[sample_no].gt = Some(gt.to_string());
                    }
                    // Obtain `GenotypeInfo::gq` from `FORMAT/GQ`.
                    ("GQ", Some(GenotypeValue::Integer(gq))) => {
                        entries[sample_no].gq = Some(*gq);
                    }
                    // Obtain `GenotypeInfo::ft` from `FORMAT/FT`.
                    ("FT", Some(GenotypeValue::String(ft))) => {
                        entries[sample_no].ft =
                            Some(ft.split(';').map(|s| s.to_string()).collect());
                    }
                    // Obtain `GenotypeInfo::{pev,pec}` from `FORMAT/PR`.
                    ("PR", Some(GenotypeValue::IntegerArray(dv))) => {
                        let ref_ = dv[0].expect("PR[0] is missing");
                        let var = dv[1].expect("PR[1] is missing");
                        entries[sample_no].pec = Some(ref_ + var);
                        entries[sample_no].pev = Some(var);
                    }
                    // Obtain `GenotypeInfo::{pev,pec}` from `FORMAT/PR`.
                    ("SR", Some(GenotypeValue::IntegerArray(sr))) => {
                        let ref_ = sr[0].expect("SR[0] is missing");
                        let var = sr[1].expect("SR[1] is missing");
                        entries[sample_no].src = Some(ref_ + var);
                        entries[sample_no].srv = Some(var);
                    }
                    // Ignore all other keys.
                    _ => (),
                }
            }
        }

        // TODO: get average mapping quality/amq from `maelstrom-core bam-collect-doc` output.

        tsv_record.genotype.entries = entries;

        Ok(())
    }

    impl VcfRecordConverter for MantaVcfRecordConverter {
        fn caller_version(&self) -> String {
            format!("MANTAv{}", self.version)
        }

        fn fill_cis(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            extract_standard_cis(vcf_record, tsv_record)
        }

        fn fill_genotypes(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            fill_genotypes_illumina_sv(&self.samples, vcf_record, tsv_record)
        }
    }

    pub struct PopdelVcfRecordConverter {
        /// The samples from the VCF file.
        pub samples: Vec<String>,
        /// The Dragen SV caller version.
        pub version: String,
    }

    impl PopdelVcfRecordConverter {
        pub fn new<T: AsRef<str>>(version: &str, samples: &[T]) -> Self {
            Self {
                samples: samples.iter().map(|s| s.as_ref().to_string()).collect(),
                version: version.to_string(),
            }
        }
    }

    impl VcfRecordConverter for PopdelVcfRecordConverter {
        fn caller_version(&self) -> String {
            format!("POPDELv{}", self.version)
        }

        fn fill_cis(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            extract_standard_cis(vcf_record, tsv_record)
        }

        fn fill_genotypes(
            &self,
            vcf_record: &VcfRecord,
            tsv_record: &mut VarFishStrucvarTsvRecord,
        ) -> Result<(), anyhow::Error> {
            let mut entries: Vec<GenotypeInfo> = vec![Default::default(); self.samples.len()];

            // Extract `FORMAT/*` values.
            for (sample_no, gt) in vcf_record.genotypes().deref().iter().enumerate() {
                entries[sample_no].name = self.samples[sample_no].clone();

                for (key, value) in gt.deref().iter() {
                    match (key.as_ref(), value) {
                        // Obtain `GenotypeInfo::gt` from `FORMAT/GT`.
                        ("GT", Some(GenotypeValue::String(gt))) => {
                            entries[sample_no].gt = Some(gt.to_string());
                        }
                        // Obtain `GenotypeInfo::gq` from `FORMAT/GQ`.
                        ("GQ", Some(GenotypeValue::Integer(gq))) => {
                            entries[sample_no].gq = Some(*gq);
                        }
                        // Obtain `GenotypeInfo::{pev,pec}` from `FORMAT/DAD[{0,3}]`.
                        ("PR", Some(GenotypeValue::IntegerArray(dad))) => {
                            let ref_ = dad[0].expect("DAD[0] is missing");
                            let var = dad[3].expect("DAD[3] is missing");
                            entries[sample_no].pec = Some(ref_ + var);
                            entries[sample_no].pev = Some(var);
                        }
                        // Ignore all other keys.
                        _ => (),
                    }
                }
            }

            // TODO: get average mapping quality/amq from `maelstrom-core bam-collect-doc` output.

            tsv_record.genotype.entries = entries;

            Ok(())
        }
    }
}

/// Construct a `VcfRecordConverter` for the given caller.
pub fn build_vcf_record_converter<T: AsRef<str>>(
    caller: &SvCaller,
    samples: &[T],
) -> Box<dyn VcfRecordConverter> {
    match caller {
        SvCaller::Delly { version } => {
            Box::new(conv::DellyVcfRecordConverter::new(version, samples))
        }
        SvCaller::Manta { version } => {
            Box::new(conv::MantaVcfRecordConverter::new(version, samples))
        }
        SvCaller::DragenSv { version } => {
            Box::new(conv::DragenSvVcfRecordConverter::new(version, samples))
        }
        SvCaller::DragenCnv { version } => {
            Box::new(conv::DragenCnvVcfRecordConverter::new(version, samples))
        }
        SvCaller::Gcnv { version } => Box::new(conv::GcnvVcfRecordConverter::new(version, samples)),
        SvCaller::Popdel { version } => {
            Box::new(conv::PopdelVcfRecordConverter::new(version, samples))
        }
    }
}

/// Convert the records in the VCF path `path_input` to the JSONL file per contig in `tmp_dir`.
///
/// Note that we will consider the "25 canonical" contigs only (chr1..chr22, chrX, chrY, chrM).
fn run_vcf_to_jsonl(
    path_input: &str,
    tmp_dir: &TempDir,
    rng: &mut StdRng,
) -> Result<(), anyhow::Error> {
    tracing::debug!("opening temporary files in {}", tmp_dir.path().display());
    let mut tmp_files = {
        let mut files = Vec::new();

        for i in 1..=25 {
            let path = tmp_dir.path().join(format!("chrom-{}.jsonl", i));
            tracing::debug!("  path = {}", path.display());
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .expect("could not open temporary file");
            files.push(file);
        }

        files
    };

    tracing::debug!("processing VCF file {}", path_input);
    let sv_caller = guess_sv_caller(path_input)?;
    tracing::debug!("guessed caller/version to be {:?}", &sv_caller);

    let mut reader = vcf::reader::Builder::default().build_from_path(path_input)?;
    let header: VcfHeader = reader.read_header()?.parse()?;

    let samples = header
        .sample_names()
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();
    let converter = build_vcf_record_converter(&sv_caller, &samples);

    let mapping = CHROM_TO_CHROM_NO.deref();
    let mut uuid_buf = [0u8; 16];

    for record in reader.records(&header) {
        rng.fill_bytes(&mut uuid_buf);
        let uuid = Uuid::from_bytes(uuid_buf);

        let record = converter.convert(&record?, uuid, GenomeRelease::Grch37)?;
        if let Some(chromosome_no) = mapping.get(&record.chromosome) {
            let out_jsonl = &mut tmp_files[*chromosome_no as usize - 1];
            jsonl::write(out_jsonl, &record)?;
        } else {
            tracing::warn!(
                "skipping record on chromosome {} (not in canonical set)",
                record.chromosome
            );
        }
    }

    Ok(())
}

/// Read through the JSONL file in `tmp_dir` for contig no. `contig_no` and cluster the records.
///
/// # Arguments
///
/// * `tmp_dir`: The temporary directory containing the JSONL files.
/// * `contig_no`: The contig number (1..25).
/// * `args`: The command line arguments.
///
/// # Returns
///
/// The clustered records for the contig.
fn read_and_cluster_for_contig(
    tmp_dir: &TempDir,
    contig_no: usize,
    args: &Args,
) -> Result<Vec<VarFishStrucvarTsvRecord>, anyhow::Error> {
    let jsonl_path = tmp_dir.path().join(format!("chrom-{}.jsonl", contig_no));
    tracing::debug!("clustering files from {}", jsonl_path.display());
    let mut reader = File::open(jsonl_path).map(BufReader::new)?;

    // The records will be first written to `records`.  The index serves as the ID.
    let mut records: Vec<VarFishStrucvarTsvRecord> = Vec::new();
    // The clusters are build while reading.  The cluster ID is the index into `clusters`.
    // Each cluster consists of the contained record IDs.
    let mut clusters: Vec<Vec<usize>> = vec![];
    // We use the dynamic inteval trees from bio-rust.
    let mut tree: IntervalTree<i32, usize> = IntervalTree::new();

    // Read through the JSONL file, store records, and create clusters.
    loop {
        match jsonl::read::<_, VarFishStrucvarTsvRecord>(&mut reader) {
            Ok(record) => {
                let slack = match record.sv_type {
                    SvType::Ins => args.slack_ins,
                    SvType::Bnd => args.slack_bnd,
                    _ => 0,
                };
                let query = match record.sv_type {
                    SvType::Ins | SvType::Bnd => {
                        let query_start = std::cmp::max(1, record.start - 1 - slack);
                        let query_end = record.start + slack;
                        query_start..query_end
                    }
                    _ => {
                        let query_start = std::cmp::max(1, record.start - 1 - slack);
                        let query_end = record.end + slack;
                        query_start..query_end
                    }
                };
                let mut found_any_cluster = false;
                for mut it_tree in tree.find_mut(&query) {
                    let cluster_idx = *it_tree.data();
                    let mut match_all_in_cluster = true;
                    for it_cluster in &clusters[cluster_idx] {
                        let record_id = it_cluster;
                        let match_this_range = match record.sv_type {
                            SvType::Bnd | SvType::Ins => true,
                            _ => {
                                let ovl = record.overlap(&records[*record_id]);
                                assert!(ovl >= 0f32);
                                ovl >= args.min_overlap
                            }
                        };
                        let match_this_sv_type = record.sv_type == records[*record_id].sv_type;
                        match_all_in_cluster =
                            match_all_in_cluster && match_this_range && match_this_sv_type;
                    }
                    if match_all_in_cluster {
                        // extend cluster
                        clusters[cluster_idx].push(records.len());
                        found_any_cluster = true;
                        break;
                    }
                }
                if !found_any_cluster {
                    // create new cluster
                    match record.sv_type {
                        SvType::Ins | SvType::Bnd => {
                            tree.insert((record.start - 1)..record.start, clusters.len());
                        }
                        _ => {
                            tree.insert((record.start - 1)..record.end, clusters.len());
                        }
                    };
                    clusters.push(vec![records.len()]);
                }
                // always register the record
                records.push(record);
            }
            Err(jsonl::ReadError::Eof) => {
                break; // all done
            }
            _ => anyhow::bail!("Problem reading JSONL file"),
        }
    }

    // Iterate over the clusters, merge into the first one found.
    let mut result = Vec::new();
    for cluster in &clusters {
        let mut record = records[cluster[0]].clone();
        for other_idx in cluster[1..].iter() {
            record.merge(&records[*other_idx]);
        }
        result.push(record);
    }

    // Finally, sort records by start position and write out.
    result.sort_by(|a, b| a.start.cmp(&b.start));

    Ok(result)
}

/// Run the annotation with the given `Write` within the `VcfWriter`.
///
/// # Arguments
///
/// * `writer`: The VCF writer.
/// * `args`: The command line arguments.
fn run_with_writer(
    writer: &mut dyn AnnotatedVcfWriter,
    args: &Args,
    pedigree: &PedigreeByName,
) -> Result<(), anyhow::Error> {
    let mut rng = if let Some(rng_seed) = args.rng_seed {
        StdRng::seed_from_u64(rng_seed)
    } else {
        StdRng::from_entropy()
    };

    // Create temporary directory.  We will create one temporary file (containing `jsonl`
    // seriealized `VarFishStrucvarTsvRecord`s) for each SV type and contig.
    let tmp_dir = TempDir::new("mehari")?;

    // Read through input VCF files and write out to temporary files.
    tracing::info!("Input VCF files to temporary files...");
    for path_input in args.path_input_vcf.iter() {
        run_vcf_to_jsonl(path_input, &tmp_dir, &mut rng)?;
    }
    tracing::info!("... done converting input files.");

    // Generate output header and write to `writer`.
    tracing::info!("Write output header...");
    let header_out = vcf_header::build(
        args.genome_release
            .expect("genome release must be known here")
            .into(),
        pedigree,
        &Utc::now().date_naive(),
    )?;
    writer.write_header(&header_out)?;

    tracing::info!("Clustering SVs to output...");
    // Read through temporary files by contig, cluster by overlap as configured, and write to `writer`.
    for contig_no in 1..=25 {
        tracing::info!("  contig: {}", CANONICAL[contig_no - 1]);
        let clusters = read_and_cluster_for_contig(&tmp_dir, contig_no, args)?;
        for record in clusters {
            writer.write_record(&record.try_into()?)?;
        }
    }
    tracing::info!("... done clustering SVs to output");

    Ok(())
}

/// Main entry point for `annotate strucvars` sub command.
pub fn run(_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("config = {:#?}", &args);
    // Load the pedigree.
    tracing::info!("Loading pedigree...");
    let pedigree = PedigreeByName::from_path(&args.path_input_ped)?;
    tracing::info!("... done loading pedigree");

    // Guess genome release from contigs in first VCF header.
    let assembly = args.genome_release.map(|gr| match gr {
        GenomeRelease::Grch37 => Assembly::Grch37p10, // has chrMT!
        GenomeRelease::Grch38 => Assembly::Grch38,
    });
    let assembly = {
        let mut reader = VariantReaderBuilder::default().build_from_path(
            args.path_input_vcf
                .first()
                .expect("must have at least input VCF"),
        )?;
        let header = reader.read_header()?;
        guess_assembly(&header, false, assembly)?
    };
    tracing::info!("Determined input assembly to be {:?}", &assembly);
    let args = Args {
        genome_release: Some(assembly.into()),
        ..args.clone()
    };

    if let Some(path_output_vcf) = &args.output.path_output_vcf {
        if path_output_vcf.ends_with(".vcf.gz") || path_output_vcf.ends_with(".vcf.bgzf") {
            let mut writer = VcfWriter::new(
                File::create(path_output_vcf)
                    .map(BufWriter::new)
                    .map(BgzfWriter::new)?,
            );
            writer.set_assembly(assembly);
            writer.set_pedigree(&pedigree);
            run_with_writer(&mut writer, &args, &pedigree)?;
        } else {
            let mut writer = VcfWriter::new(File::create(path_output_vcf).map(BufWriter::new)?);
            writer.set_assembly(assembly);
            writer.set_pedigree(&pedigree);
            run_with_writer(&mut writer, &args, &pedigree)?;
        }
    } else {
        let path_output_tsv = args
            .output
            .path_output_tsv
            .as_ref()
            .expect("tsv path must be set; vcf and tsv are mutually exclusive, vcf unset");
        let mut writer = VarFishStrucvarTsvWriter::with_path(path_output_tsv);
        writer.set_assembly(assembly);
        writer.set_pedigree(&pedigree);

        run_with_writer(&mut writer, &args, &pedigree)?;
    }

    Ok(())
}

/// Parsing and representation of breakends.
///
/// From the VCF4.2 docs (Section 5.4). Single breake-ends are currently not supported.
///
/// ```text
/// 2  321681 bndW G G]17:198982] 6 PASS SVTYPE=BND 5to5    leading       left_open
/// 17 198982 bndY A A]2:321681]  6 PASS SVTYPE=BND 5to5    leading       left_open
///
/// 2  321682 bndV T ]13:123456]T 6 PASS SVTYPE=BND 3to5   !leading       left_open
/// 13 123456 bndU C C[2:321682[  6 PASS SVTYPE=BND 5to3    leading      !left_open
///
/// 13 123457 bndX A [17:198983[A 6 PASS SVTYPE=BND 3to3   !leading      !left_open
/// 17 198983 bndZ C [13:123457[C 6 PASS SVTYPE=BND 3to3   !leading      !left_open
/// ```
pub mod bnd {
    use super::PeOrientation;

    #[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
    pub struct Breakend {
        /// Name of the chromosome.
        pub chrom: String,
        /// 1-based position on the chromosome.
        pub pos: i32,
        /// The reference base from the VCF field.
        pub ref_base: String,
        /// The alternative base from the VCF field.
        pub alt_base: String,
        /// Whether the alternative base is leading.
        pub leading_base: bool,
        /// Whether the square brackets are left open.
        pub left_open: bool,
        /// The derived paired-end orientation.
        pub pe_orientation: PeOrientation,
    }

    impl Breakend {
        /// Construct from reference and alternative bases
        pub fn from_ref_alt_str(reference: &str, alternative: &str) -> Result<Self, anyhow::Error> {
            let mut alt_chars = alternative.chars();
            let first = alt_chars.next().ok_or(anyhow::anyhow!(
                "ALT string must not be empty: {}",
                reference
            ))?;
            let last = alt_chars.next_back().ok_or(anyhow::anyhow!(
                "ALT string must have at least two characters: {}",
                reference
            ))?;

            let leading_bracket = "[]".contains(first);
            let trailing_bracket = "[]".contains(last);
            if leading_bracket == trailing_bracket {
                anyhow::bail!(
                    "Cannot have trailing AND leading bracket in ALT: {}",
                    alternative
                );
            }
            // Determine alternative base and ensure that `alt_chars` only contains the chrom/pos pair.
            let left_open;
            let alt_base;
            if leading_bracket {
                left_open = first == ']';
                alt_base = last.to_string();
                alt_chars.next_back();
            } else {
                left_open = last == ']';
                alt_base = first.to_string();
                alt_chars.next();
            }

            // Split chrom/pos and parse integer position.
            let chrom_pos = alt_chars.as_str();
            let mut chrom_pos = chrom_pos.split(':');
            let chrom = chrom_pos
                .next()
                .ok_or(anyhow::anyhow!(
                    "Chromosome must be present in ALT: {}",
                    alternative
                ))?
                .to_string();
            let pos = chrom_pos.next().ok_or(anyhow::anyhow!(
                "Position must be present in ALT: {}",
                alternative
            ))?;
            let pos: i32 = pos.parse()?;

            // Determine PE orientation.
            let pe_orientation = match (!leading_bracket, left_open) {
                (true, true) => PeOrientation::FiveToFive,
                (true, false) => PeOrientation::FiveToThree,
                (false, true) => PeOrientation::ThreeToFive,
                (false, false) => PeOrientation::ThreeToThree,
            };

            Ok(Self {
                chrom,
                pos,
                ref_base: reference.to_string(),
                alt_base,
                leading_base: !leading_bracket,
                left_open,
                pe_orientation,
            })
        }
    }
}

#[cfg(test)]
mod test {
    use std::fs::File;

    use chrono::NaiveDate;
    use hgvs::static_data::Assembly;
    use linked_hash_map::LinkedHashMap;
    use noodles::vcf;
    use pretty_assertions::assert_eq;
    use temp_testdir::TempDir;
    use uuid::Uuid;

    use crate::{
        annotate::{
            seqvars::AnnotatedVcfWriter,
            strucvars::{
                GenotypeCalls, GenotypeInfo, InfoRecord, PeOrientation, SvCaller, SvSubType,
                SvType, VarFishStrucvarTsvRecord,
            },
        },
        common::GenomeRelease,
        ped::{Disease, Individual, PedigreeByName, Sex},
    };

    use super::{
        bnd::Breakend,
        build_vcf_record_converter,
        conv::{
            DellyVcfRecordConverter, DragenCnvVcfRecordConverter, DragenSvVcfRecordConverter,
            GcnvVcfRecordConverter, MantaVcfRecordConverter, PopdelVcfRecordConverter,
        },
        guess_sv_caller, vcf_header, VarFishStrucvarTsvWriter, VcfHeader, VcfRecord,
        VcfRecordConverter,
    };

    /// Test for the parsing of breakend alleles.
    #[test]
    fn parse_bnd() -> Result<(), anyhow::Error> {
        assert_eq!(
            Breakend {
                chrom: String::from("17"),
                pos: 198982,
                ref_base: String::from("G"),
                alt_base: String::from("A"),
                leading_base: true,
                left_open: true,
                pe_orientation: PeOrientation::FiveToFive,
            },
            Breakend::from_ref_alt_str("G", "A]17:198982]")?,
        );
        assert_eq!(
            Breakend {
                chrom: String::from("17"),
                pos: 198982,
                ref_base: String::from("G"),
                alt_base: String::from("A"),
                leading_base: true,
                left_open: false,
                pe_orientation: PeOrientation::FiveToThree,
            },
            Breakend::from_ref_alt_str("G", "A[17:198982[")?,
        );
        assert_eq!(
            Breakend {
                chrom: String::from("17"),
                pos: 198982,
                ref_base: String::from("G"),
                alt_base: String::from("A"),
                leading_base: false,
                left_open: true,
                pe_orientation: PeOrientation::ThreeToFive,
            },
            Breakend::from_ref_alt_str("G", "]17:198982]A")?,
        );
        assert_eq!(
            Breakend {
                chrom: String::from("17"),
                pos: 198982,
                ref_base: String::from("G"),
                alt_base: String::from("A"),
                leading_base: false,
                left_open: false,
                pe_orientation: PeOrientation::ThreeToThree,
            },
            Breakend::from_ref_alt_str("G", "[17:198982[A")?,
        );

        Ok(())
    }

    /// Helper function that tests the VCF to JSONL conversion.
    fn run_test_vcf_to_jsonl(
        path_input_vcf: &str,
        path_expected_jsonl: &str,
        converter: &dyn VcfRecordConverter,
    ) -> Result<(), anyhow::Error> {
        let out_file_name = "out.jsonl";

        let temp = TempDir::default();
        let out_jsonl = File::create(temp.join(out_file_name))?;

        let mut reader = vcf::reader::Builder::default().build_from_path(path_input_vcf)?;
        let header_in = reader.read_header()?.parse()?;

        // Setup deterministic bytes for UUID generation.
        let mut bytes = [
            0xa1, 0xa2, 0xa3, 0xa4, 0xb1, 0xb2, 0xc1, 0xc2, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6,
            0xd7, 0xd8,
        ];

        // Convert all VCF records to JSONL, incrementing the first byte for deterministic UUID
        // generation.
        let mut records = reader.records(&header_in);
        loop {
            if let Some(record) = records.next() {
                let uuid = Uuid::from_bytes(bytes);
                let record = converter.convert(&record?, uuid, GenomeRelease::Grch37)?;
                jsonl::write(&out_jsonl, &record)?;
            } else {
                break; // all done
            }

            bytes[0] += 1;
        }
        drop(out_jsonl);

        let expected = std::fs::read_to_string(path_expected_jsonl)?;
        let actual = std::fs::read_to_string(temp.join(out_file_name))?;

        assert_eq!(expected, actual);

        Ok(())
    }

    /// Helper that returns sample names from VCF.
    fn vcf_samples(path: &str) -> Result<Vec<String>, anyhow::Error> {
        let mut reader = vcf::reader::Builder::default().build_from_path(path)?;
        let header: VcfHeader = reader.read_header()?.parse()?;
        Ok(header
            .sample_names()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>())
    }

    #[test]
    fn vcf_to_jsonl_delly_min() -> Result<(), anyhow::Error> {
        let path_input_vcf = "tests/data/annotate/strucvars/delly2-min.vcf";
        let samples = vcf_samples(path_input_vcf)?;

        run_test_vcf_to_jsonl(
            path_input_vcf,
            "tests/data/annotate/strucvars/delly2-min.out.jsonl",
            &DellyVcfRecordConverter::new("1.1.3", &samples),
        )
    }

    #[test]
    fn vcf_to_jsonl_dragen_sv_min() -> Result<(), anyhow::Error> {
        let path_input_vcf = "tests/data/annotate/strucvars/dragen-sv-min.vcf";
        let samples = vcf_samples(path_input_vcf)?;

        run_test_vcf_to_jsonl(
            path_input_vcf,
            "tests/data/annotate/strucvars/dragen-sv-min.out.jsonl",
            &DragenSvVcfRecordConverter::new("07.021.624.3.10.4", &samples),
        )
    }

    #[test]
    fn vcf_to_jsonl_dragen_cnv_min() -> Result<(), anyhow::Error> {
        let path_input_vcf = "tests/data/annotate/strucvars/dragen-cnv-min.vcf";
        let samples = vcf_samples(path_input_vcf)?;

        run_test_vcf_to_jsonl(
            path_input_vcf,
            "tests/data/annotate/strucvars/dragen-cnv-min.out.jsonl",
            &DragenCnvVcfRecordConverter::new("07.021.624.3.10.4", &samples),
        )
    }

    #[test]
    fn vcf_to_jsonl_gcnv_min() -> Result<(), anyhow::Error> {
        let path_input_vcf = "tests/data/annotate/strucvars/gcnv-min.vcf";
        let samples = vcf_samples(path_input_vcf)?;

        run_test_vcf_to_jsonl(
            path_input_vcf,
            "tests/data/annotate/strucvars/gcnv-min.out.jsonl",
            &GcnvVcfRecordConverter::new("4.3.0.0", &samples),
        )
    }

    #[test]
    fn vcf_to_jsonl_manta_min() -> Result<(), anyhow::Error> {
        let path_input_vcf = "tests/data/annotate/strucvars/manta-min.vcf";
        let samples = vcf_samples(path_input_vcf)?;

        run_test_vcf_to_jsonl(
            path_input_vcf,
            "tests/data/annotate/strucvars/manta-min.out.jsonl",
            &MantaVcfRecordConverter::new("1.6.0", &samples),
        )
    }

    #[test]
    fn vcf_to_jsonl_popdel_min() -> Result<(), anyhow::Error> {
        let path_input_vcf = "tests/data/annotate/strucvars/popdel-min.vcf";
        let samples = vcf_samples(path_input_vcf)?;

        run_test_vcf_to_jsonl(
            path_input_vcf,
            "tests/data/annotate/strucvars/popdel-min.out.jsonl",
            &PopdelVcfRecordConverter::new("1.1.2", &samples),
        )
    }

    #[test]
    fn vcf_to_jsonl_with_detection() -> Result<(), anyhow::Error> {
        let keys = &[
            "delly2",
            "dragen-cnv",
            "dragen-sv",
            "gcnv",
            "manta",
            "popdel",
        ];

        for key in keys {
            let path_input_vcf = format!("tests/data/annotate/strucvars/{}-min.vcf", key);
            let path_expected_txt = format!("tests/data/annotate/strucvars/{}-min.out.jsonl", key);
            let samples = vcf_samples(&path_input_vcf)?;
            let sv_caller = guess_sv_caller(&path_input_vcf)?;
            let converter = build_vcf_record_converter(&sv_caller, &samples);

            run_test_vcf_to_jsonl(&path_input_vcf, &path_expected_txt, converter.as_ref())?;
        }

        Ok(())
    }

    #[test]
    fn guess_sv_caller_delly() -> Result<(), anyhow::Error> {
        let sv_caller = guess_sv_caller("tests/data/annotate/strucvars/delly2-min.vcf")?;

        assert_eq!(
            sv_caller,
            SvCaller::Delly {
                version: String::from("1.1.3")
            }
        );

        Ok(())
    }

    #[test]
    fn guess_sv_caller_dragen_sv() -> Result<(), anyhow::Error> {
        let sv_caller = guess_sv_caller("tests/data/annotate/strucvars/dragen-sv-min.vcf")?;

        assert_eq!(
            sv_caller,
            SvCaller::DragenSv {
                version: String::from("07.021.624.3.10.4")
            }
        );

        Ok(())
    }

    #[test]
    fn guess_sv_caller_dragen_cnv() -> Result<(), anyhow::Error> {
        let sv_caller = guess_sv_caller("tests/data/annotate/strucvars/dragen-cnv-min.vcf")?;

        assert_eq!(
            sv_caller,
            SvCaller::DragenCnv {
                version: String::from("07.021.624.3.10.4")
            }
        );

        Ok(())
    }

    #[test]
    fn guess_sv_caller_gcnv() -> Result<(), anyhow::Error> {
        let sv_caller = guess_sv_caller("tests/data/annotate/strucvars/gcnv-min.vcf")?;

        assert_eq!(
            sv_caller,
            SvCaller::Gcnv {
                version: String::from("4.3.0.0")
            }
        );

        Ok(())
    }

    #[test]
    fn guess_sv_caller_manta() -> Result<(), anyhow::Error> {
        let sv_caller = guess_sv_caller("tests/data/annotate/strucvars/manta-min.vcf")?;

        assert_eq!(
            sv_caller,
            SvCaller::Manta {
                version: String::from("1.6.0")
            }
        );

        Ok(())
    }

    #[test]
    fn guess_sv_caller_popdel() -> Result<(), anyhow::Error> {
        let sv_caller = guess_sv_caller("tests/data/annotate/strucvars/popdel-min.vcf")?;

        assert_eq!(
            sv_caller,
            SvCaller::Popdel {
                version: String::from("1.1.2")
            }
        );

        Ok(())
    }

    #[test]
    fn build_vcf_header_37_no_pedigree() -> Result<(), anyhow::Error> {
        let header = vcf_header::build(
            Assembly::Grch37p10,
            &Default::default(),
            &NaiveDate::from_ymd_opt(2015, 3, 14).unwrap(),
        )?;

        let mut writer = vcf::Writer::new(Vec::new());
        writer.write_header(&header)?;
        let actual = std::str::from_utf8(&writer.get_ref()[..])?;

        let expected =
            std::fs::read_to_string("tests/data/annotate/strucvars/header-grch37-noped.vcf")?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn build_vcf_header_37_trio() -> Result<(), anyhow::Error> {
        let header = vcf_header::build(
            Assembly::Grch37p10,
            &example_trio(),
            &NaiveDate::from_ymd_opt(2015, 3, 14).unwrap(),
        )?;

        let mut writer = vcf::Writer::new(Vec::new());
        writer.write_header(&header)?;
        let actual = std::str::from_utf8(&writer.get_ref()[..])?;

        let expected =
            std::fs::read_to_string("tests/data/annotate/strucvars/header-grch37-trio.vcf")?;

        assert_eq!(actual, expected);

        Ok(())
    }

    /// Generate example trio data.
    fn example_trio() -> PedigreeByName {
        let individuals = LinkedHashMap::from_iter(
            vec![
                (
                    String::from("index"),
                    Individual {
                        family: String::from("FAM"),
                        name: String::from("index"),
                        father: Some(String::from("father")),
                        mother: Some(String::from("mother")),
                        sex: Sex::Female,
                        disease: Disease::Affected,
                    },
                ),
                (
                    String::from("father"),
                    Individual {
                        family: String::from("FAM"),
                        name: String::from("father"),
                        father: None,
                        mother: None,
                        sex: Sex::Male,
                        disease: Disease::Unaffected,
                    },
                ),
                (
                    String::from("mother"),
                    Individual {
                        family: String::from("FAM"),
                        name: String::from("mother"),
                        father: None,
                        mother: None,
                        sex: Sex::Female,
                        disease: Disease::Unaffected,
                    },
                ),
            ]
            .into_iter(),
        );
        PedigreeByName { individuals }
    }

    #[test]
    fn build_vcf_header_38_no_pedigree() -> Result<(), anyhow::Error> {
        let header = vcf_header::build(
            Assembly::Grch38,
            &Default::default(),
            &NaiveDate::from_ymd_opt(2015, 3, 14).unwrap(),
        )?;

        let mut writer = vcf::Writer::new(Vec::new());
        writer.write_header(&header)?;
        let actual = std::str::from_utf8(&writer.get_ref()[..])?;

        let expected =
            std::fs::read_to_string("tests/data/annotate/strucvars/header-grch38-noped.vcf")?;

        assert_eq!(actual, expected);

        Ok(())
    }

    /// Example record for writing out.
    fn example_records() -> Vec<VarFishStrucvarTsvRecord> {
        // Setup deterministic bytes for UUID generation.
        let bytes = [
            0xa1, 0xa2, 0xa3, 0xa4, 0xb1, 0xb2, 0xc1, 0xc2, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6,
            0xd7, 0xd8,
        ];
        let mut bytes2 = bytes;
        bytes2[0] += 1;

        vec![
            VarFishStrucvarTsvRecord {
                release: String::from("GRCh37"),
                chromosome: String::from("1"),
                chromosome_no: 1,
                bin: 512,
                chromosome2: String::from("1"),
                chromosome_no2: 1,
                bin2: 512,
                pe_orientation: PeOrientation::FiveToFive,
                start: 1042,
                end: 2042,
                start_ci_left: 0,
                start_ci_right: 0,
                end_ci_left: 0,
                end_ci_right: 0,
                case_id: 0,
                set_id: 0,
                sv_uuid: Uuid::from_bytes(bytes),
                callers: vec![String::from("MANTAv1.1.2")],
                sv_type: SvType::Inv,
                sv_sub_type: SvSubType::Inv,
                info: Default::default(),
                genotype: GenotypeCalls {
                    entries: vec![
                        GenotypeInfo {
                            name: String::from("index"),
                            gt: Some(String::from("0/1")),
                            ft: Some(vec![String::from("PASS")]),
                            gq: Some(99),
                            pec: Some(143),
                            pev: Some(43),
                            src: Some(143),
                            srv: Some(43),
                            amq: Some(99),
                            cn: Some(1),
                            anc: Some(0.5),
                            pc: Some(10),
                        },
                        GenotypeInfo {
                            name: String::from("father"),
                            gt: Some(String::from("0/0")),
                            ft: Some(vec![String::from("PASS")]),
                            gq: Some(98),
                            pec: Some(43),
                            pev: Some(0),
                            src: Some(44),
                            srv: Some(0),
                            amq: Some(98),
                            cn: Some(2),
                            anc: Some(1.0),
                            pc: Some(10),
                        },
                        GenotypeInfo {
                            name: String::from("mother"),
                            gt: Some(String::from("0/0")),
                            ft: Some(vec![String::from("PASS")]),
                            gq: Some(97),
                            pec: Some(32),
                            pev: Some(0),
                            src: Some(33),
                            srv: Some(0),
                            amq: Some(97),
                            cn: Some(2),
                            anc: Some(1.0),
                            pc: Some(10),
                        },
                    ],
                },
            },
            VarFishStrucvarTsvRecord {
                release: String::from("GRCh37"),
                chromosome: String::from("1"),
                chromosome_no: 1,
                bin: 512,
                chromosome2: String::from("17"),
                chromosome_no2: 1,
                bin2: 1234,
                pe_orientation: PeOrientation::FiveToFive,
                start: 10_1042,
                end: 198_982,
                start_ci_left: 0,
                start_ci_right: 0,
                end_ci_left: 0,
                end_ci_right: 0,
                case_id: 0,
                set_id: 0,
                sv_uuid: Uuid::from_bytes(bytes2),
                callers: vec![String::from("MANTAv1.1.2")],
                sv_type: SvType::Bnd,
                sv_sub_type: SvSubType::Bnd,
                info: InfoRecord {
                    alt: Some(String::from("]17:198982]A")),
                },
                genotype: GenotypeCalls {
                    entries: vec![
                        GenotypeInfo {
                            name: String::from("index"),
                            gt: Some(String::from("0/1")),
                            ft: Some(vec![String::from("PASS")]),
                            gq: Some(99),
                            pec: Some(143),
                            pev: Some(43),
                            src: Some(143),
                            srv: Some(43),
                            amq: Some(99),
                            cn: Some(1),
                            anc: Some(0.5),
                            pc: Some(10),
                        },
                        GenotypeInfo {
                            name: String::from("father"),
                            gt: Some(String::from("0/1")),
                            ft: Some(vec![String::from("PASS")]),
                            gq: Some(99),
                            pec: Some(143),
                            pev: Some(43),
                            src: Some(143),
                            srv: Some(43),
                            amq: Some(99),
                            cn: Some(1),
                            anc: Some(0.5),
                            pc: Some(10),
                        },
                        GenotypeInfo {
                            name: String::from("mother"),
                            gt: Some(String::from("0/1")),
                            ft: Some(vec![String::from("PASS")]),
                            gq: Some(99),
                            pec: Some(143),
                            pev: Some(43),
                            src: Some(143),
                            srv: Some(43),
                            amq: Some(99),
                            cn: Some(1),
                            anc: Some(0.5),
                            pc: Some(10),
                        },
                    ],
                },
            },
        ]
    }

    #[test]
    fn write_vcf_from_varfish_records() -> Result<(), anyhow::Error> {
        let header = vcf_header::build(
            Assembly::Grch38,
            &example_trio(),
            &NaiveDate::from_ymd_opt(2015, 3, 14).unwrap(),
        )?;

        let mut writer = vcf::Writer::new(Vec::new());
        writer.write_header(&header)?;

        for varfish_record in example_records() {
            let vcf_record: VcfRecord = varfish_record.try_into()?;
            writer.write_record(&vcf_record)?;
        }

        let actual = std::str::from_utf8(&writer.get_ref()[..])?;

        let expected = std::fs::read_to_string("tests/data/annotate/strucvars/example-grch38.vcf")?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn write_tsv_from_varfish_records() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();

        // scope for writer
        {
            let mut writer = VarFishStrucvarTsvWriter::with_path(temp.join("out.tsv"));
            writer.set_assembly(Assembly::Grch37p10);
            writer.set_pedigree(&example_trio());

            for varfish_record in example_records() {
                let vcf_record: VcfRecord = varfish_record.try_into()?;
                writer.write_record(&vcf_record)?;
            }
        }

        let expected = std::fs::read_to_string("tests/data/annotate/strucvars/example-grch38.tsv")?;
        let actual = std::fs::read_to_string(temp.join("out.tsv"))?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
