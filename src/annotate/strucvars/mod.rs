//! Annotation of structural variant VCF files.

use std::io::Write;
use std::path::Path;
use std::str::FromStr;
use std::{fs::File, io::BufWriter};

use crate::annotate::seqvars::HgncRecord;
use crate::common::GenomeRelease;
use crate::ped::PedigreeByName;
use clap::{Args as ClapArgs, Parser};
use flate2::write::GzEncoder;
use flate2::Compression;
use hgvs::static_data::Assembly;
use noodles::bgzf::Writer as BgzfWriter;
use noodles::vcf::record::alternate_bases::Allele;
use noodles::vcf::record::Position;
use noodles::vcf::{
    self, header::info::key::Key as InfoKey, header::info::key::Standard as InfoKeyStandard,
    Header as VcfHeader, Record as VcfRecord, Writer as VcfWriter,
};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use strum::{Display, EnumIter};
use tempdir::TempDir;
use uuid::Uuid;

use super::seqvars::{binning, AnnotatedVcfWriter, CHROM_TO_CHROM_NO};

/// Command line arguments for `annotate strucvars` sub command.
#[derive(Parser, Debug)]
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
    pub path_input_ped: Option<String>,
    /// Path to the input VCF files.
    #[arg(long)]
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
}

/// Command line arguments to enforce either `--path-output-vcf` or `--path-output-tsv`.
#[derive(Debug, ClapArgs)]
#[group(required = true, multiple = false)]
pub struct PathOutput {
    /// Path to the output VCF file.
    #[arg(long)]
    pub path_output_vcf: Option<String>,

    /// Path to the output TSV file (for import into VarFish).
    #[arg(long, requires = "path-input-ped")]
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

fn build_header() -> VcfHeader {
    let mut builder = vcf::Header::builder();

    todo!();

    builder.build()
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
    /// Mapping from HGNC id to record with gene symbol and gene identifiers.
    hgnc_map: Option<FxHashMap<String, HgncRecord>>,
}

/// Per-genotype call information.
#[derive(Debug, Default, Serialize, Deserialize, PartialEq, Clone)]
pub struct GenotypeInfo {
    /// Sample name.
    pub name: String,
    /// Genotype value.
    pub gt: Option<String>,
    /// Per-genotype filter values.
    pub ft: Option<Vec<String>>,
    /// Genotype quality.
    pub gq: Option<f64>,
    /// Paired-end coverage.
    pub pec: Option<i32>,
    /// Paired-end variant support.
    pub pev: Option<i32>,
    /// Split-read coverage.
    pub src: Option<i32>,
    /// Split-read variant support.
    pub srv: Option<i32>,
    /// Average mapping quality.
    pub amq: Option<i32>,
    /// Copy number.
    pub cn: Option<i32>,
    /// Average normalized coverage.
    pub anc: Option<f64>,
    /// Point count (windows/targets/probes).
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

            if prev {
                result.push(',');
            }
            if let Some(ft) = &entry.ft {
                prev = true;
                result.push_str(&format!("\"\"\"ft\"\"\":["));
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

            if prev {
                result.push(',');
            }
            if let Some(gq) = &entry.gq {
                prev = true;
                result.push_str(&format!("\"\"\"ad\"\"\":{}", gq));
            }

            if prev {
                result.push(',');
            }
            if let Some(pec) = &entry.pec {
                prev = true;
                result.push_str(&format!("\"\"\"pec\"\"\":{}", pec));
            }

            if prev {
                result.push(',');
            }
            if let Some(pev) = &entry.pev {
                prev = true;
                result.push_str(&format!("\"\"\"pev\"\"\":{}", pev));
            }

            if prev {
                result.push(',');
            }
            if let Some(src) = &entry.src {
                prev = true;
                result.push_str(&format!("\"\"\"src\"\"\":{}", src));
            }

            if prev {
                result.push(',');
            }
            if let Some(srv) = &entry.srv {
                prev = true;
                result.push_str(&format!("\"\"\"srv\"\"\":{}", srv));
            }

            if prev {
                result.push(',');
            }
            if let Some(amq) = &entry.amq {
                prev = true;
                result.push_str(&format!("\"\"\"amq\"\"\":{}", amq));
            }

            if prev {
                result.push(',');
            }
            if let Some(cn) = &entry.cn {
                prev = true;
                result.push_str(&format!("\"\"\"cn\"\"\":{}", cn));
            }

            if prev {
                result.push(',');
            }
            if let Some(anc) = &entry.anc {
                prev = true;
                result.push_str(&format!("\"\"\"anc\"\"\":{}", anc));
            }

            if prev {
                result.push(',');
            }
            if let Some(pc) = &entry.pc {
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

        // if !self.fill_coords(
        //     self.assembly.expect("assembly must have been set"),
        //     record,
        //     &mut tsv_record,
        // )? {
        //     // Record was not on canonical chromosome and should not be written out.
        //     return Ok(());
        // }
        // self.fill_genotype_and_freqs(record, &mut tsv_record)?;
        // self.fill_bg_freqs(record, &mut tsv_record)?;
        // self.fill_clinvar(record, &mut tsv_record)?;
        // self.expand_refseq_ensembl_and_write(record, &mut tsv_record)
        todo!()
    }

    fn set_hgnc_map(&mut self, hgnc_map: FxHashMap<String, HgncRecord>) {
        self.hgnc_map = Some(hgnc_map)
    }

    fn set_assembly(&mut self, assembly: &Assembly) {
        self.assembly = Some(*assembly)
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
            inner: if p.as_ref().extension().unwrap() == "gz" {
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
    // pub info: String,
    /// Genotype call information.
    pub genotype: GenotypeCalls,
}

/// Enumeration for the supported variant callers.
pub enum SvCaller {
    Manta { version: String },
    Delly { version: String },
    Gcnv { version: String },
    DragenSv { version: String },
    DragenCnv { version: String },
    Popdel { version: String },
}

/// Trait that allows the conversion from a VCF record into a `VarFishStrucvarTsvRecord`.
pub trait VcfRecordConverter {
    /// Convert the VCF record into a `VarFishStrucvarTsvRecord`.
    fn convert(
        &self,
        record: &VcfRecord,
        genome_release: GenomeRelease,
    ) -> Result<VarFishStrucvarTsvRecord, anyhow::Error>;

    /// Fill the coordinates of the given record.
    ///
    /// In the case of breakends, the second chromosome and end position will be parsed from
    /// the (only by assumption) alternate allele in `vcf_record`.
    ///
    /// # Preconditions
    ///
    /// - `tsv_record.sv_type` must have been set
    /// - `tsv_record.sv_sub_type` must have been set
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
        // end position.
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
            let bnd = bnd::Breakend::from_ref_alt_str(&reference, &bnd_string)?;

            tsv_record.chromosome2 = bnd.chrom.clone();
            end = Some(bnd.pos);
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
}

/// Conversion from VCF records to `VarFishStrucvarTsvRecord`.
mod conv {
    use std::str::FromStr;

    use crate::common::GenomeRelease;

    use super::SvSubType;
    use super::VarFishStrucvarTsvRecord;
    use super::VcfRecordConverter;
    use noodles::vcf::{
        header::info::key::Key as InfoKey, header::info::key::Standard as InfoKeyStandard,
        record::info::field::value::Value as InfoValue, Record as VcfRecord,
    };

    pub struct MantaVcfRecordConverter {
        pub version: String,
    }

    impl MantaVcfRecordConverter {
        pub fn new(version: &str) -> Self {
            Self {
                version: version.to_string(),
            }
        }
    }

    impl VcfRecordConverter for MantaVcfRecordConverter {
        fn convert(
            &self,
            record: &VcfRecord,
            genome_release: GenomeRelease,
        ) -> Result<VarFishStrucvarTsvRecord, anyhow::Error> {
            let mut result = VarFishStrucvarTsvRecord::default();

            let sv_type = record
                .info()
                .get(&InfoKey::Standard(InfoKeyStandard::SvType))
                .ok_or_else(|| anyhow::anyhow!("SVTYPE not found"))?
                .ok_or_else(|| anyhow::anyhow!("SVTYPE empty"))?;
            if let InfoValue::String(value) = sv_type {
                result.sv_sub_type = SvSubType::from_str(value)?;
            } else {
                anyhow::bail!("SVTYPE is not a string");
            }
            result.sv_type = result.sv_sub_type.into();

            self.fill_coords(&record, genome_release, &mut result)?;
            Ok(result)
        }
    }
}

/// Construct a `VcfRecordConverter` for the given caller.
pub fn build_vcf_record_converter(caller: &SvCaller) -> Box<dyn VcfRecordConverter> {
    match caller {
        SvCaller::Manta { version } => Box::new(conv::MantaVcfRecordConverter::new(version)),
        _ => todo!(),
    }
}

/// Convert the records in the VCF path `path_input` to the JSONL file per contig in `tmp_dir`.
///
/// Note that we will consider the "25 canonical" contigs only (chr1..chr22, chrX, chrY, chrM).
fn run_vcf_to_jsonl(path_input: &str, tmp_dir: &TempDir) -> Result<(), anyhow::Error> {
    todo!();

    Ok(())
}

/// Read through the JSONL file in `tmp_dir` for contig no. `i` and cluster the records.
fn read_and_cluster(
    tmp_dir: &TempDir,
    i: i32,
    args: &Args,
) -> Result<Vec<VarFishStrucvarTsvRecord>, anyhow::Error> {
    todo!()
}

/// Write the clusters to the VCF `writer`.
fn write_for_contig(
    writer: &mut dyn AnnotatedVcfWriter,
    clusters: Vec<VarFishStrucvarTsvRecord>,
) -> Result<(), anyhow::Error> {
    todo!()
}

/// Run the annotation with the given `Write` within the `VcfWriter`.
fn run_with_writer(writer: &mut dyn AnnotatedVcfWriter, args: &Args) -> Result<(), anyhow::Error> {
    // Create temporary directory.  We will create one temporary file (containing `jsonl`
    // seriealized `VarFishStrucvarTsvRecord`s) for each SV type and contig.
    let tmp_dir = TempDir::new("mehari")?;

    // Read through input VCF files and write out to temporary files.
    for path_input in &args.path_input_vcf {
        run_vcf_to_jsonl(path_input, &tmp_dir)?;
    }

    // Generate output header and write to `writer`.
    let header_out = build_header();
    writer.write_header(&header_out)?;

    // Read through temporary files by contig, cluster by overlap as configured, and write to `writer`.
    for i in 1..=25 {
        let clusters = read_and_cluster(&tmp_dir, i, args)?;
        write_for_contig(writer, clusters)?;
    }

    // tracing::info!("Open VCF and read header");
    // let mut reader = VariantReaderBuilder::default().build_from_path(&args.path_input_vcf)?;
    // let header_in = reader.read_header()?;
    // let header_out = build_header(&header_in);

    // // Guess genome release from contigs in VCF header.
    // let genome_release = args.genome_release.map(|gr| match gr {
    //     GenomeRelease::Grch37 => Assembly::Grch37p10, // has chrMT!
    //     GenomeRelease::Grch38 => Assembly::Grch38,
    // });
    // let assembly = guess_assembly(&header_in, false, genome_release)?;
    // writer.set_assembly(&assembly);
    // tracing::info!("Determined input assembly to be {:?}", &assembly);

    // // Open the serialized transcripts.
    // tracing::info!("Opening transcript database");
    // let tx_db = load_tx_db(
    //     &format!(
    //         "{}/seqvars/{}/txs.bin",
    //         &args.path_db,
    //         path_component(assembly)
    //     ),
    //     args.max_fb_tables,
    // )?;
    // tracing::info!("Building transcript interval trees ...");
    // let provider = Rc::new(MehariProvider::new(tx_db, assembly));
    // tracing::info!("... done building transcript interval trees");

    // // Perform the VCF annotation (mostly merging for SVs).
    // tracing::info!("Annotating VCF ...");
    // let start = Instant::now();
    // let mut prev = Instant::now();
    // let mut total_written = 0usize;

    // writer.write_header(&header_out)?;
    // let mut records = reader.records(&header_in);
    // loop {
    //     if let Some(rec,ord) = records.next() {
    //         // TODO
    //     } else {
    //         break; // all done
    //     }

    //     total_written += 1;
    //     if let Some(max_var_count) = args.max_var_count {
    //         if total_written >= max_var_count {
    //             tracing::warn!(
    //                 "Stopping after {} records as requested by --max-var-count",
    //                 total_written
    //             );
    //             break;
    //         }
    //     }
    // }
    // tracing::info!(
    //     "... annotated {} records in {:?}",
    //     total_written.separate_with_commas(),
    //     start.elapsed()
    // );

    Ok(())
}

/// Main entry point for `annotate strucvars` sub command.
pub fn run(_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("config = {:#?}", &args);
    if let Some(path_output_vcf) = &args.output.path_output_vcf {
        if path_output_vcf.ends_with(".vcf.gz") || path_output_vcf.ends_with(".vcf.bgzf") {
            let mut writer = VcfWriter::new(
                File::create(path_output_vcf)
                    .map(BufWriter::new)
                    .map(BgzfWriter::new)?,
            );
            run_with_writer(&mut writer, args)?;
        } else {
            let mut writer = VcfWriter::new(File::create(path_output_vcf).map(BufWriter::new)?);
            run_with_writer(&mut writer, args)?;
        }
    } else {
        // Load the HGNC xlink map.
        let hgnc_map = {
            tracing::info!("Loading HGNC map ...");
            let mut result = FxHashMap::default();

            let tsv_file = File::open(&format!("{}/hgnc.tsv", &args.path_db,))?;
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
        let mut writer = VarFishStrucvarTsvWriter::with_path(path_output_tsv);

        // Load the pedigree.
        tracing::info!("Loading pedigree...");
        writer.set_pedigree(&PedigreeByName::from_path(
            &args.path_input_ped.as_ref().unwrap(),
        )?);
        tracing::info!("... done loading pedigree");

        writer.set_hgnc_map(hgnc_map);
        run_with_writer(&mut writer, args)?;
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
            let mut chrom_pos = chrom_pos.split(":");
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

    use noodles::vcf;
    use pretty_assertions::assert_eq;
    use temp_testdir::TempDir;

    use crate::{annotate::strucvars::PeOrientation, common::GenomeRelease};

    use super::{bnd::Breakend, conv::MantaVcfRecordConverter, VcfRecordConverter};

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

    /// Conversion of minimal Manta VCF to JSONL.
    #[test]
    fn vcf_to_jsonl_manta_min() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let out_jsonl = File::create(temp.join("out.jsonl"))?;

        let mut reader = vcf::reader::Builder::default()
            .build_from_path("tests/data/annotate/strucvars/manta-min.vcf")?;
        let header_in = reader.read_header()?.parse()?;

        let converter = MantaVcfRecordConverter::new("1.6.0");
        let mut records = reader.records(&header_in);
        loop {
            if let Some(record) = records.next() {
                let record = converter.convert(&record?, GenomeRelease::Grch37)?;
                jsonl::write(&out_jsonl, &record)?;
            } else {
                break; // all done
            }
        }
        drop(out_jsonl);

        let expected = std::fs::read_to_string("tests/data/annotate/strucvars/manta-min.jsonl")?;
        let actual = std::fs::read_to_string(temp.join("out.jsonl"))?;

        assert_eq!(expected, actual);

        Ok(())
    }
}
