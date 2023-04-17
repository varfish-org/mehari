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
use noodles::vcf::{Header as VcfHeader, Record as VcfRecord, Writer as VcfWriter};
use rustc_hash::FxHashMap;
use strum::{Display, EnumIter};
use uuid::Uuid;

use super::seqvars::AnnotatedVcfWriter;

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
    /// Path to the input VCF file.
    #[arg(long)]
    pub path_input_vcf: String,
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
#[derive(Debug, Default)]
struct GenotypeInfo {
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
            result.push_str(&format!("\"\"\"{}\"\"\":{", entry.name));

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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Display, Default)]
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
#[derive(EnumIter, PartialEq, Eq, Ord, PartialOrd, Hash, Debug, Clone, Copy, Default, Display)]
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
#[derive(PartialEq, Debug, Clone, Copy, Default, Display)]
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
#[derive(Debug, Default)]
pub struct VarFishStrucvarTsvRecord {
    pub release: String,
    pub chromosome: String,
    pub chromosome_no: u32,
    pub bin: u32,
    pub chromosome2: String,
    pub chromosome_no2: u32,
    pub bin2: u32,
    pub pe_orientation: PeOrientation,

    pub start: usize,
    pub end: usize,
    pub start_ci_left: i32,
    pub start_ci_right: i32,
    pub end_ci_left: i32,
    pub end_ci_right: i32,

    pub case_id: usize,
    pub set_id: usize,
    pub sv_uuid: Uuid,
    pub callers: Vec<String>,
    pub sv_type: SvType,
    pub sv_sub_type: SvSubType,

    // pub info: String,
    pub genotype: GenotypeCalls,
}

/// Run the annotation with the given `Write` within the `VcfWriter`.
fn run_with_writer(writer: &mut dyn AnnotatedVcfWriter, args: &Args) -> Result<(), anyhow::Error> {
    todo!()
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

#[cfg(test)]
mod test {}
