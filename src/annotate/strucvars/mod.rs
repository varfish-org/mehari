//! Annotation of structural variant VCF files.

use std::io::Write;
use std::path::Path;
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
    inner: Box<dyn Write>,
    assembly: Option<Assembly>,
    pedigree: Option<PedigreeByName>,
    header: Option<VcfHeader>,
    hgnc_map: Option<FxHashMap<String, HgncRecord>>,
}

/// TODO XXX
#[derive(Debug, Default)]
struct GenotypeInfo {
    pub name: String,
    pub gt: Option<String>,
    // XXX
}

#[derive(Debug, Default)]
struct GenotypeCalls {
    pub entries: Vec<GenotypeInfo>,
}

impl GenotypeCalls {
    pub fn for_tsv(&self) -> String {
        todo!()
    }
}

/// Implement `AnnotatedVcfWriter` for `VarFishTsvWriter`.
impl AnnotatedVcfWriter for VarFishStrucvarTsvWriter {
    fn write_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error> {
        todo!()
    }

    fn write_record(&mut self, record: &VcfRecord) -> Result<(), anyhow::Error> {
        todo!()
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

            // Load the pedigree.
            tracing::info!("Loading pedigree...");
            writer.set_pedigree(&PedigreeByName::from_path(
                &args.path_input_ped.as_ref().unwrap(),
            )?);
            tracing::info!("... done loading pedigree");

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
        writer.set_hgnc_map(hgnc_map);
        run_with_writer(&mut writer, args)?;
    }

    Ok(())
}
