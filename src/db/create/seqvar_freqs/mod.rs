//! Population frequencies for sequence variants.

pub mod reading;
pub mod serialized;

use std::{fs::File, io::BufReader};

use clap::Parser;
use noodles::{bgzf, vcf};

use serialized::*;

/// Select the genome release to use.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum GenomeRelease {
    Grch37,
    Grch38,
}

/// Command line arguments for `db create seqvar-freqs` sub command.
#[derive(Parser, Debug)]
#[command(about = "Construct mehari sequence variant frequencies database", long_about = None)]
pub struct Args {
    /// Genome release to use, default is to auto-detect.
    #[arg(long, value_enum)]
    pub genome_release: Option<GenomeRelease>,
    /// Path to the output database to build.
    #[arg(long)]
    pub path_output_db: String,

    /// Path(s) to the gnomAD exomes VCF file(s).
    #[arg(long)]
    pub path_gnomad_exomes: Option<Vec<String>>,
    /// Path(s) to the gnomAD genomes VCF file(s).
    #[arg(long)]
    pub path_gnomad_genomes: Option<Vec<String>>,
    /// Path(s) to the gnomAD mtDNA VCF file(s).
    #[arg(long)]
    pub path_gnomad_mtdna: Option<String>,
    /// Path(s) to the HelixMtDb TSV file.
    #[arg(long)]
    pub path_helix_mtdb: Option<String>,

    /// For debug purposes, maximal number of variants to import.
    #[arg(long)]
    pub max_var_count: Option<usize>,
}

/// Helper for reading through gnomAD mtDNA and HelixMtDb data;
struct MtReader {
    /// CSV reader for the gnomad mitochondrial records.
    gnomad_reader: Option<vcf::Reader<BufReader<bgzf::reader::Reader<File>>>>,
    /// CSV reader for the HelixMtDb records.
    helix_reader: Option<csv::Reader<BufReader<File>>>,
}

impl MtReader {
    pub fn new(
        path_gnomad: &Option<String>,
        path_helix: &Option<String>,
    ) -> Result<Self, anyhow::Error> {
        let gnomad_reader = if let Some(path_gnomad) = path_gnomad {
            tracing::info!("Opening gnomAD chrMT file {}", &path_gnomad);
            let bgzf = bgzf::reader::Builder::default().build_from_path(path_gnomad)?;
            let buf_reader = BufReader::new(bgzf);
            let vcf_reader = vcf::Reader::new(buf_reader);

            Some(vcf_reader)
        } else {
            None
        };

        let helix_reader = if let Some(path_helix) = path_helix {
            tracing::info!("Opening HelixMtDb chrMT file {}", &path_helix);
            let f = File::open(path_helix)?;
            let buf_reader = BufReader::new(f);
            let csv_reader = csv::ReaderBuilder::new()
                .has_headers(true)
                .delimiter(b'\t')
                .from_reader(buf_reader);

            Some(csv_reader)
        } else {
            None
        };

        Ok(Self {
            gnomad_reader,
            helix_reader,
        })
    }

    pub fn run<F>(&mut self, mut func: F) -> Result<(), anyhow::Error>
    where
        F: FnMut(Option<(VcfVar, MtCounts)>, Option<(VcfVar, MtCounts)>),
    {
        // // Prepare VCF.
        // let vcf_header = self
        //     .gnomad_reader
        //     .as_mut()
        //     .map(|gnomad_reader| {
        //         Rc::new(
        //             gnomad_reader
        //                 .read_header()
        //                 .expect("could not read vcf header")
        //                 .parse::<vcf::Header>()
        //                 .expect("could not parse VCF header"),
        //         )
        //     })
        //     .unwrap_or(Rc::default());
        // let vcf_records = self
        //     .gnomad_reader
        //     .as_mut()
        //     .map(|gnomad_reader| gnomad_reader.records(vcf_header.as_ref()));

        // // Prepare TSV.
        // let csv_xs: Option<csv::DeserializeRecordsIter<_, HelixRecord>> = self
        //     .helix_reader
        //     .as_mut()
        //     .map(|helix_reader| helix_reader.deserialize());

        // // Read initial "next" records into buffers.
        // let mut next_gnomad = vcf_records.map(|mut vcf_records| -> Option<_> {
        //     if let Some(Ok(vcf_record)) = vcf_records.next() {
        //         let chrom = match vcf_record.chromosome() {
        //             vcf::record::Chromosome::Name(s) | vcf::record::Chromosome::Symbol(s) => {
        //                 s.clone()
        //             }
        //         };
        //         let pos: usize = vcf_record.position().into();
        //         let reference = vcf_record.reference_bases().to_string();
        //         let alternative = vcf_record.alternate_bases().deref()[0].to_string();

        //         let ac_hom = match vcf_record
        //             .info()
        //             .get(&Key::Other("AC_hom".to_string()))
        //             .unwrap()
        //             .unwrap()
        //         {
        //             vcf::record::info::field::Value::Integer(ac_hom) => *ac_hom as u32,
        //             _ => panic!("invalid type for AC_hom"),
        //         };
        //         let ac_het = match vcf_record
        //             .info()
        //             .get(&Key::Other("AC_het".to_string()))
        //             .unwrap()
        //             .unwrap()
        //         {
        //             vcf::record::info::field::Value::Integer(ac_het) => *ac_het as u32,
        //             _ => panic!("invalid type for AC_het"),
        //         };
        //         let an = match vcf_record
        //             .info()
        //             .get(&Key::TotalAlleleCount)
        //             .unwrap()
        //             .unwrap()
        //         {
        //             vcf::record::info::field::Value::Integer(an) => *an as u32,
        //             _ => panic!("invalid type for AN"),
        //         };

        //         let var = VcfVar {
        //             chrom,
        //             pos: pos as i32,
        //             reference,
        //             alternative,
        //         };
        //         let mt_counts = MtCounts { ac_hom, ac_het, an };
        //         tracing::info!("Read VCF variant {:?} -- {:?}", &var, &mt_counts);

        //         Some(0)
        //     } else {
        //         None
        //     }
        // });

        todo!()
    }
}

/// Check genome release from gnomAD and compare with the provided one from arguments.
fn check_release_with_guessed_from_gnomad(args: &Args) -> Result<GenomeRelease, anyhow::Error> {
    let release = if args.path_gnomad_exomes.is_some() || args.path_gnomad_genomes.is_some() {
        tracing::info!("Genome release evaluation for gnomAD files...");
        let mut paths = args
            .path_gnomad_exomes
            .as_ref()
            .map(|s| s.clone())
            .unwrap_or(Vec::new());
        paths.append(
            &mut args
                .path_gnomad_genomes
                .as_ref()
                .map(|s| s.clone())
                .unwrap_or(Vec::new()),
        );

        let mut current = None;
        let mut current_from = String::new();
        // for path in &paths {
        //     let guessed = gnomad_nuclear::guess_genome_release(path)?;
        //     if let Some(current) = &current {
        //         if current != &guessed {
        //             let msg = format!(
        //                 "Inconsistent gnomAD releases. Earlier, we guessed {:?} from {} but \
        //                 now we guessed {:?} from {}",
        //                 current, &current_from, guessed, &path
        //             );
        //             tracing::error!(msg);
        //             return Err(anyhow::anyhow!(msg));
        //         }
        //     } else {
        //         current = Some(guessed);
        //         current_from = path.clone();
        //     }
        // }

        current
    } else {
        args.genome_release
    };

    if let Some(release) = release {
        Ok(release)
    } else {
        let msg = "Could not determine genome release from gnomAD and none given on command line";
        tracing::error!(msg);
        Err(anyhow::anyhow!(msg))
    }
}

/// Main entry point for `db create seqvar_freqs` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!(
        "Building sequence variant frequencies table\ncommon args: {:#?}\nargs: {:#?}",
        common,
        args
    );

    // Guess genome release from paths.
    let release = check_release_with_guessed_from_gnomad(args)?;

    // Open the output database, obtain column family handles, and write out meta data.
    tracing::info!("Opening output database");
    let mut options = rocksdb::Options::default();
    options.create_if_missing(true);
    options.create_missing_column_families(true);
    options.prepare_for_bulk_load();
    let mut db = rocksdb::DB::open_cf(
        &options,
        &args.path_output_db,
        &["meta", "nuclear", "mitochondrial"],
    )?;

    let cf_meta = db.cf_handle("meta").unwrap();
    let cf_ndna = db.cf_handle("nuclear").unwrap();
    let cf_mtdna = db.cf_handle("mitochondrial").unwrap();

    tracing::info!("Writing meta data to database");
    db.put_cf(cf_meta, "genome-release", format!("{:?}", release))?;

    // Import gnomAD variants in a chromosome-wise fashion.
    tracing::info!("Processing gnomAD nuclear variant data ...");
    tracing::info!("Opening gnomAD exomes file(s)");
    tracing::info!("Opening gnomAD genomes file(s)");

    // Import chrMT variants.
    tracing::info!("Processing chrMT data ...");
    MtReader::new(&args.path_gnomad_mtdna, &args.path_helix_mtdb)?
        .run(|gnomad: Option<(VcfVar, MtCounts)>, helix: Option<(VcfVar, MtCounts)>| {})?;

    tracing::info!("Opening gnomAD mtDNA file(s)");
    tracing::info!("Opening HelixMtDb file(s)");

    // Finally, compact manually.
    tracing::info!("Enforcing manual compaction");
    db.compact_range_cf(cf_ndna, None::<&[u8]>, None::<&[u8]>);
    db.compact_range_cf(cf_mtdna, None::<&[u8]>, None::<&[u8]>);

    tracing::info!("Done building sequence variant frequency table");
    Ok(())
}
