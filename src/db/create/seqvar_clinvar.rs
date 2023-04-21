//! Creation of mehari internal sequence variant ClinVar database.

use std::{fmt::Display, str::FromStr};
use std::{
    io::{BufRead, BufReader},
    time::Instant,
};

use bgzip::BGZFReader;
use clap::Parser;
use hgvs::static_data::Assembly;
use rocksdb::{DBWithThreadMode, SingleThreaded};
use serde::{Deserialize, Serialize};
use thousands::Separable;

use crate::common::GenomeRelease;

use super::seqvar_freqs::serialized::vcf::Var;

/// Enumeration for ClinVar pathogenicity.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Pathogenicity {
    /// Pathogenic.
    Pathogenic,
    /// Likely pathogenic.
    LikelyPathogenic,
    /// Uncertain significance.
    UncertainSignificance,
    /// Likely benign.
    LikelyBenign,
    /// Benign.
    Benign,
}

impl Display for Pathogenicity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Pathogenicity::Pathogenic => write!(f, "pathogenic"),
            Pathogenicity::LikelyPathogenic => write!(f, "likely pathogenic"),
            Pathogenicity::UncertainSignificance => write!(f, "uncertain significance"),
            Pathogenicity::LikelyBenign => write!(f, "likely benign"),
            Pathogenicity::Benign => write!(f, "benign"),
        }
    }
}

impl FromStr for Pathogenicity {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "pathogenic" => Ok(Pathogenicity::Pathogenic),
            "likely pathogenic" => Ok(Pathogenicity::LikelyPathogenic),
            "uncertain significance" => Ok(Pathogenicity::UncertainSignificance),
            "likely benign" => Ok(Pathogenicity::LikelyBenign),
            "benign" => Ok(Pathogenicity::Benign),
            _ => anyhow::bail!("Unknown pathogenicity: {}", s),
        }
    }
}

impl Serialize for Pathogenicity {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for Pathogenicity {
    fn deserialize<D>(d: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(d)?;
        Self::from_str(&s).map_err(serde::de::Error::custom)
    }
}

/// Reading `clinvar-tsv` sequence variant files.
pub mod reading {
    use serde::{Deserialize, Serialize};

    use super::Pathogenicity;

    /// Representation of a record from the `clinvar-tsv` output.
    ///
    /// Note that the pathogenicity and review status are available in two fashions.  The first is
    /// "ClinVar style" and attempts to follow the ClinVar approach.  Here, variant assessments
    /// with a higher star rating override those a lowe rone.  This is what most users want.
    /// The assessment "paranoid" uses all assessments, including those without a star rating,
    /// on the same level.
    #[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Eq)]
    pub struct Record {
        /// Genome release.
        pub release: String,
        /// Chromosome name.
        pub chromosome: String,
        /// 1-based start position.
        pub start: u32,
        /// 1-based end position.
        pub end: u32,
        /// Reference allele bases in VCF notation.
        pub reference: String,
        /// Alternative allele bases in VCF notation.
        pub alternative: String,
        /// VCV accession identifier.
        pub vcv: String,
        /// Pathogenicity summary for the variant (ClinVar style).
        #[serde(deserialize_with = "deserialize_pathogenicity")]
        pub summary_clinvar_pathogenicity: Vec<Pathogenicity>,
        /// Pathogenicity gold stars (ClinVar style).
        pub summary_clinvar_gold_stars: u32,
        /// Pathogenicity summary for the variant ("paranoid" style).
        #[serde(deserialize_with = "deserialize_pathogenicity")]
        pub summary_paranoid_pathogenicity: Vec<Pathogenicity>,
        /// Pathogenicity gold stars ("paranoid" style).
        pub summary_paranoid_gold_stars: u32,
    }

    fn deserialize_pathogenicity<'de, D>(deserializer: D) -> Result<Vec<Pathogenicity>, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s: &str = Deserialize::deserialize(deserializer)?;
        let s = s.replace('{', "[").replace('}', "]");
        serde_json::from_str(&s).map_err(serde::de::Error::custom)
    }
}

pub mod serialize {
    use serde::{Deserialize, Serialize};

    use super::Pathogenicity;

    #[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Eq)]
    pub struct Record {
        /// Genome release.
        pub release: String,
        /// Chromosome name.
        pub chromosome: String,
        /// 1-based start position.
        pub start: u32,
        /// 1-based end position.
        pub end: u32,
        /// Reference allele bases in VCF notation.
        pub reference: String,
        /// Alternative allele bases in VCF notation.
        pub alternative: String,
        /// VCV accession identifier.
        pub vcv: String,
        /// Pathogenicity summary for the variant (ClinVar style).
        pub summary_clinvar_pathogenicity: Vec<Pathogenicity>,
        /// Pathogenicity gold stars (ClinVar style).
        pub summary_clinvar_gold_stars: u32,
        /// Pathogenicity summary for the variant ("paranoid" style).
        pub summary_paranoid_pathogenicity: Vec<Pathogenicity>,
        /// Pathogenicity gold stars ("paranoid" style).
        pub summary_paranoid_gold_stars: u32,
    }

    impl From<super::reading::Record> for Record {
        fn from(r: super::reading::Record) -> Self {
            Self {
                release: r.release,
                chromosome: r.chromosome,
                start: r.start,
                end: r.end,
                reference: r.reference,
                alternative: r.alternative,
                vcv: r.vcv,
                summary_clinvar_pathogenicity: r.summary_clinvar_pathogenicity,
                summary_clinvar_gold_stars: r.summary_clinvar_gold_stars,
                summary_paranoid_pathogenicity: r.summary_paranoid_pathogenicity,
                summary_paranoid_gold_stars: r.summary_paranoid_gold_stars,
            }
        }
    }
}

#[derive(Parser, Debug)]
#[command(about = "Construct mehari sequence variant ClinVar database", long_about = None)]
pub struct Args {
    /// Genome release to use, default is to auto-detect.
    #[arg(long, value_enum)]
    pub genome_release: Option<GenomeRelease>,
    /// Path to the output database to build.
    #[arg(long)]
    pub path_output_db: String,

    /// Path to the sequence variant ClinVar TSV file from `clinvar-tsv`.
    #[arg(long)]
    pub path_clinvar_tsv: String,

    /// For debug purposes, maximal number of variants to import.
    #[arg(long)]
    pub max_var_count: Option<usize>,
}

/// Convert release string to Assembly.
fn str_to_assembly(s: &str) -> Result<Assembly, anyhow::Error> {
    match s {
        "GRCh37" => Ok(Assembly::Grch37p10),
        "GRCh38" => Ok(Assembly::Grch38),
        _ => anyhow::bail!("Unknown genome release: {}", s),
    }
}

/// Guess genome release from TSV file.
fn guess_genome_release(reader: Box<dyn BufRead>) -> Result<Assembly, anyhow::Error> {
    let reader = std::io::BufReader::new(reader);
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(reader);
    let record: reading::Record = reader.deserialize().next().unwrap()?;
    str_to_assembly(record.release.as_str())
}

fn open_tsv(args: &Args) -> Result<Box<dyn BufRead>, anyhow::Error> {
    Ok(if args.path_clinvar_tsv.ends_with(".gz") {
        Box::new(BGZFReader::new(std::fs::File::open(
            &args.path_clinvar_tsv,
        )?)?)
    } else {
        Box::new(BufReader::new(std::fs::File::open(&args.path_clinvar_tsv)?))
    })
}

/// Import ClinVar TSV file into RocksDB column family.
fn import_clinvar_seqvars(
    args: &Args,
    genome_release: Option<Assembly>,
    db: &DBWithThreadMode<SingleThreaded>,
    cf_clinvar: &rocksdb::ColumnFamily,
) -> Result<(), anyhow::Error> {
    let start = Instant::now();

    // Read first record from TSV file to guess genome release.
    let reader = open_tsv(args)?;
    let guessed_release = guess_genome_release(reader)?;
    // Check genome release guessed vs the one passed as argument.
    let genome_release = if let Some(genome_release) = genome_release {
        if genome_release != guessed_release {
            anyhow::bail!(
                "Genome release mismatch: {:?} vs. {:?}",
                genome_release,
                guessed_release
            );
        } else {
            genome_release
        }
    } else {
        guessed_release
    };

    // Open TSV file again and read all records, writing them to db/cf_clinvar.
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(open_tsv(args)?);
    let mut prev = Instant::now();
    let mut records_written = 0;
    for record in reader.deserialize() {
        let mut record: reading::Record = record?;

        record.summary_clinvar_pathogenicity.sort();
        record.summary_paranoid_pathogenicity.sort();

        let var = Var {
            chrom: record.chromosome.clone(),
            pos: record.start,
            reference: record.reference.clone(),
            alternative: record.alternative.clone(),
        };
        tracing::trace!("now at {:?}", &var);
        if var.pos == 865567 {
            tracing::info!("at {:?}", &var);
        }

        if prev.elapsed().as_secs() >= 60 {
            tracing::info!("at {:?}", &var);
            prev = Instant::now();
        }

        // Validate record's genome release.
        let record_release = str_to_assembly(record.release.as_str())?;
        if record_release != genome_release {
            anyhow::bail!(
                "Genome release mismatch: {:?} vs. {:?}",
                genome_release,
                record_release
            );
        }

        // Serialize record to `Vec<u8>` using `bincode`.
        let record: serialize::Record = record.into();
        let value = bincode::serialize(&record)?;
        // Derive key from record.
        let key: Vec<u8> = var.into();

        db.put_cf(cf_clinvar, key, value)?;

        records_written += 1;

        if let Some(max_var_count) = args.max_var_count {
            if records_written >= max_var_count {
                tracing::warn!("Stopping after {} records as requested", max_var_count);
                break;
            }
        }
    }

    tracing::info!(
        "  wrote {} ClinVar seqvars records in {:?}",
        records_written.separate_with_commas(),
        start.elapsed()
    );

    Ok(())
}

pub fn rocksdb_tuning(options: rocksdb::Options) -> rocksdb::Options {
    let mut options = options;

    // compress all files with Zstandard
    options.set_compression_per_level(&[]);
    options.set_compression_type(rocksdb::DBCompressionType::Zstd);
    // We only want to set level to 2 but have to set the rest as well using the Rust interface.
    // The (default) values for the other levels were taken from the output of a RocksDB
    // output folder created with default settings.
    options.set_compression_options(-14, 2, 0, 0);
    // options.set_zstd_max_train_bytes(100 * 1024);

    options
}

/// Main entry point for `db create seqvar-clinvar` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!(
        "Building sequence variant ClinVar database\ncommon args: {:#?}\nargs: {:#?}",
        common,
        args
    );

    // Guess genome release from paths.
    let genome_release = args.genome_release.map(|gr| match gr {
        GenomeRelease::Grch37 => Assembly::Grch37p10, // has chrMT!
        GenomeRelease::Grch38 => Assembly::Grch38,
    });

    // Open the output database, obtain column family handles, and write out meta data.
    tracing::info!("Opening output database");
    let mut options = rocksdb_tuning(rocksdb::Options::default());
    options.create_if_missing(true);
    options.create_missing_column_families(true);
    options.prepare_for_bulk_load();
    options.set_disable_auto_compactions(true);
    let db = rocksdb::DB::open_cf(&options, &args.path_output_db, ["meta", "clinvar_seqvars"])?;

    let cf_meta = db.cf_handle("meta").unwrap();
    let cf_clinvar = db.cf_handle("clinvar_seqvars").unwrap();

    tracing::info!("Writing meta data to database");
    db.put_cf(cf_meta, "genome-release", format!("{genome_release:?}"))?;

    // Import the ClinVar data TSV file.
    tracing::info!("Processing ClinVar TSV ...");
    import_clinvar_seqvars(args, genome_release, &db, cf_clinvar)?;

    // Finally, compact manually.
    tracing::info!("Enforcing manual compaction");
    db.compact_range_cf(cf_meta, None::<&[u8]>, None::<&[u8]>);
    db.compact_range_cf(cf_clinvar, None::<&[u8]>, None::<&[u8]>);

    let compaction_start = Instant::now();
    let mut last_printed = compaction_start;
    while db
        .property_int_value(rocksdb::properties::COMPACTION_PENDING)?
        .unwrap()
        > 0
        || db
            .property_int_value(rocksdb::properties::NUM_RUNNING_COMPACTIONS)?
            .unwrap()
            > 0
    {
        std::thread::sleep(std::time::Duration::from_millis(100));
        if last_printed.elapsed() > std::time::Duration::from_millis(1000) {
            log::info!(
                "... waiting for compaction for {:?}",
                compaction_start.elapsed()
            );
            last_printed = Instant::now();
        }
    }

    tracing::info!("Done building sequence variant ClinVar database");
    Ok(())
}
