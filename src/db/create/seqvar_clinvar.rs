//! Creation of mehari internal sequence variant ClinVar database.

use std::time::Instant;

use clap::Parser;
use hgvs::static_data::Assembly;

use crate::common::GenomeRelease;

/// Reading `clinvar-tsv` sequence variant files.
pub mod reading {
    use serde::{Serialize, Deserialize};

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
        release: String,
        /// Chromosome name.
        chromosome: String,
        /// 1-based start position.
        start: i32,
        /// 1-based end position.
        end: i32,
        /// Reference allele bases in VCF notation.
        reference: String,
        /// Alternative allele bases in VCF notation.
        alternative: String,
        /// VCV accession identifier.
        vcv: String,
        /// Pathogenicity summary for the variant (ClinVar style).
        #[serde(rename = "summary_clinvar_pathogenicity_label")]
        summary_clinvar_pathogenicity: String,
        /// Pathogenicity review status for the variant (ClinVar style).
        summary_clinvar_review_status: String,
        /// Pathogenicity gold stars (ClinVar style).
        summary_clinvar_gold_stars: i32,
        /// Pathogenicity summary for the variant ("paranoid" style).
        #[serde(rename = "summary_paranoid_pathogenicity_label")]
        summary_paranoid_pathogenicity: String,
        /// Pathogenicity review status for the variant ("paranoid" style).
        summary_paranoid_review_status: String,
        /// Pathogenicity gold stars ("paranoid" style).
        summary_paranoid_gold_stars: i32,
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
    let mut options = rocksdb::Options::default();
    options.create_if_missing(true);
    options.create_missing_column_families(true);
    options.prepare_for_bulk_load();
    options.set_disable_auto_compactions(true);
    let db = rocksdb::DB::open_cf(
        &options,
        &args.path_output_db,
        ["meta", "clinvar_seqvars"],
    )?;

    let cf_meta = db.cf_handle("meta").unwrap();
    let cf_clinvar = db.cf_handle("clinvar_seqvars").unwrap();

    tracing::info!("Writing meta data to database");
    db.put_cf(cf_meta, "genome-release", format!("{genome_release:?}"))?;

    // Import the ClinVar data TSV file.
    tracing::info!("Processing ClinVar TSV ...");
    todo!();

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
