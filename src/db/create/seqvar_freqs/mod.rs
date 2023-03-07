//! Population frequencies for sequence variants.

pub mod reading;
pub mod serialized;

use clap::Parser;
use hgvs::static_data::Assembly;
use noodles::vcf::Record as VcfRecord;

use serialized::*;

use self::reading::MultiVcfReader;

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
pub struct MtReader {
    /// CSV reader for the gnomAD mitochondrial records.
    gnomad_reader: Option<MultiVcfReader>,
    /// Next variant from gnomAD.
    gnomad_next: Option<VcfRecord>,
    /// CSV reader for the HelixMtDb records.
    helix_reader: Option<MultiVcfReader>,
    /// Next variant from gnomAD.
    helix_next: Option<VcfRecord>,
}

impl MtReader {
    /// Construct new reader with optional paths to HelixMtDb and gnomAD mtDNA.
    ///
    /// Optionally, you can provide an assembly to validate the VCF contigs against.
    pub fn new(
        path_gnomad: Option<&str>,
        path_helix: Option<&str>,
        assembly: Option<Assembly>,
    ) -> Result<Self, anyhow::Error> {
        let mut gnomad_reader = path_gnomad
            .as_ref()
            .map(|path_gnomad| {
                tracing::info!("Opening gnomAD chrMT file {}", &path_gnomad);
                MultiVcfReader::new(&[path_gnomad], assembly)
            })
            .transpose()?;
        let mut helix_reader = path_helix
            .as_ref()
            .map(|path_helix| {
                tracing::info!("Opening HelixMtDb chrMT file {}", &path_helix);
                MultiVcfReader::new(&[path_helix], assembly)
            })
            .transpose()?;

        let gnomad_next = if let Some(gnomad_reader) = gnomad_reader.as_mut() {
            gnomad_reader.pop()?.0
        } else {
            None
        };
        let helix_next = if let Some(helix_reader) = helix_reader.as_mut() {
            helix_reader.pop()?.0
        } else {
            None
        };

        tracing::info!("... done opening and popping");

        Ok(Self {
            gnomad_reader,
            helix_reader,
            gnomad_next,
            helix_next,
        })
    }

    /// Run the reading of the chrMT frequencies.
    ///
    /// Returns whether there is a next value.
    pub fn run<F>(&mut self, mut func: F) -> Result<bool, anyhow::Error>
    where
        F: FnMut(VcfVar, MtCounts, MtCounts) -> Result<(), anyhow::Error>,
    {
        match (&self.gnomad_next, &self.helix_next) {
            (None, Some(helix)) => {
                func(
                    VcfVar::from_vcf(helix),
                    MtCounts::default(),
                    MtCounts::from_vcf(helix),
                )?;
                self.helix_next = self.helix_reader.as_mut().unwrap().pop()?.0;
            }
            (Some(gnomad), None) => {
                func(
                    VcfVar::from_vcf(gnomad),
                    MtCounts::from_vcf(gnomad),
                    MtCounts::default(),
                )?;
                self.gnomad_next = self.gnomad_reader.as_mut().unwrap().pop()?.0;
            }
            (Some(gnomad), Some(helix)) => {
                let var_gnomad = VcfVar::from_vcf(gnomad);
                let var_helix = VcfVar::from_vcf(helix);
                match var_gnomad.cmp(&var_helix) {
                    std::cmp::Ordering::Less => {
                        func(var_gnomad, MtCounts::from_vcf(gnomad), MtCounts::default())?;
                        self.gnomad_next = self.gnomad_reader.as_mut().unwrap().pop()?.0;
                    }
                    std::cmp::Ordering::Equal => {
                        func(
                            var_gnomad,
                            MtCounts::from_vcf(gnomad),
                            MtCounts::from_vcf(helix),
                        )?;
                        self.helix_next = self.helix_reader.as_mut().unwrap().pop()?.0;
                        self.gnomad_next = self.gnomad_reader.as_mut().unwrap().pop()?.0;
                    }
                    std::cmp::Ordering::Greater => {
                        func(var_helix, MtCounts::default(), MtCounts::from_vcf(helix))?;
                        self.helix_next = self.helix_reader.as_mut().unwrap().pop()?.0;
                    }
                }
            }
            (None, None) => (),
        }

        Ok(self.gnomad_next.is_some() || self.helix_next.is_some())
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
    let db = rocksdb::DB::open_cf(
        &options,
        &args.path_output_db,
        ["meta", "nuclear", "mitochondrial"],
    )?;

    let _cf_meta = db.cf_handle("meta").unwrap();
    let cf_ndna = db.cf_handle("nuclear").unwrap();
    let cf_mtdna = db.cf_handle("mitochondrial").unwrap();

    tracing::info!("Writing meta data to database");
    // db.put_cf(cf_meta, "genome-release", format!("{:?}", release))?;

    // Import gnomAD variants in a chromosome-wise fashion.
    tracing::info!("Processing gnomAD nuclear variant data ...");
    tracing::info!("Opening gnomAD exomes file(s)");
    tracing::info!("Opening gnomAD genomes file(s)");

    // Import chrMT variants.
    tracing::info!("Processing chrMT data ...");
    let mut chrmt_written = 0usize;
    let mut mt_reader = MtReader::new(
        args.path_gnomad_mtdna.as_deref(),
        args.path_helix_mtdb.as_deref(),
        genome_release,
    )?;
    let mut has_next = true;
    while has_next {
        has_next = mt_reader.run(|variant, gnomad_mtdna, helix_mtdb| {
            tracing::trace!(
                "at {:?} | {:?} | {:?}",
                &variant,
                &gnomad_mtdna,
                &helix_mtdb
            );
            let key: Vec<u8> = variant.into();
            let mut value = [0u8; 24];
            MtRecord {
                gnomad_mtdna,
                helix_mtdb,
            }
            .to_buf(&mut value);
            db.put_cf(cf_mtdna, key, value)?;

            chrmt_written += 1;

            Ok(())
        })?;
    }
    tracing::info!("  wrote {} chrMT records", chrmt_written);

    tracing::info!("Opening gnomAD mtDNA file(s)");
    tracing::info!("Opening HelixMtDb file(s)");

    // Finally, compact manually.
    tracing::info!("Enforcing manual compaction");
    db.compact_range_cf(cf_ndna, None::<&[u8]>, None::<&[u8]>);
    db.compact_range_cf(cf_mtdna, None::<&[u8]>, None::<&[u8]>);

    tracing::info!("Done building sequence variant frequency table");
    Ok(())
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use super::*;

    #[test]
    fn test_mt_reader() -> Result<(), anyhow::Error> {
        let mut reader = MtReader::new(
            Some("tests/data/db/create/seqvar_freqs/gnomad.chrM.vcf"),
            Some("tests/data/db/create/seqvar_freqs/helix.chrM.vcf"),
            Some(Assembly::Grch37p10),
        )?;

        {
            let mut called = false;
            let res = reader.run(|var, gnomad, helix| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("chrM"),
                        pos: 3,
                        reference: String::from("T"),
                        alternative: String::from("C")
                    }
                );
                assert_eq!(
                    gnomad,
                    MtCounts {
                        an: 56434,
                        ac_hom: 19,
                        ac_het: 1
                    }
                );
                assert_eq!(helix, MtCounts::default());
                called = true;
                Ok(())
            })?;
            assert!(res);
            assert!(called);
        }

        {
            let mut called = false;
            let res = reader.run(|var, gnomad, helix| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("chrM"),
                        pos: 5,
                        reference: String::from("A"),
                        alternative: String::from("C")
                    }
                );
                assert_eq!(
                    helix,
                    MtCounts {
                        an: 196554,
                        ac_hom: 1,
                        ac_het: 0
                    }
                );
                assert_eq!(gnomad, MtCounts::default());
                called = true;
                Ok(())
            })?;
            assert!(res);
            assert!(called);
        }

        {
            let mut called = false;
            let res = reader.run(|var, gnomad, helix| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("chrM"),
                        pos: 6,
                        reference: String::from("C"),
                        alternative: String::from("CCTCAA")
                    }
                );
                assert_eq!(
                    gnomad,
                    MtCounts {
                        an: 56433,
                        ac_hom: 0,
                        ac_het: 0
                    }
                );
                assert_eq!(helix, MtCounts::default());
                called = true;
                Ok(())
            })?;
            assert!(res);
            assert!(called);
        }

        {
            let mut called = false;
            let res = reader.run(|var, gnomad, helix| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("chrM"),
                        pos: 10,
                        reference: String::from("T"),
                        alternative: String::from("C")
                    }
                );
                assert_eq!(
                    helix,
                    MtCounts {
                        an: 196554,
                        ac_hom: 7,
                        ac_het: 1
                    }
                );
                assert_eq!(
                    gnomad,
                    MtCounts {
                        an: 56433,
                        ac_hom: 0,
                        ac_het: 0
                    }
                );
                called = true;
                Ok(())
            })?;
            assert!(res);
            assert!(called);
        }

        {
            let mut called = false;
            let res = reader.run(|var, gnomad, helix| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("chrM"),
                        pos: 11,
                        reference: String::from("C"),
                        alternative: String::from("T")
                    }
                );
                assert_eq!(
                    helix,
                    MtCounts {
                        an: 196554,
                        ac_hom: 0,
                        ac_het: 1
                    }
                );
                assert_eq!(gnomad, MtCounts::default());
                called = true;
                Ok(())
            })?;
            assert!(!res);
            assert!(called);
        }

        {
            let mut called = false;
            let res = reader.run(|_var, _gnomad, _helix| {
                called = true;
                Ok(())
            })?;
            assert!(!res);
            assert!(!called);
        }

        Ok(())
    }
}
