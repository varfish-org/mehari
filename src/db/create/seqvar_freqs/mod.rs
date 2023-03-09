//! Population frequencies for sequence variants.

pub mod reading;
pub mod serialized;

use std::time::Instant;

use clap::Parser;
use hgvs::static_data::Assembly;

use rocksdb::{DBWithThreadMode, SingleThreaded};

use self::serialized::auto::Record as AutoRecord;
use self::serialized::mt::Record as MtRecord;
use self::serialized::xy::Record as XyRecord;

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

    /// Path(s) to the autosomal gnomAD exomes VCF file(s).
    #[arg(long)]
    pub path_gnomad_exomes_auto: Option<Vec<String>>,
    /// Path(s) to the autosomal gnomAD genomes VCF file(s).
    #[arg(long)]
    pub path_gnomad_genomes_auto: Option<Vec<String>>,
    /// Path(s) to the gonosomal gnomAD exomes VCF file(s).
    #[arg(long)]
    pub path_gnomad_exomes_xy: Option<Vec<String>>,
    /// Path(s) to the gonosomal gnomAD genomes VCF file(s).
    #[arg(long)]
    pub path_gnomad_genomes_xy: Option<Vec<String>>,
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

/// Reading of autosomal records.
pub mod auto {
    use hgvs::static_data::Assembly;
    use noodles::vcf::Record as VcfRecord;

    use super::serialized::auto::Counts as AutoCounts;
    use super::serialized::vcf::Var as VcfVar;

    use super::reading::MultiVcfReader;

    /// Helper for reading through gnomAD exomes and genomes data;
    pub struct Reader {
        /// CSV reader for the gnomAD mitochondrial records.
        gnomad_exomes_reader: Option<MultiVcfReader>,
        /// Next variant from gnomAD.
        gnomad_exomes_next: Option<VcfRecord>,
        /// CSV reader for the HelixMtDb records.
        gnomad_genomes_reader: Option<MultiVcfReader>,
        /// Next variant from gnomAD.
        gnomad_genomes_next: Option<VcfRecord>,
    }

    impl Reader {
        /// Construct new reader with optional paths to gnomAD genomes and exomes data.
        ///
        /// Optionally, you can provide an assembly to validate the VCF contigs against.
        pub fn new(
            path_gnomad_exomes: Option<&[&str]>,
            path_gnomad_genomes: Option<&[&str]>,
            assembly: Option<Assembly>,
        ) -> Result<Self, anyhow::Error> {
            let mut gnomad_exomes_reader = path_gnomad_exomes
                .map(|path_gnomad_exomes| {
                    tracing::info!("Opening gnomAD exomes file {:?}", &path_gnomad_exomes);
                    MultiVcfReader::new(path_gnomad_exomes, assembly)
                })
                .transpose()?;

            let gnomad_exomes_next =
                if let Some(gnomad_exomes_reader) = gnomad_exomes_reader.as_mut() {
                    gnomad_exomes_reader.pop()?.0
                } else {
                    None
                };

            let mut gnomad_genomes_reader = path_gnomad_genomes
                .map(|path_gnomad_genomes| {
                    tracing::info!("Opening gnomAD genomes file {:?}", &path_gnomad_genomes);
                    MultiVcfReader::new(path_gnomad_genomes, assembly)
                })
                .transpose()?;

            let gnomad_genomes_next =
                if let Some(gnomad_genomes_reader) = gnomad_genomes_reader.as_mut() {
                    gnomad_genomes_reader.pop()?.0
                } else {
                    None
                };

            Ok(Self {
                gnomad_exomes_reader,
                gnomad_genomes_reader,
                gnomad_exomes_next,
                gnomad_genomes_next,
            })
        }

        /// Run the reading of the chrMT frequencies.
        ///
        /// Returns whether there is a next value.
        pub fn run<F>(&mut self, mut func: F) -> Result<bool, anyhow::Error>
        where
            F: FnMut(VcfVar, AutoCounts, AutoCounts) -> Result<(), anyhow::Error>,
        {
            match (&self.gnomad_exomes_next, &self.gnomad_genomes_next) {
                (None, Some(gnomad_genomes)) => {
                    func(
                        VcfVar::from_vcf(gnomad_genomes),
                        AutoCounts::default(),
                        AutoCounts::from_vcf(gnomad_genomes),
                    )?;
                    self.gnomad_genomes_next =
                        self.gnomad_genomes_reader.as_mut().unwrap().pop()?.0;
                }
                (Some(gnomad_exomes), None) => {
                    func(
                        VcfVar::from_vcf(gnomad_exomes),
                        AutoCounts::from_vcf(gnomad_exomes),
                        AutoCounts::default(),
                    )?;
                    self.gnomad_exomes_next = self.gnomad_exomes_reader.as_mut().unwrap().pop()?.0;
                }
                (Some(gnomad_exomes), Some(gnomad_genomes)) => {
                    let var_gnomad_exomes = VcfVar::from_vcf(gnomad_exomes);
                    let var_gnomad_genomes = VcfVar::from_vcf(gnomad_genomes);
                    match var_gnomad_exomes.cmp(&var_gnomad_genomes) {
                        std::cmp::Ordering::Less => {
                            func(
                                var_gnomad_exomes,
                                AutoCounts::from_vcf(gnomad_exomes),
                                AutoCounts::default(),
                            )?;
                            self.gnomad_exomes_next =
                                self.gnomad_exomes_reader.as_mut().unwrap().pop()?.0;
                        }
                        std::cmp::Ordering::Equal => {
                            func(
                                var_gnomad_exomes,
                                AutoCounts::from_vcf(gnomad_exomes),
                                AutoCounts::from_vcf(gnomad_genomes),
                            )?;
                            self.gnomad_exomes_next =
                                self.gnomad_exomes_reader.as_mut().unwrap().pop()?.0;
                            self.gnomad_genomes_next =
                                self.gnomad_genomes_reader.as_mut().unwrap().pop()?.0;
                        }
                        std::cmp::Ordering::Greater => {
                            func(
                                var_gnomad_genomes,
                                AutoCounts::default(),
                                AutoCounts::from_vcf(gnomad_genomes),
                            )?;
                            self.gnomad_genomes_next =
                                self.gnomad_genomes_reader.as_mut().unwrap().pop()?.0;
                        }
                    }
                }
                (None, None) => (),
            }

            Ok(self.gnomad_exomes_next.is_some() || self.gnomad_genomes_next.is_some())
        }
    }
}

/// Reading of gonosomal records.
pub mod xy {
    use hgvs::static_data::Assembly;
    use noodles::vcf::Record as VcfRecord;

    use super::serialized::vcf::Var as VcfVar;
    use super::serialized::xy::Counts as XyCounts;

    use super::reading::MultiVcfReader;

    /// Helper for reading through gnomAD mtDNA and HelixMtDb data;
    pub struct Reader {
        /// CSV reader for the gnomAD mitochondrial records.
        gnomad_exomes_reader: Option<MultiVcfReader>,
        /// Next variant from gnomAD.
        gnomad_exomes_next: Option<VcfRecord>,
        /// CSV reader for the HelixMtDb records.
        gnomad_genomes_reader: Option<MultiVcfReader>,
        /// Next variant from gnomAD.
        gnomad_genomes_next: Option<VcfRecord>,
    }

    impl Reader {
        /// Construct new reader with optional paths to gnomAD genomes and exomes data.
        ///
        /// Optionally, you can provide an assembly to validate the VCF contigs against.
        pub fn new(
            path_gnomad_exomes: Option<&[&str]>,
            path_gnomad_genomes: Option<&[&str]>,
            assembly: Option<Assembly>,
        ) -> Result<Self, anyhow::Error> {
            let mut gnomad_exomes_reader = path_gnomad_exomes
                .map(|path_gnomad_exomes| {
                    tracing::info!("Opening gnomAD exomes file {:?}", &path_gnomad_exomes);
                    MultiVcfReader::new(path_gnomad_exomes, assembly)
                })
                .transpose()?;

            let gnomad_exomes_next =
                if let Some(gnomad_exomes_reader) = gnomad_exomes_reader.as_mut() {
                    gnomad_exomes_reader.pop()?.0
                } else {
                    None
                };

            let mut gnomad_genomes_reader = path_gnomad_genomes
                .map(|path_gnomad_genomes| {
                    tracing::info!("Opening gnomAD genomes file {:?}", &path_gnomad_genomes);
                    MultiVcfReader::new(path_gnomad_genomes, assembly)
                })
                .transpose()?;

            let gnomad_genomes_next =
                if let Some(gnomad_genomes_reader) = gnomad_genomes_reader.as_mut() {
                    gnomad_genomes_reader.pop()?.0
                } else {
                    None
                };

            Ok(Self {
                gnomad_exomes_reader,
                gnomad_genomes_reader,
                gnomad_exomes_next,
                gnomad_genomes_next,
            })
        }

        /// Run the reading of the chrMT frequencies.
        ///
        /// Returns whether there is a next value.
        pub fn run<F>(&mut self, mut func: F) -> Result<bool, anyhow::Error>
        where
            F: FnMut(VcfVar, XyCounts, XyCounts) -> Result<(), anyhow::Error>,
        {
            match (&self.gnomad_exomes_next, &self.gnomad_genomes_next) {
                (None, Some(gnomad_genomes)) => {
                    func(
                        VcfVar::from_vcf(gnomad_genomes),
                        XyCounts::default(),
                        XyCounts::from_vcf(gnomad_genomes),
                    )?;
                    self.gnomad_genomes_next =
                        self.gnomad_genomes_reader.as_mut().unwrap().pop()?.0;
                }
                (Some(gnomad_exomes), None) => {
                    func(
                        VcfVar::from_vcf(gnomad_exomes),
                        XyCounts::from_vcf(gnomad_exomes),
                        XyCounts::default(),
                    )?;
                    self.gnomad_exomes_next = self.gnomad_exomes_reader.as_mut().unwrap().pop()?.0;
                }
                (Some(gnomad_exomes), Some(gnomad_genomes)) => {
                    let var_gnomad_exomes = VcfVar::from_vcf(gnomad_exomes);
                    let var_gnomad_genomes = VcfVar::from_vcf(gnomad_genomes);
                    match var_gnomad_exomes.cmp(&var_gnomad_genomes) {
                        std::cmp::Ordering::Less => {
                            func(
                                var_gnomad_exomes,
                                XyCounts::from_vcf(gnomad_exomes),
                                XyCounts::default(),
                            )?;
                            self.gnomad_exomes_next =
                                self.gnomad_exomes_reader.as_mut().unwrap().pop()?.0;
                        }
                        std::cmp::Ordering::Equal => {
                            func(
                                var_gnomad_exomes,
                                XyCounts::from_vcf(gnomad_exomes),
                                XyCounts::from_vcf(gnomad_genomes),
                            )?;
                            self.gnomad_exomes_next =
                                self.gnomad_exomes_reader.as_mut().unwrap().pop()?.0;
                            self.gnomad_genomes_next =
                                self.gnomad_genomes_reader.as_mut().unwrap().pop()?.0;
                        }
                        std::cmp::Ordering::Greater => {
                            func(
                                var_gnomad_genomes,
                                XyCounts::default(),
                                XyCounts::from_vcf(gnomad_genomes),
                            )?;
                            self.gnomad_genomes_next =
                                self.gnomad_genomes_reader.as_mut().unwrap().pop()?.0;
                        }
                    }
                }
                (None, None) => (),
            }

            Ok(self.gnomad_exomes_next.is_some() || self.gnomad_genomes_next.is_some())
        }
    }
}

/// Reading of chrMT records.
pub mod mt {
    use hgvs::static_data::Assembly;
    use noodles::vcf::Record as VcfRecord;

    use super::serialized::mt::Counts as MtCounts;
    use super::serialized::vcf::Var as VcfVar;

    use super::reading::MultiVcfReader;

    /// Helper for reading through gnomAD mtDNA and HelixMtDb data;
    pub struct Reader {
        /// CSV reader for the gnomAD mitochondrial records.
        gnomad_reader: Option<MultiVcfReader>,
        /// Next variant from gnomAD.
        gnomad_next: Option<VcfRecord>,
        /// CSV reader for the HelixMtDb records.
        helix_reader: Option<MultiVcfReader>,
        /// Next variant from gnomAD.
        helix_next: Option<VcfRecord>,
    }

    impl Reader {
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
}

/// Import autosomal data.
fn import_autosomal(
    args: &Args,
    genome_release: Option<Assembly>,
    db: &DBWithThreadMode<SingleThreaded>,
    cf_gonosomal: &rocksdb::ColumnFamily,
) -> Result<(), anyhow::Error> {
    tracing::info!("Processing chr1, ..., chr22 data ...");
    let start = Instant::now();
    let mut chrauto_written = 0usize;
    let path_gnomad_exomes: Option<Vec<&str>> = args
        .path_gnomad_exomes_auto
        .as_ref()
        .map(|paths| paths.iter().map(|s| s.as_str()).collect());
    let path_gnomad_genomes: Option<Vec<&str>> = args
        .path_gnomad_genomes_auto
        .as_ref()
        .map(|paths| paths.iter().map(|s| s.as_str()).collect());
    let mut auto_reader = auto::Reader::new(
        path_gnomad_exomes.as_deref(),
        path_gnomad_genomes.as_deref(),
        genome_release,
    )?;

    let mut has_next = true;
    while has_next {
        has_next = auto_reader.run(|variant, gnomad_exomes, gnomad_genomes| {
            tracing::trace!(
                "at {:?} | {:?} | {:?}",
                &variant,
                &gnomad_exomes,
                &gnomad_genomes
            );

            let key: Vec<u8> = variant.into();
            let mut value = [0u8; 32];
            AutoRecord {
                gnomad_exomes,
                gnomad_genomes,
            }
            .to_buf(&mut value);
            db.put_cf(cf_gonosomal, key, value)?;

            chrauto_written += 1;

            Ok(())
        })?;
    }

    tracing::info!(
        "  wrote {} chr1, ..., chr22 records in {:?}",
        chrauto_written,
        start.elapsed()
    );
    Ok(())
}

/// Import gonomosomal data.
fn import_gonomosomal(
    args: &Args,
    genome_release: Option<Assembly>,
    db: &DBWithThreadMode<SingleThreaded>,
    cf_gonosomal: &rocksdb::ColumnFamily,
) -> Result<(), anyhow::Error> {
    tracing::info!("Processing chrX and chrY data ...");
    let start = Instant::now();
    let mut chrxy_written = 0usize;
    let path_gnomad_exomes: Option<Vec<&str>> = args
        .path_gnomad_exomes_xy
        .as_ref()
        .map(|paths| paths.iter().map(|s| s.as_str()).collect());
    let path_gnomad_genomes: Option<Vec<&str>> = args
        .path_gnomad_genomes_xy
        .as_ref()
        .map(|paths| paths.iter().map(|s| s.as_str()).collect());
    let mut xy_reader = xy::Reader::new(
        path_gnomad_exomes.as_deref(),
        path_gnomad_genomes.as_deref(),
        genome_release,
    )?;

    let mut has_next = true;
    while has_next {
        has_next = xy_reader.run(|variant, gnomad_exomes, gnomad_genomes| {
            tracing::trace!(
                "at {:?} | {:?} | {:?}",
                &variant,
                &gnomad_exomes,
                &gnomad_genomes
            );

            let key: Vec<u8> = variant.into();
            let mut value = [0u8; 32];
            XyRecord {
                gnomad_exomes,
                gnomad_genomes,
            }
            .to_buf(&mut value);
            db.put_cf(cf_gonosomal, key, value)?;

            chrxy_written += 1;

            Ok(())
        })?;
    }

    tracing::info!(
        "  wrote {} chrX and chrY records in {:?}",
        chrxy_written,
        start.elapsed()
    );
    Ok(())
}

/// Import chrMT frequencies.
fn import_chrmt(
    args: &Args,
    genome_release: Option<Assembly>,
    db: &DBWithThreadMode<SingleThreaded>,
    cf_mtdna: &rocksdb::ColumnFamily,
) -> Result<(), anyhow::Error> {
    tracing::info!("Processing chrMT data ...");
    let start = Instant::now();

    let mut chrmt_written = 0usize;
    let mut mt_reader = mt::Reader::new(
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

    tracing::info!(
        "  wrote {} chrMT records in {:?}",
        chrmt_written,
        start.elapsed()
    );
    Ok(())
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
        ["meta", "autosomal", "gonosomal", "mitochondrial"],
    )?;

    let cf_meta = db.cf_handle("meta").unwrap();
    let cf_autosomal = db.cf_handle("autosomal").unwrap();
    let cf_gonosomal = db.cf_handle("gonosomal").unwrap();
    let cf_mtdna = db.cf_handle("mitochondrial").unwrap();

    tracing::info!("Writing meta data to database");
    db.put_cf(cf_meta, "genome-release", format!("{genome_release:?}"))?;

    // Import gnomAD variants in a chromosome-wise fashion.
    tracing::info!("Processing autosomal variant data ...");

    import_autosomal(args, genome_release, &db, cf_autosomal)?;
    import_gonomosomal(args, genome_release, &db, cf_gonosomal)?;
    import_chrmt(args, genome_release, &db, cf_mtdna)?;

    // Finally, compact manually.
    tracing::info!("Enforcing manual compaction");
    db.compact_range_cf(cf_autosomal, None::<&[u8]>, None::<&[u8]>);
    db.compact_range_cf(cf_gonosomal, None::<&[u8]>, None::<&[u8]>);
    db.compact_range_cf(cf_mtdna, None::<&[u8]>, None::<&[u8]>);

    tracing::info!("Done building sequence variant frequency table");
    Ok(())
}

#[cfg(test)]
mod test {
    use hgvs::static_data::Assembly;
    use pretty_assertions::assert_eq;

    use super::mt::Reader as MtReader;
    use super::serialized::mt::Counts as MtCounts;
    use super::serialized::vcf::Var as VcfVar;
    use super::serialized::xy::Counts as XyCounts;
    use super::xy::Reader as XyReader;

    #[test]
    fn test_mt_reader() -> Result<(), anyhow::Error> {
        let mut reader = MtReader::new(
            Some("tests/data/db/create/seqvar_freqs/mt/gnomad.chrM.vcf"),
            Some("tests/data/db/create/seqvar_freqs/mt/helix.chrM.vcf"),
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

    #[test]
    fn test_auto_reader() -> Result<(), anyhow::Error> {
        let mut reader = XyReader::new(
            Some(&[
                "tests/data/db/create/seqvar_freqs/12-37/gnomad.exomes.r2.1.1.sites.chr1.vcf",
                "tests/data/db/create/seqvar_freqs/12-37/gnomad.exomes.r2.1.1.sites.chr2.vcf",
            ]),
            Some(&[
                "tests/data/db/create/seqvar_freqs/12-37/gnomad.genomes.r2.1.1.sites.chr1.vcf",
                "tests/data/db/create/seqvar_freqs/12-37/gnomad.genomes.r2.1.1.sites.chr2.vcf",
            ]),
            Some(Assembly::Grch37p10),
        )?;

        {
            let mut called = false;
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("1"),
                        pos: 60003,
                        reference: String::from("A"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(exomes, XyCounts::default());
                assert_eq!(
                    genomes,
                    XyCounts {
                        an: 2488,
                        ac_hom: 0,
                        ac_het: 10,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("1"),
                        pos: 60008,
                        reference: String::from("T"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(exomes, XyCounts::default());
                assert_eq!(
                    genomes,
                    XyCounts {
                        an: 8576,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("1"),
                        pos: 60009,
                        reference: String::from("A"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(exomes, XyCounts::default());
                assert_eq!(
                    genomes,
                    XyCounts {
                        an: 3214,
                        ac_hom: 1,
                        ac_het: 3,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("1"),
                        pos: 200808,
                        reference: String::from("C"),
                        alternative: String::from("A")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 249706,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("1"),
                        pos: 200808,
                        reference: String::from("C"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 249706,
                        ac_hom: 0,
                        ac_het: 1,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("1"),
                        pos: 200808,
                        reference: String::from("C"),
                        alternative: String::from("T")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 249706,
                        ac_hom: 0,
                        ac_het: 1,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("2"),
                        pos: 2654979,
                        reference: String::from("G"),
                        alternative: String::from("C")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 67174,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 2
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("2"),
                        pos: 2655003,
                        reference: String::from("C"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 67694,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 1
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("2"),
                        pos: 2655009,
                        reference: String::from("A"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 67729,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 1
                    }
                );
                called = true;
                Ok(())
            })?;
            assert!(!res);
            assert!(called);
        }

        {
            let mut called = false;
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("2"),
                        pos: 2654979,
                        reference: String::from("G"),
                        alternative: String::from("C")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 67174,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 2
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("2"),
                        pos: 2655003,
                        reference: String::from("C"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 67694,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 1
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("2"),
                        pos: 2655009,
                        reference: String::from("A"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 67729,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 1
                    }
                );
                called = true;
                Ok(())
            })?;
            assert!(!res);
            assert!(called);
        }

        {
            let mut called = false;
            let res = reader.run(|_var, _exomes, _genomes| {
                called = true;
                Ok(())
            })?;
            assert!(!res);
            assert!(!called);
        }

        Ok(())
    }

    #[test]
    fn test_xy_reader() -> Result<(), anyhow::Error> {
        let mut reader = XyReader::new(
            Some(&[
                "tests/data/db/create/seqvar_freqs/xy-37/gnomad.exomes.r2.1.1.sites.chrX.vcf",
                "tests/data/db/create/seqvar_freqs/xy-37/gnomad.exomes.r2.1.1.sites.chrY.vcf",
            ]),
            Some(&["tests/data/db/create/seqvar_freqs/xy-37/gnomad.genomes.r2.1.1.sites.chrX.vcf"]),
            Some(Assembly::Grch37p10),
        )?;

        {
            let mut called = false;
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("X"),
                        pos: 60003,
                        reference: String::from("A"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(exomes, XyCounts::default());
                assert_eq!(
                    genomes,
                    XyCounts {
                        an: 2488,
                        ac_hom: 0,
                        ac_het: 10,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("X"),
                        pos: 60008,
                        reference: String::from("T"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(exomes, XyCounts::default());
                assert_eq!(
                    genomes,
                    XyCounts {
                        an: 8576,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("X"),
                        pos: 60009,
                        reference: String::from("A"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(exomes, XyCounts::default());
                assert_eq!(
                    genomes,
                    XyCounts {
                        an: 3214,
                        ac_hom: 1,
                        ac_het: 3,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("X"),
                        pos: 200808,
                        reference: String::from("C"),
                        alternative: String::from("A")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 249706,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("X"),
                        pos: 200808,
                        reference: String::from("C"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 249706,
                        ac_hom: 0,
                        ac_het: 1,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("X"),
                        pos: 200808,
                        reference: String::from("C"),
                        alternative: String::from("T")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 249706,
                        ac_hom: 0,
                        ac_het: 1,
                        ac_hemi: 0
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("Y"),
                        pos: 2654979,
                        reference: String::from("G"),
                        alternative: String::from("C")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 67174,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 2
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("Y"),
                        pos: 2655003,
                        reference: String::from("C"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 67694,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 1
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
            let res = reader.run(|var, exomes, genomes| {
                assert_eq!(
                    var,
                    VcfVar {
                        chrom: String::from("Y"),
                        pos: 2655009,
                        reference: String::from("A"),
                        alternative: String::from("G")
                    }
                );
                assert_eq!(genomes, XyCounts::default());
                assert_eq!(
                    exomes,
                    XyCounts {
                        an: 67729,
                        ac_hom: 0,
                        ac_het: 0,
                        ac_hemi: 1
                    }
                );
                called = true;
                Ok(())
            })?;
            assert!(!res);
            assert!(called);
        }

        {
            let mut called = false;
            let res = reader.run(|_var, _exomes, _genomes| {
                called = true;
                Ok(())
            })?;
            assert!(!res);
            assert!(!called);
        }

        Ok(())
    }
}
