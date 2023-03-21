//! Annotation of sequence variants.

pub mod ann;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use clap::Parser;
use hgvs::static_data::Assembly;
use memmap2::Mmap;
use noodles::bgzf::Writer as BgzfWriter;
use noodles::vcf::header::{
    record::value::map::{info::Type, Info},
    Number,
};
use noodles::vcf::record::info::field::Value;
use noodles::vcf::{header::record::value::map::Map, Header as VcfHeader, Writer as VcfWriter};
use noodles_util::variant::reader::Builder as VariantReaderBuilder;
use rocksdb::ThreadMode;
use thousands::Separable;

use crate::common::GenomeRelease;
use crate::db::create::seqvar_freqs::reading::guess_assembly;
use crate::db::create::seqvar_freqs::serialized::vcf::Var as VcfVar;
use crate::db::create::seqvar_freqs::serialized::{
    auto::Record as AutoRecord, mt::Record as MtRecord, xy::Record as XyRecord,
};

use crate::world_flatbuffers::mehari::TxSeqDatabase;

/// Command line arguments for `annotate seqvars` sub command.
#[derive(Parser, Debug)]
#[command(about = "Annotate sequence variant VCF files", long_about = None)]
pub struct Args {
    /// Path to the mehari database folder.
    #[arg(long)]
    pub path_db: String,

    /// Genome release to use, default is to auto-detect.
    #[arg(long, value_enum)]
    pub genome_release: Option<GenomeRelease>,
    /// Path to the input VCF file.
    #[arg(long)]
    pub path_input_vcf: String,
    /// Path to the output VCF file.
    #[arg(long)]
    pub path_output_vcf: String,

    /// For debug purposes, maximal number of variants to annotate.
    #[arg(long)]
    pub max_var_count: Option<usize>,
    /// Maximal number of flatbuffers tables, should not need tweaking.
    #[arg(long, default_value_t = 5_000_000)]
    pub max_fb_tables: usize,
}

pub mod keys {
    use std::str::FromStr;

    use noodles::vcf::{
        header::info::key::Key as InfoKey, header::info::key::Other as InfoKeyOther,
    };

    lazy_static::lazy_static! {
        pub static ref GNOMAD_EXOMES_AN: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_exomes_an").unwrap());
        pub static ref GNOMAD_EXOMES_HOM: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_exomes_hom").unwrap());
        pub static ref GNOMAD_EXOMES_HET: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_exomes_het").unwrap());
        pub static ref GNOMAD_EXOMES_HEMI: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_exomes_hemi").unwrap());

        pub static ref GNOMAD_GENOMES_AN: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_genomes_an").unwrap());
        pub static ref GNOMAD_GENOMES_HOM: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_genomes_hom").unwrap());
        pub static ref GNOMAD_GENOMES_HET: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_genomes_het").unwrap());
        pub static ref GNOMAD_GENOMES_HEMI: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_genomes_hemi").unwrap());

        pub static ref GNOMAD_MTDNA_AN: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_mtdna_an").unwrap());
        pub static ref GNOMAD_MTDNA_HOM: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_mtdna_hom").unwrap());
        pub static ref GNOMAD_MTDNA_HET: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_mtdna_het").unwrap());
        pub static ref GNOMAD_MTDNA_HEMI: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_mtdna_hemi").unwrap());

        pub static ref HELIX_AN: InfoKey = InfoKey::Other(InfoKeyOther::from_str("helix_an").unwrap());
        pub static ref HELIX_HOM: InfoKey = InfoKey::Other(InfoKeyOther::from_str("helix_hom").unwrap());
        pub static ref HELIX_HET: InfoKey = InfoKey::Other(InfoKeyOther::from_str("helix_het").unwrap());
    }
}

fn build_header(header_in: &VcfHeader) -> VcfHeader {
    let mut header_out = header_in.clone();

    header_out.infos_mut().insert(
        keys::GNOMAD_EXOMES_AN.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of samples in gnomAD exomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_EXOMES_HOM.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of hom. alt. carriers in gnomAD exomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_EXOMES_HET.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of het. alt. carriers in gnomAD exomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_EXOMES_HEMI.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of hemi. alt. carriers in gnomAD exomes",
        ),
    );

    header_out.infos_mut().insert(
        keys::GNOMAD_GENOMES_AN.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of samples in gnomAD genomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_GENOMES_HOM.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of hom. alt. carriers in gnomAD genomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_GENOMES_HET.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of het. alt. carriers in gnomAD genomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_GENOMES_HEMI.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of hemi. alt. carriers in gnomAD genomes",
        ),
    );

    header_out.infos_mut().insert(
        keys::HELIX_AN.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of samples in HelixMtDb",
        ),
    );
    header_out.infos_mut().insert(
        keys::HELIX_HOM.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of hom. alt. carriers in HelixMtDb",
        ),
    );
    header_out.infos_mut().insert(
        keys::HELIX_HET.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of het. alt. carriers in HelixMtDb",
        ),
    );
    header_out
}

/// Annotate record on autosomal chromosome with gnomAD exomes/genomes.
fn annotate_record_auto<T>(
    db: &rocksdb::DBWithThreadMode<T>,
    cf: &rocksdb::ColumnFamily,
    vcf_var: VcfVar,
    vcf_record: &mut noodles::vcf::Record,
) -> Result<(), anyhow::Error>
where
    T: ThreadMode,
{
    let key: Vec<u8> = vcf_var.into();
    if let Some(freq) = db.get_cf(cf, key)? {
        let auto_record = AutoRecord::from_buf(&freq);

        vcf_record.info_mut().insert(
            keys::GNOMAD_EXOMES_AN.clone(),
            Some(Value::Integer(auto_record.gnomad_exomes.an as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_EXOMES_HOM.clone(),
            Some(Value::Integer(auto_record.gnomad_exomes.ac_hom as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_EXOMES_HET.clone(),
            Some(Value::Integer(auto_record.gnomad_exomes.ac_het as i32)),
        );

        vcf_record.info_mut().insert(
            keys::GNOMAD_GENOMES_AN.clone(),
            Some(Value::Integer(auto_record.gnomad_genomes.an as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_GENOMES_HOM.clone(),
            Some(Value::Integer(auto_record.gnomad_genomes.ac_hom as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_GENOMES_HET.clone(),
            Some(Value::Integer(auto_record.gnomad_genomes.ac_het as i32)),
        );
    };
    Ok(())
}

/// Annotate record on gonomosomal chromosome with gnomAD exomes/genomes.
fn annotate_record_xy<T>(
    db: &rocksdb::DBWithThreadMode<T>,
    cf: &rocksdb::ColumnFamily,
    vcf_var: VcfVar,
    vcf_record: &mut noodles::vcf::Record,
) -> Result<(), anyhow::Error>
where
    T: ThreadMode,
{
    let key: Vec<u8> = vcf_var.into();
    if let Some(freq) = db.get_cf(cf, key)? {
        let auto_record = XyRecord::from_buf(&freq);

        vcf_record.info_mut().insert(
            keys::GNOMAD_EXOMES_AN.clone(),
            Some(Value::Integer(auto_record.gnomad_exomes.an as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_EXOMES_HOM.clone(),
            Some(Value::Integer(auto_record.gnomad_exomes.ac_hom as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_EXOMES_HET.clone(),
            Some(Value::Integer(auto_record.gnomad_exomes.ac_het as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_EXOMES_HEMI.clone(),
            Some(Value::Integer(auto_record.gnomad_exomes.ac_hemi as i32)),
        );

        vcf_record.info_mut().insert(
            keys::GNOMAD_GENOMES_AN.clone(),
            Some(Value::Integer(auto_record.gnomad_genomes.an as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_GENOMES_HOM.clone(),
            Some(Value::Integer(auto_record.gnomad_genomes.ac_hom as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_GENOMES_HET.clone(),
            Some(Value::Integer(auto_record.gnomad_genomes.ac_het as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_GENOMES_HEMI.clone(),
            Some(Value::Integer(auto_record.gnomad_genomes.ac_hemi as i32)),
        );
    };
    Ok(())
}

/// Annotate record on mitochondrial genome with gnomAD mtDNA and HelixMtDb.
fn annotate_record_mt<T>(
    db: &rocksdb::DBWithThreadMode<T>,
    cf: &rocksdb::ColumnFamily,
    vcf_var: VcfVar,
    vcf_record: &mut noodles::vcf::Record,
) -> Result<(), anyhow::Error>
where
    T: ThreadMode,
{
    let key: Vec<u8> = vcf_var.into();
    if let Some(freq) = db.get_cf(cf, key)? {
        let mt_record = MtRecord::from_buf(&freq);

        vcf_record.info_mut().insert(
            keys::HELIX_AN.clone(),
            Some(Value::Integer(mt_record.helix_mtdb.an as i32)),
        );
        vcf_record.info_mut().insert(
            keys::HELIX_HOM.clone(),
            Some(Value::Integer(mt_record.helix_mtdb.ac_hom as i32)),
        );
        vcf_record.info_mut().insert(
            keys::HELIX_HET.clone(),
            Some(Value::Integer(mt_record.helix_mtdb.ac_het as i32)),
        );

        vcf_record.info_mut().insert(
            keys::GNOMAD_GENOMES_AN.clone(),
            Some(Value::Integer(mt_record.gnomad_mtdna.an as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_GENOMES_HOM.clone(),
            Some(Value::Integer(mt_record.gnomad_mtdna.ac_hom as i32)),
        );
        vcf_record.info_mut().insert(
            keys::GNOMAD_GENOMES_HET.clone(),
            Some(Value::Integer(mt_record.gnomad_mtdna.ac_het as i32)),
        );
    };
    Ok(())
}

lazy_static::lazy_static! {
    static ref CHROM_MT: HashSet<&'static str> = HashSet::from_iter(["M", "MT", "chrM", "chrMT"].into_iter());
    static ref CHROM_XY: HashSet<&'static str> = HashSet::from_iter(["M", "MT", "chrM", "chrMT"].into_iter());
    static ref CHROM_AUTO: HashSet<&'static str> = HashSet::from_iter([
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
        "19", "20", "21", "22", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
        "chr21", "chr22"
    ].into_iter());
}

/// Return path component for the assembly.
fn path_component(assembly: Assembly) -> &'static str {
    match assembly {
        Assembly::Grch37 | Assembly::Grch37p10 => "grch37",
        Assembly::Grch38 => "grch38",
    }
}

type IntervalTree = ArrayBackedIntervalTree<i32, u32>;

pub struct TxIntervalTrees {
    /// Mapping from contig accession to index in `trees`.
    pub contig_to_idx: HashMap<String, usize>,
    /// Interval tree to index in `TxSeqDatabase::tx_db::transcripts`, for each contig.
    pub trees: Vec<IntervalTree>,
}

impl TxIntervalTrees {
    pub fn new(db: &TxSeqDatabase) -> Self {
        let (contig_to_idx, trees) = Self::build_indices(db);
        Self {
            contig_to_idx,
            trees,
        }
    }

    fn build_indices(db: &TxSeqDatabase) -> (HashMap<String, usize>, Vec<IntervalTree>) {
        let mut contig_to_idx = HashMap::new();
        let mut trees: Vec<IntervalTree> = Vec::new();

        let mut txs = 0;

        for (tx_id, tx) in db
            .tx_db()
            .unwrap()
            .transcripts()
            .unwrap()
            .iter()
            .enumerate()
        {
            for genome_alignment in tx.genome_alignments().unwrap() {
                let contig = genome_alignment.contig().unwrap();
                let contig_idx = *contig_to_idx
                    .entry(contig.to_string())
                    .or_insert(trees.len());
                if contig_idx >= trees.len() {
                    trees.push(IntervalTree::new());
                }
                let mut start = std::i32::MAX;
                let mut stop = std::i32::MIN;
                for exon in genome_alignment.exons().unwrap().iter() {
                    start = std::cmp::min(start, exon.alt_start_i() - 1);
                    stop = std::cmp::max(stop, exon.alt_end_i());
                }
                trees[contig_idx].insert(start..stop, tx_id as u32);
            }

            txs += 1;
        }

        tracing::debug!("Loaded {} transcript", txs);
        trees.iter_mut().for_each(|t| t.index());

        (contig_to_idx, trees)
    }
}

/// Run the annotation with the given `Write` within the `VcfWriter`.
fn run_with_writer<Inner: Write>(
    mut writer: VcfWriter<Inner>,
    args: &Args,
) -> Result<(), anyhow::Error> {
    tracing::info!("Open VCF and read header");
    let mut reader = VariantReaderBuilder::default().build_from_path(&args.path_input_vcf)?;
    let header_in = reader.read_header()?;
    let header_out = build_header(&header_in);

    // Guess genome release from paths.
    let genome_release = args.genome_release.map(|gr| match gr {
        GenomeRelease::Grch37 => Assembly::Grch37p10, // has chrMT!
        GenomeRelease::Grch38 => Assembly::Grch38,
    });
    let assembly = guess_assembly(&header_in, false, genome_release)?;
    tracing::info!("Determined input assembly to be {:?}", &assembly);

    // Open the RocksDB database in read only mode.
    tracing::info!("Opening frequency database");
    let options = rocksdb::Options::default();
    let db = rocksdb::DB::open_cf_for_read_only(
        &options,
        &format!(
            "{}/seqvars/{}/freqs",
            &args.path_db,
            path_component(assembly)
        ),
        ["meta", "autosomal", "gonosomal", "mitochondrial"],
        false,
    )?;

    let cf_autosomal = db.cf_handle("autosomal").unwrap();
    let cf_gonosomal = db.cf_handle("gonosomal").unwrap();
    let cf_mtdna = db.cf_handle("mitochondrial").unwrap();

    // Open the transcript flatbuffer.
    tracing::info!("Opening transcript database");
    let tx_path = format!(
        "{}/seqvars/{}/txs.bin",
        &args.path_db,
        path_component(assembly)
    );
    let tx_file = File::open(tx_path)?;
    let tx_mmap = unsafe { Mmap::map(&tx_file)? };
    let mut fb_opts = flatbuffers::VerifierOptions::default();
    fb_opts.max_tables = args.max_fb_tables;
    let tx_db = flatbuffers::root_with_opts::<TxSeqDatabase>(&fb_opts, &tx_mmap)?;
    tracing::info!("Building transcript interval trees ...");
    let _tx_trees = TxIntervalTrees::new(&tx_db);
    tracing::info!("... done building transcript interval trees");

    // Perform the VCf annotation.
    tracing::info!("Annotating VCF ...");
    let start = Instant::now();
    let mut total_written = 0usize;

    writer.write_header(&header_out)?;
    let mut records = reader.records(&header_in);
    loop {
        if let Some(record) = records.next() {
            let mut vcf_record = record?;
            let vcf_var = VcfVar::from_vcf(&vcf_record);

            if CHROM_AUTO.contains(vcf_var.chrom.as_str()) {
                annotate_record_auto(&db, cf_autosomal, vcf_var, &mut vcf_record)?;
            } else if CHROM_XY.contains(vcf_var.chrom.as_str()) {
                annotate_record_xy(&db, cf_gonosomal, vcf_var, &mut vcf_record)?;
            } else if CHROM_MT.contains(vcf_var.chrom.as_str()) {
                annotate_record_mt(&db, cf_mtdna, vcf_var, &mut vcf_record)?;
            } else {
                tracing::trace!(
                    "Record @{:?} on non-canonical chromosomoe, skipping.",
                    &vcf_var
                );
            }

            writer.write_record(&vcf_record)?;
        } else {
            break; // all done
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
    Ok(())
}

/// Main entry point for `annotate seqvars` sub command.
pub fn run(_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    if args.path_output_vcf.ends_with(".vcf.gz") || args.path_output_vcf.ends_with(".vcf.bgzf") {
        let writer = VcfWriter::new(
            File::create(&args.path_output_vcf)
                .map(BufWriter::new)
                .map(BgzfWriter::new)?,
        );
        run_with_writer(writer, args)?;
    } else {
        let writer = VcfWriter::new(File::create(&args.path_output_vcf).map(BufWriter::new)?);
        run_with_writer(writer, args)?;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use clap_verbosity_flag::Verbosity;
    use pretty_assertions::assert_eq;
    use temp_testdir::TempDir;

    use super::{run, Args};

    #[test]
    fn smoke_test() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.vcf");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            genome_release: None,
            path_db: String::from("tests/data/annotate/db"),
            path_input_vcf: String::from(
                "tests/data/db/create/seqvar_freqs/db-rs1263393206/input.vcf",
            ),
            path_output_vcf: path_out.into_os_string().into_string().unwrap(),
            max_var_count: None,
            max_fb_tables: 5_000_000,
        };

        run(&args_common, &args)?;

        let actual = std::fs::read_to_string(args.path_output_vcf)?;
        let expected = std::fs::read_to_string(
            "tests/data/db/create/seqvar_freqs/db-rs1263393206/output.vcf",
        )?;
        assert_eq!(&expected, &actual);

        Ok(())
    }
}
