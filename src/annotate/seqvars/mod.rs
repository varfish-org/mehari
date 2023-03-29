//! Annotation of sequence variants.

pub mod ann;
pub mod csq;
pub mod provider;

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::rc::Rc;
use std::time::Instant;

use clap::Parser;
use enumset::EnumSet;
use flatbuffers::VerifierOptions;
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

use crate::annotate::seqvars::csq::{ConsequencePredictor, VcfVariant};
use crate::annotate::seqvars::provider::MehariProvider;
use crate::common::GenomeRelease;
use crate::db::create::seqvar_clinvar::serialize::Record as ClinvarRecord;
use crate::db::create::seqvar_freqs::reading::guess_assembly;
use crate::db::create::seqvar_freqs::serialized::vcf::Var as VcfVar;
use crate::db::create::seqvar_freqs::serialized::{
    auto::Record as AutoRecord, mt::Record as MtRecord, xy::Record as XyRecord,
};

use crate::db::create::txs::data::{
    ExonAlignment, GeneToTxId, GenomeAlignment, GenomeBuild, SequenceDb, Strand, Transcript,
    TranscriptBiotype, TranscriptDb, TranscriptTag, TxSeqDatabase,
};
use crate::world_flatbuffers::mehari::{
    GenomeBuild as FlatGenomeBuild, Strand as FlatStrand,
    TranscriptBiotype as FlatTranscriptBiotype, TxSeqDatabase as FlatTxSeqDatabase,
};

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
    /// Maximal number of flatbuffers tables, should not need tweaking.  However, if you see a
    /// "too many tables" error output, increase this value.
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

        pub static ref HELIX_AN: InfoKey = InfoKey::Other(InfoKeyOther::from_str("helix_an").unwrap());
        pub static ref HELIX_HOM: InfoKey = InfoKey::Other(InfoKeyOther::from_str("helix_hom").unwrap());
        pub static ref HELIX_HET: InfoKey = InfoKey::Other(InfoKeyOther::from_str("helix_het").unwrap());

        pub static ref ANN: InfoKey = InfoKey::Other(InfoKeyOther::from_str("ANN").unwrap());

        pub static ref CLINVAR_PATHO: InfoKey = InfoKey::Other(InfoKeyOther::from_str("clinvar_patho").unwrap());
        pub static ref CLINVAR_VCV: InfoKey = InfoKey::Other(InfoKeyOther::from_str("clinvar_vcv").unwrap());
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

    header_out.infos_mut().insert(
        keys::ANN.clone(),
        Map::<Info>::new(
            Number::Unknown,
            Type::String,
            "Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | \
            Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | \
            cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | \
            ERRORS / WARNINGS / INFO'",
        ),
    );

    header_out.infos_mut().insert(
        keys::CLINVAR_PATHO.clone(),
        Map::<Info>::new(Number::Count(1), Type::String, "ClinVar pathogenicity"),
    );
    header_out.infos_mut().insert(
        keys::CLINVAR_VCV.clone(),
        Map::<Info>::new(Number::Count(1), Type::String, "ClinVar VCV accession"),
    );

    header_out
}

/// Annotate record on autosomal chromosome with gnomAD exomes/genomes.
fn annotate_record_auto<T>(
    db: &rocksdb::DBWithThreadMode<T>,
    cf: &rocksdb::ColumnFamily,
    key: &Vec<u8>,
    vcf_record: &mut noodles::vcf::Record,
) -> Result<(), anyhow::Error>
where
    T: ThreadMode,
{
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
    key: &Vec<u8>,
    vcf_record: &mut noodles::vcf::Record,
) -> Result<(), anyhow::Error>
where
    T: ThreadMode,
{
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
    key: &Vec<u8>,
    vcf_record: &mut noodles::vcf::Record,
) -> Result<(), anyhow::Error>
where
    T: ThreadMode,
{
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

/// Annotate record with ClinVar information.
fn annotate_record_clinvar<T>(
    db: &rocksdb::DBWithThreadMode<T>,
    cf: &rocksdb::ColumnFamily,
    key: &Vec<u8>,
    vcf_record: &mut noodles::vcf::Record,
) -> Result<(), anyhow::Error>
where
    T: ThreadMode,
{
    if let Some(clinvar_anno) = db.get_cf(cf, key)? {
        let clinvar_record: ClinvarRecord = bincode::deserialize(&clinvar_anno)?;

        let ClinvarRecord {
            summary_clinvar_pathogenicity,
            vcv,
            ..
        } = clinvar_record;

        vcf_record.info_mut().insert(
            keys::CLINVAR_PATHO.clone(),
            Some(Value::String(
                summary_clinvar_pathogenicity.first().unwrap().to_string(),
            )),
        );
        vcf_record
            .info_mut()
            .insert(keys::CLINVAR_VCV.clone(), Some(Value::String(vcv.clone())));
    }

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

/// Load flatbuffers transcripts.
pub fn load_tx_db(tx_path: &str, max_fb_tables: usize) -> Result<TxSeqDatabase, anyhow::Error> {
    let tx_file = File::open(tx_path)?;
    let tx_mmap = unsafe { Mmap::map(&tx_file)? };
    let fb_opts = VerifierOptions {
        max_tables: max_fb_tables,
        ..Default::default()
    };
    let fb_tx_db = flatbuffers::root_with_opts::<FlatTxSeqDatabase>(&fb_opts, &tx_mmap)?;

    let transcripts = fb_tx_db
        .tx_db()
        .unwrap()
        .transcripts()
        .unwrap()
        .into_iter()
        .map(|rec| {
            let mut tags = EnumSet::new();
            let flat_tags = rec.biotype().0;
            if flat_tags & 1 != 0 {
                tags.insert(TranscriptTag::Basic);
            }
            if flat_tags & 2 != 0 {
                tags.insert(TranscriptTag::EnsemblCanonical);
            }
            if flat_tags & 4 != 0 {
                tags.insert(TranscriptTag::ManeSelect);
            }
            if flat_tags & 8 != 0 {
                tags.insert(TranscriptTag::ManePlusClinical);
            }
            if flat_tags & 16 != 0 {
                tags.insert(TranscriptTag::RefSeqSelect);
            }

            let genome_alignments = rec
                .genome_alignments()
                .unwrap()
                .iter()
                .map(|rec| {
                    let exons = rec
                        .exons()
                        .unwrap()
                        .iter()
                        .map(|rec| ExonAlignment {
                            alt_start_i: rec.alt_start_i(),
                            alt_end_i: rec.alt_end_i(),
                            ord: rec.ord(),
                            alt_cds_start_i: if rec.alt_cds_start_i() >= 0 {
                                Some(rec.alt_cds_start_i())
                            } else {
                                None
                            },
                            alt_cds_end_i: if rec.alt_cds_end_i() >= 0 {
                                Some(rec.alt_cds_end_i())
                            } else {
                                None
                            },
                            cigar: rec.cigar().unwrap().to_string(),
                        })
                        .collect();

                    GenomeAlignment {
                        genome_build: match rec.genome_build() {
                            FlatGenomeBuild::Grch37 => GenomeBuild::Grch37,
                            FlatGenomeBuild::Grch38 => GenomeBuild::Grch38,
                            _ => panic!("Invalid genome build: {:?}", rec.genome_build()),
                        },
                        contig: rec.contig().unwrap().to_string(),
                        cds_start: if rec.cds_start() >= 0 {
                            Some(rec.cds_start())
                        } else {
                            None
                        },
                        cds_end: if rec.cds_end() >= 0 {
                            Some(rec.cds_end())
                        } else {
                            None
                        },
                        strand: match rec.strand() {
                            FlatStrand::Plus => Strand::Plus,
                            FlatStrand::Minus => Strand::Minus,
                            _ => panic!("invalid strand: {:?}", rec.strand()),
                        },
                        exons,
                    }
                })
                .collect();

            Transcript {
                id: rec.id().unwrap().to_string(),
                gene_name: rec.gene_name().unwrap().to_string(),
                gene_id: rec.gene_id().unwrap().to_string(),
                biotype: match rec.biotype() {
                    FlatTranscriptBiotype::Coding => TranscriptBiotype::Coding,
                    FlatTranscriptBiotype::NonCoding => TranscriptBiotype::NonCoding,
                    _ => panic!("Invalid biotype: {:?}", rec.biotype()),
                },
                tags,
                protein: if rec.protein().is_none() || rec.protein().unwrap().is_empty() {
                    None
                } else {
                    Some(rec.protein().unwrap().to_string())
                },
                start_codon: if rec.start_codon() >= 0 {
                    Some(rec.start_codon())
                } else {
                    None
                },
                stop_codon: if rec.stop_codon() >= 0 {
                    Some(rec.stop_codon())
                } else {
                    None
                },
                genome_alignments,
            }
        })
        .collect();

    let gene_to_tx = fb_tx_db
        .tx_db()
        .unwrap()
        .gene_to_tx()
        .unwrap()
        .into_iter()
        .map(|rec| GeneToTxId {
            gene_name: rec.gene_name().unwrap().to_string(),
            tx_ids: rec
                .tx_ids()
                .unwrap()
                .into_iter()
                .map(|s| s.to_string())
                .collect(),
        })
        .collect();

    let tx_db = TranscriptDb {
        transcripts,
        gene_to_tx,
    };

    let seq_db = SequenceDb {
        aliases: fb_tx_db
            .seq_db()
            .unwrap()
            .aliases()
            .unwrap()
            .into_iter()
            .map(|alias| alias.to_string())
            .collect(),
        aliases_idx: fb_tx_db
            .seq_db()
            .unwrap()
            .aliases_idx()
            .unwrap()
            .into_iter()
            .collect(),
        seqs: fb_tx_db
            .seq_db()
            .unwrap()
            .seqs()
            .unwrap()
            .into_iter()
            .map(|seq| seq.to_string())
            .collect(),
    };

    Ok(TxSeqDatabase { tx_db, seq_db })
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

    // Open the frequency RocksDB database in read only mode.
    tracing::info!("Opening frequency database");
    let rocksdb_path = format!(
        "{}/seqvars/{}/freqs",
        &args.path_db,
        path_component(assembly)
    );
    tracing::debug!("RocksDB path = {}", &rocksdb_path);
    let options = rocksdb::Options::default();
    let db_freq = rocksdb::DB::open_cf_for_read_only(
        &options,
        &rocksdb_path,
        ["meta", "autosomal", "gonosomal", "mitochondrial"],
        false,
    )?;

    let cf_autosomal = db_freq.cf_handle("autosomal").unwrap();
    let cf_gonosomal = db_freq.cf_handle("gonosomal").unwrap();
    let cf_mtdna = db_freq.cf_handle("mitochondrial").unwrap();

    // Open the ClinVar RocksDB database in read only mode.
    tracing::info!("Opening ClinVar database");
    let rocksdb_path = format!(
        "{}/seqvars/{}/clinvar",
        &args.path_db,
        path_component(assembly)
    );
    tracing::debug!("RocksDB path = {}", &rocksdb_path);
    let options = rocksdb::Options::default();
    let db_clinvar = rocksdb::DB::open_cf_for_read_only(
        &options,
        &rocksdb_path,
        ["meta", "clinvar_seqvars"],
        false,
    )?;

    let cf_clinvar = db_clinvar.cf_handle("clinvar_seqvars").unwrap();

    // Open the transcript flatbuffer.
    tracing::info!("Opening transcript database");
    let tx_db = load_tx_db(
        &format!(
            "{}/seqvars/{}/txs.bin",
            &args.path_db,
            path_component(assembly)
        ),
        args.max_fb_tables,
    )?;
    tracing::info!("Building transcript interval trees ...");
    let provider = Rc::new(MehariProvider::new(tx_db, assembly));
    let predictor = ConsequencePredictor::new(provider, assembly);
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

            // TODO: ignores all but the first alternative allele!
            let vcf_var = VcfVar::from_vcf(&vcf_record);

            if vcf_var.pos == 874778 {
                println!("HAHAHAHA");
            }

            // Build key for RocksDB database from `vcf_var`.
            let key: Vec<u8> = vcf_var.clone().into();

            // Annotate with frequency.
            if CHROM_AUTO.contains(vcf_var.chrom.as_str()) {
                annotate_record_auto(&db_freq, cf_autosomal, &key, &mut vcf_record)?;
            } else if CHROM_XY.contains(vcf_var.chrom.as_str()) {
                annotate_record_xy(&db_freq, cf_gonosomal, &key, &mut vcf_record)?;
            } else if CHROM_MT.contains(vcf_var.chrom.as_str()) {
                annotate_record_mt(&db_freq, cf_mtdna, &key, &mut vcf_record)?;
            } else {
                tracing::trace!(
                    "Record @{:?} on non-canonical chromosome, skipping.",
                    &vcf_var
                );
            }

            // Annotate with ClinVar.
            annotate_record_clinvar(&db_clinvar, cf_clinvar, &key, &mut vcf_record)?;

            let VcfVar {
                chrom,
                pos,
                reference,
                alternative,
            } = vcf_var;

            // Annotate with variant effect.
            if let Some(ann_fields) = predictor.predict(&VcfVariant {
                chromosome: chrom,
                position: pos as i32,
                reference,
                alternative,
            })? {
                if !ann_fields.is_empty() {
                    vcf_record.info_mut().insert(
                        keys::ANN.clone(),
                        Some(Value::StringArray(
                            ann_fields.iter().map(|ann| Some(ann.to_string())).collect(),
                        )),
                    );
                }
            }

            // Annotate with ClinVar annotation.

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
    tracing::info!("config = {:#?}", &args);
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
