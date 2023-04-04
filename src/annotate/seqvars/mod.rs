//! Annotation of sequence variants.

pub mod ann;
pub mod csq;
pub mod provider;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::ops::Deref;
use std::path::Path;
use std::rc::Rc;
use std::str::FromStr;
use std::time::Instant;

use clap::{Args as ClapArgs, Parser};
use enumset::EnumSet;
use flatbuffers::VerifierOptions;
use flate2::write::GzEncoder;
use flate2::Compression;
use hgvs::static_data::Assembly;
use memmap2::Mmap;
use noodles::bgzf::Writer as BgzfWriter;
use noodles::vcf::header::format::key::{
    CONDITIONAL_GENOTYPE_QUALITY, GENOTYPE, READ_DEPTH, READ_DEPTHS,
};
use noodles::vcf::header::{
    record::value::map::{info::Type, Info},
    Number,
};
use noodles::vcf::record::info::field::Value;
use noodles::vcf::record::Chromosome;
use noodles::vcf::{
    header::format::key::Key as FormatKey, header::format::key::Other as FormatKeyOther,
};
use noodles::vcf::{
    header::record::value::map::Map, Header as VcfHeader, Record as VcfRecord, Writer as VcfWriter,
};
use noodles_util::variant::reader::Builder as VariantReaderBuilder;
use rocksdb::ThreadMode;
use thousands::Separable;

use crate::annotate::seqvars::csq::{ConsequencePredictor, VcfVariant};
use crate::annotate::seqvars::provider::MehariProvider;
use crate::common::GenomeRelease;
use crate::db::create::seqvar_clinvar::serialize::Record as ClinvarRecord;
use crate::db::create::seqvar_freqs::reading::{guess_assembly, is_canonical};
use crate::db::create::seqvar_freqs::serialized::vcf::Var as VcfVar;
use crate::db::create::seqvar_freqs::serialized::{
    auto::Record as AutoRecord, mt::Record as MtRecord, xy::Record as XyRecord,
};

use crate::db::create::txs::data::{
    ExonAlignment, GeneToTxId, GenomeAlignment, GenomeBuild, SequenceDb, Strand, Transcript,
    TranscriptBiotype, TranscriptDb, TranscriptTag, TxSeqDatabase,
};
use crate::ped::{PedigreeByName, Sex};
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
    /// Path to the input PED file.
    #[arg(long)]
    pub path_input_ped: String,
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
    #[arg(long)]
    pub path_output_tsv: Option<String>,
}

pub mod keys {
    use std::str::FromStr;

    use noodles::vcf::{
        header::format::key::Key as FormatKey, header::format::key::Other as FormatKeyOther,
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
            .insert(keys::CLINVAR_VCV.clone(), Some(Value::String(vcv)));
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

    static ref CHROM_TO_CHROM_NO: HashMap<String, u32> = {
        let mut m = HashMap::new();

        for i in 1..=22 {
            m.insert(format!("chr{}", i), i);
            m.insert(format!("{}", i), i);
        }
        m.insert(String::from("X"), 23);
        m.insert(String::from("chrX"), 23);
        m.insert(String::from("Y"), 23);
        m.insert(String::from("chrY"), 23);
        m.insert(String::from("M"), 23);
        m.insert(String::from("chrM"), 23);
        m.insert(String::from("MT"), 24);
        m.insert(String::from("chrMT"), 24);

        m
    };
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

/// Mehari-local trait for writing out annotated VCF records as VCF or VarFish TSV.
pub trait AnnotatedVcfWriter {
    fn write_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error>;
    fn write_record(&mut self, record: &VcfRecord) -> Result<(), anyhow::Error>;
    fn set_assembly(&mut self, assembly: &Assembly) { /* nop */
    }
    fn set_pedigree(&mut self, pedigree: &PedigreeByName) { /* nop */
    }
}

/// Implement `AnnotatedVcfWriter` for `VcfWriter`.
impl<Inner: Write> AnnotatedVcfWriter for VcfWriter<Inner> {
    fn write_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error> {
        self.write_header(header)
            .map_err(|e| anyhow::anyhow!("Error writing VCF header: {}", e))
    }

    fn write_record(&mut self, record: &VcfRecord) -> Result<(), anyhow::Error> {
        self.write_record(record)
            .map_err(|e| anyhow::anyhow!("Error writing VCF record: {}", e))
    }
}

/// Writing of VarFish TSV files.
struct VarFishTsvWriter {
    inner: Box<dyn Write>,
    assembly: Option<Assembly>,
    pedigree: Option<PedigreeByName>,
    header: Option<VcfHeader>,
}

// Offsets for UCSC binning.
//
// cf. http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
static BIN_OFFSETS: &[i32] = &[512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0];
const BIN_FIRST_SHIFT: i32 = 17;
const BIN_NEXT_SHIFT: i32 = 3;

// Compute UCSC bin from 0-based half-open interval.
pub fn bin_from_range(begin: i32, end: i32) -> Result<i32, anyhow::Error> {
    let mut begin_bin = begin >> BIN_FIRST_SHIFT;
    let mut end_bin = std::cmp::max(begin, end - 1) >> BIN_FIRST_SHIFT;

    for offset in BIN_OFFSETS {
        if begin_bin == end_bin {
            return Ok(offset + begin_bin);
        }
        begin_bin >>= BIN_NEXT_SHIFT;
        end_bin >>= BIN_NEXT_SHIFT;
    }

    anyhow::bail!(
        "begin {}, end {} out of range in bin_from_range (max is 512M",
        begin,
        end
    );
}

/// Entry with genotype (`gt`), coverage (`dp`), allele depth (`ad`) and
/// genotype quality (`gq`).
#[derive(Debug, Default)]
struct GenotypeInfo {
    pub name: String,
    pub gt: Option<String>,
    pub dp: Option<i32>,
    pub ad: Option<i32>,
    pub gq: Option<i32>,
}

#[derive(Debug, Default)]
struct GenotypeCalls {
    pub entries: Vec<GenotypeInfo>,
}

impl GenotypeCalls {
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
                result.push_str(&format!("\"\"\"gt\"\"\"\":\"\"\"\"{}\"\"\"\"", gt));
            }

            if prev {
                result.push(',');
            }
            if let Some(ad) = &entry.ad {
                prev = true;
                result.push_str(&format!("\"\"\"ad\"\"\"\":{}", ad));
            }

            if prev {
                result.push(',');
            }
            if let Some(dp) = &entry.dp {
                prev = true;
                result.push_str(&format!("\"\"\"dp\"\"\"\":{}", dp));
            }

            if prev {
                result.push(',');
            }
            if let Some(gq) = &entry.gq {
                prev = true;
                result.push_str(&format!("\"\"\"gq\"\"\"\":{}", gq));
            }

            result.push('}');
        }

        result.push('}');
        result
    }
}

impl VarFishTsvWriter {
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
            assembly: None,
            pedigree: None,
            header: None,
        }
    }

    /// Fill `record` coordinate fields.
    ///
    /// # Returns
    ///
    /// Returns `true` if the record is to be written to the TSV file and `false` if
    /// it is to be skipped because it was not on a canonical chromosome (chr1..chr22,
    /// chrX, chrY, chrM/chrMT).
    fn fill_coords(
        &self,
        assembly: Assembly,
        record: &VcfRecord,
        tsv_record: &mut VarFishTsvRecord,
    ) -> Result<bool, anyhow::Error> {
        tsv_record.release = match assembly {
            Assembly::Grch37 | Assembly::Grch37p10 => String::from("GRCh37"),
            Assembly::Grch38 => String::from("GRCh38"),
        };
        if let Chromosome::Name(name) = record.chromosome() {
            // strip "chr" prefix from `name`
            let name = if name.starts_with("chr") {
                &name[3..]
            } else {
                name
            };
            // add back "chr" prefix if necessary (for GRCh38)
            tsv_record.chromosome = if assembly == Assembly::Grch38 {
                format!("chr{}", name)
            } else {
                name.to_string()
            };
        } else {
            anyhow::bail!("Cannot handle chromosome: {:?}", record.chromosome())
        }
        if let Some(chromosome_no) = CHROM_TO_CHROM_NO.get(&tsv_record.chromosome) {
            tsv_record.chromosome_no = *chromosome_no;
        } else {
            return Ok(false);
        }

        tsv_record.reference = record.reference_bases().to_string();
        tsv_record.alternative = record.alternate_bases().deref()[0].to_string();

        tsv_record.start = record.position().into();
        tsv_record.end = tsv_record.start + tsv_record.reference.len() - 1;
        tsv_record.bin = bin_from_range(tsv_record.start as i32 - 1, tsv_record.end as i32)? as u32;

        tsv_record.var_type = if tsv_record.reference.len() == tsv_record.alternative.len() {
            if tsv_record.reference.len() == 1 {
                String::from("snv")
            } else {
                String::from("mnv")
            }
        } else {
            String::from("indel")
        };

        Ok(true)
    }

    /// Fill `record` genotype and family-local allele frequencies.
    fn fill_genotype_and_freqs(
        &self,
        record: &VcfRecord,
        tsv_record: &mut VarFishTsvRecord,
    ) -> Result<(), anyhow::Error> {
        // Extract genotype information.
        let hdr = self
            .header
            .as_ref()
            .expect("VCF header must be set/written");
        let mut gt_calls = GenotypeCalls::default();
        for (sample, genotype) in hdr.sample_names().iter().zip(record.genotypes().iter()) {
            let mut gt_info = GenotypeInfo {
                name: sample.to_string(),
                ..Default::default()
            };

            if let Some(gt) = genotype
                .get(&GENOTYPE)
                .map(|value| match value {
                    Some(noodles::vcf::record::genotypes::genotype::field::Value::String(s)) => {
                        Ok(s.to_owned())
                    }
                    _ => anyhow::bail!("invalid GT value"),
                })
                .transpose()?
            {
                let individual = self
                    .pedigree
                    .as_ref()
                    .expect("pedigree must be set")
                    .individuals
                    .get(sample)
                    .expect(&format!("individual {} not found in pedigree", sample));
                // Update per-family counts.
                if ["X", "chrX"]
                    .iter()
                    .any(|c| *c == tsv_record.chromosome.as_str())
                {
                    match individual.sex {
                        Sex::Male => {
                            if gt.contains("1") {
                                tsv_record.num_hemi_alt += 1;
                            } else {
                                tsv_record.num_hemi_ref += 1;
                            }
                        }
                        Sex::Female | Sex::Unknown => {
                            // assume diploid/female if unknown
                            let matches = gt.matches("1").count();
                            if matches == 0 {
                                tsv_record.num_hom_ref += 1;
                            } else if matches == 1 {
                                tsv_record.num_het += 1;
                            } else {
                                tsv_record.num_hom_alt += 1;
                            }
                        }
                    }
                } else if ["Y", "chrY"]
                    .iter()
                    .any(|c| *c == tsv_record.chromosome.as_str())
                {
                    if individual.sex == Sex::Male {
                        if gt.contains("1") {
                            tsv_record.num_hemi_alt += 1;
                        } else {
                            tsv_record.num_hemi_ref += 1;
                        }
                    }
                } else {
                    let matches = gt.matches("1").count();
                    if matches == 0 {
                        tsv_record.num_hom_ref += 1;
                    } else if matches == 1 {
                        tsv_record.num_het += 1;
                    } else {
                        tsv_record.num_hom_alt += 1;
                    }
                }

                // Store genotype value.
                gt_info.gt = Some(gt);
            }

            if let Some(dp) = genotype
                .get(&READ_DEPTH)
                .map(|value| match value {
                    Some(noodles::vcf::record::genotypes::genotype::field::Value::Integer(i)) => {
                        Ok(*i)
                    }
                    _ => anyhow::bail!("invalid DP value"),
                })
                .transpose()?
            {
                gt_info.dp = Some(dp);
            }

            if let Some(ad) = genotype
                .get(&READ_DEPTHS)
                .map(|value| match value {
                    Some(
                        noodles::vcf::record::genotypes::genotype::field::Value::IntegerArray(arr),
                    ) => Ok(arr[0].expect("missing AD value")),
                    _ => anyhow::bail!("invalid AD value"),
                })
                .transpose()?
            {
                gt_info.ad = Some(ad);
            }

            if let Some(gq) = genotype
                .get(&CONDITIONAL_GENOTYPE_QUALITY)
                .map(|value| match value {
                    Some(noodles::vcf::record::genotypes::genotype::field::Value::Integer(i)) => {
                        Ok(*i)
                    }
                    _ => anyhow::bail!("invalid GQ value"),
                })
                .transpose()?
            {
                gt_info.gq = Some(gq);
            }

            let x = genotype.get(&FormatKey::Other(FormatKeyOther::from_str("SQ").unwrap()));

            if let Some(sq) = genotype
                .get(&FormatKey::Other(FormatKeyOther::from_str("SQ").unwrap()))
                .map(|value| match value {
                    Some(noodles::vcf::record::genotypes::genotype::field::Value::Float(f)) => {
                        Ok(*f)
                    }
                    _ => anyhow::bail!("invalid SQ value"),
                })
                .transpose()?
            {
                gt_info.gq = Some(sq as i32);
            }

            gt_calls.entries.push(gt_info);
        }
        tsv_record.genotype = gt_calls.for_tsv();

        // Extract family counts.

        // Extract gnomAD frequencies.

        let gnomad_exomes_an = record
            .info()
            .get(&keys::GNOMAD_EXOMES_AN)
            .unwrap()
            .map(|v| match v {
                Value::Integer(value) => *value,
                _ => panic!("Unexpected value type for GNOMAD_EXOMES_AN"),
            })
            .unwrap_or_default();
        if gnomad_exomes_an > 0 {
            tsv_record.gnomad_exomes_homozygous = record
                .info()
                .get(&keys::GNOMAD_EXOMES_HOM)
                .unwrap_or(Some(&Value::Integer(0)))
                .map(|v| match v {
                    Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_EXOMES_HOM"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_exomes_heterozygous = record
                .info()
                .get(&keys::GNOMAD_EXOMES_HET)
                .unwrap_or(Some(&Value::Integer(0)))
                .map(|v| match v {
                    Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_EXOMES_HET"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_exomes_hemizygous = record
                .info()
                .get(&keys::GNOMAD_EXOMES_HEMI)
                .unwrap_or(Some(&Value::Integer(0)))
                .map(|v| match v {
                    Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_EXOMES_HEMI"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_exomes_frequency = (tsv_record.gnomad_exomes_hemizygous
                + tsv_record.gnomad_exomes_heterozygous
                + tsv_record.gnomad_exomes_homozygous * 2)
                as f64
                / gnomad_exomes_an as f64;
        }

        let gnomad_genomes_an = record
            .info()
            .get(&keys::GNOMAD_GENOMES_AN)
            .unwrap()
            .map(|v| match v {
                Value::Integer(value) => *value,
                _ => panic!("Unexpected value type for GNOMAD_GENOMES_AN"),
            })
            .unwrap_or_default();
        if gnomad_genomes_an > 0 {
            tsv_record.gnomad_genomes_homozygous = record
                .info()
                .get(&keys::GNOMAD_GENOMES_HOM)
                .unwrap_or(Some(&Value::Integer(0)))
                .map(|v| match v {
                    Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_GENOMES_HOM"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_genomes_heterozygous = record
                .info()
                .get(&keys::GNOMAD_GENOMES_HET)
                .unwrap_or(Some(&Value::Integer(0)))
                .map(|v| match v {
                    Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_GENOMES_HET"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_genomes_hemizygous = record
                .info()
                .get(&keys::GNOMAD_GENOMES_HEMI)
                .unwrap_or(Some(&Value::Integer(0)))
                .map(|v| match v {
                    Value::Integer(value) => *value,
                    _ => panic!("Unexpected value type for GNOMAD_GENOMES_HEMI"),
                })
                .unwrap_or_default();
            tsv_record.gnomad_genomes_frequency = (tsv_record.gnomad_genomes_hemizygous
                + tsv_record.gnomad_genomes_heterozygous
                + tsv_record.gnomad_genomes_homozygous * 2)
                as f64
                / gnomad_genomes_an as f64;
        }

        Ok(())
    }

    /// Fill `record` background frequencies.
    fn fill_bg_freqs(
        &self,
        record: &VcfRecord,
        tsv_record: &mut VarFishTsvRecord,
    ) -> Result<(), anyhow::Error> {
        Ok(())
    }

    /// Fill `record` RefSeq fields.
    fn fill_refseq(
        &self,
        record: &VcfRecord,
        tsv_record: &mut VarFishTsvRecord,
    ) -> Result<(), anyhow::Error> {
        Ok(())
    }

    /// Fill `record` ENSEMBL fields.
    fn fill_ensembl(
        &self,
        record: &VcfRecord,
        tsv_record: &mut VarFishTsvRecord,
    ) -> Result<(), anyhow::Error> {
        Ok(())
    }
}

/// A record, as written out to a VarFish TSV file.
#[derive(Debug, Default)]
pub struct VarFishTsvRecord {
    pub release: String,
    pub chromosome: String,
    pub chromosome_no: u32,
    pub start: usize,
    pub end: usize,
    pub bin: u32,
    pub reference: String,
    pub alternative: String,
    pub var_type: String,

    // Writing out case_id and set_info is not used anyway.
    // pub case_id: String,
    // pub set_id: String,

    // The info field is not populated anyway.
    // pub info: String,

    // TODO: lookup Java code for supported fields in genotype.
    pub genotype: String,

    pub num_hom_alt: u32,
    pub num_hom_ref: u32,
    pub num_het: u32,
    pub num_hemi_alt: u32,
    pub num_hemi_ref: u32,

    pub in_clinvar: bool,

    // NB: ExAc and 1000 Genomes are not written out anymore.
    // pub exac_frequency: String,
    // pub exac_homozygous: String,
    // pub exac_heterozygous: String,
    // pub exac_hemizygous: String,
    // pub thousand_genomes_frequency: String,
    // pub thousand_genomes_homozygous: String,
    // pub thousand_genomes_heterozygous: String,
    // pub thousand_genomes_hemizygous: String,
    pub gnomad_exomes_frequency: f64,
    pub gnomad_exomes_homozygous: i32,
    pub gnomad_exomes_heterozygous: i32,
    pub gnomad_exomes_hemizygous: i32,
    pub gnomad_genomes_frequency: f64,
    pub gnomad_genomes_homozygous: i32,
    pub gnomad_genomes_heterozygous: i32,
    pub gnomad_genomes_hemizygous: i32,

    pub refseq_gene_id: String,
    pub refseq_transcript_id: String,
    pub refseq_transcript_coding: bool,
    pub refseq_hgvs_c: String,
    pub refseq_hgvs_p: String,
    pub refseq_effect: Vec<String>,
    pub refseq_exon_dist: i32,

    pub ensembl_gene_id: String,
    pub ensembl_transcript_id: String,
    pub ensembl_transcript_coding: bool,
    pub ensembl_hgvs_c: String,
    pub ensembl_hgvs_p: String,
    pub ensembl_effect: Vec<String>,
    pub ensembl_exon_dist: i32,
}

impl VarFishTsvRecord {
    pub fn to_tsv(&self) -> Vec<String> {
        vec![
            self.release.clone(),
            self.chromosome.clone(),
            format!("{}", self.chromosome_no),
            format!("{}", self.start),
            format!("{}", self.end),
            format!("{}", self.bin),
            self.reference.clone(),
            self.alternative.clone(),
            self.var_type.clone(),
            String::from("."),
            String::from("."),
            String::from("{}"),
            self.genotype.clone(),
            format!("{}", self.num_hom_alt),
            format!("{}", self.num_hom_ref),
            format!("{}", self.num_het),
            format!("{}", self.num_hemi_alt),
            format!("{}", self.num_hemi_ref),
            if self.in_clinvar { "TRUE" } else { "FALSE" }.to_string(),
            // exac
            String::from("0"),
            String::from("0"),
            String::from("0"),
            String::from("0"),
            // thousand genomes
            String::from("0"),
            String::from("0"),
            String::from("0"),
            String::from("0"),
            format!("{}", self.gnomad_exomes_frequency),
            format!("{}", self.gnomad_exomes_homozygous),
            format!("{}", self.gnomad_exomes_heterozygous),
            format!("{}", self.gnomad_exomes_hemizygous),
            format!("{}", self.gnomad_genomes_frequency),
            format!("{}", self.gnomad_genomes_homozygous),
            format!("{}", self.gnomad_genomes_heterozygous),
            format!("{}", self.gnomad_genomes_hemizygous),
            self.refseq_gene_id.clone(),
            self.refseq_transcript_id.clone(),
            if self.refseq_transcript_coding {
                "TRUE"
            } else {
                "FALSE"
            }
            .to_string(),
            self.refseq_hgvs_c.clone(),
            self.refseq_hgvs_p.clone(),
            format!("{{{}}}", self.refseq_effect.join(",")),
            format!("{}", self.refseq_exon_dist),
            self.ensembl_gene_id.clone(),
            self.ensembl_transcript_id.clone(),
            if self.ensembl_transcript_coding {
                "TRUE"
            } else {
                "FALSE"
            }
            .to_string(),
            self.ensembl_hgvs_c.clone(),
            self.ensembl_hgvs_p.clone(),
            format!("{{{}}}", self.ensembl_effect.join(",")),
            format!("{}", self.ensembl_exon_dist),
        ]
    }
}

/// Implement `AnnotatedVcfWriter` for `VarFishTsvWriter`.
impl AnnotatedVcfWriter for VarFishTsvWriter {
    fn set_assembly(&mut self, assembly: &Assembly) {
        self.assembly = Some(*assembly)
    }

    fn set_pedigree(&mut self, pedigree: &PedigreeByName) {
        self.pedigree = Some(pedigree.clone())
    }

    fn write_header(&mut self, header: &VcfHeader) -> Result<(), anyhow::Error> {
        self.header = Some(header.clone());
        let header = &[
            "release",
            "chromosome",
            "chromosome_no",
            "start",
            "end",
            "bin",
            "reference",
            "alternative",
            "var_type",
            "case_id",
            "set_id",
            "info",
            "genotype",
            "num_hom_alt",
            "num_hom_ref",
            "num_het",
            "num_hemi_alt",
            "num_hemi_ref",
            "in_clinvar",
            "exac_frequency",
            "exac_homozygous",
            "exac_heterozygous",
            "exac_hemizygous",
            "thousand_genomes_frequency",
            "thousand_genomes_homozygous",
            "thousand_genomes_heterozygous",
            "thousand_genomes_hemizygous",
            "gnomad_exomes_frequency",
            "gnomad_exomes_homozygous",
            "gnomad_exomes_heterozygous",
            "gnomad_exomes_hemizygous",
            "gnomad_genomes_frequency",
            "gnomad_genomes_homozygous",
            "gnomad_genomes_heterozygous",
            "gnomad_genomes_hemizygous",
            "refseq_gene_id",
            "refseq_transcript_id",
            "refseq_transcript_coding",
            "refseq_hgvs_c",
            "refseq_hgvs_p",
            "refseq_effect",
            "refseq_exon_dist",
            "ensembl_gene_id",
            "ensembl_transcript_id",
            "ensembl_transcript_coding",
            "ensembl_hgvs_c",
            "ensembl_hgvs_p",
            "ensembl_effect",
            "ensembl_exon_dist",
        ];
        writeln!(self.inner, "{}", header.join("\t"))
            .map_err(|e| anyhow::anyhow!("Error writing VarFish TSV header: {}", e))
    }

    fn write_record(&mut self, record: &VcfRecord) -> Result<(), anyhow::Error> {
        let mut tsv_record = VarFishTsvRecord::default();

        if !self.fill_coords(
            self.assembly.expect("assembly must have been set"),
            &record,
            &mut tsv_record,
        )? {
            // Record was not on canonical chromosome and should not be written out.
            return Ok(());
        }
        self.fill_genotype_and_freqs(&record, &mut tsv_record)?;
        self.fill_bg_freqs(&record, &mut tsv_record)?;
        self.fill_refseq(&record, &mut tsv_record)?;
        self.fill_ensembl(&record, &mut tsv_record)?;

        tracing::info!("writing {:?}", tsv_record);

        writeln!(self.inner, "{}", tsv_record.to_tsv().join("\t"))
            .map_err(|e| anyhow::anyhow!("Error writing VarFish TSV record: {}", e))
    }
}

/// Run the annotation with the given `Write` within the `VcfWriter`.
fn run_with_writer(writer: &mut dyn AnnotatedVcfWriter, args: &Args) -> Result<(), anyhow::Error> {
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
    writer.set_assembly(&assembly);
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

    // Load the pedigree.
    tracing::info!("Loading pedigree...");
    writer.set_pedigree(&PedigreeByName::from_path(&args.path_input_ped)?);
    tracing::info!("... done loading pedigree");

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
    let mut prev = Instant::now();
    let mut total_written = 0usize;

    writer.write_header(&header_out)?;
    let mut records = reader.records(&header_in);
    loop {
        if let Some(record) = records.next() {
            let mut vcf_record = record?;

            // TODO: ignores all but the first alternative allele!
            let vcf_var = VcfVar::from_vcf(&vcf_record);

            if prev.elapsed().as_secs() >= 60 {
                tracing::info!("at {:?}", &vcf_var);
                prev = Instant::now();
            }

            // Only attempt lookups into RocksDB for canonical contigs.
            if is_canonical(vcf_var.chrom.as_str()) {
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

                // Annotate with ClinVar information.
                annotate_record_clinvar(&db_clinvar, cf_clinvar, &key, &mut vcf_record)?;
            }

            tracing::trace!("var = {:?}", &vcf_var);
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

            // Write out the record.
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
    if let Some(path_output_vcf) = &args.output.path_output_vcf {
        if path_output_vcf.ends_with(".vcf.gz") || path_output_vcf.ends_with(".vcf.bgzf") {
            let mut writer = VcfWriter::new(
                File::create(&path_output_vcf)
                    .map(BufWriter::new)
                    .map(BgzfWriter::new)?,
            );
            run_with_writer(&mut writer, args)?;
        } else {
            let mut writer = VcfWriter::new(File::create(&path_output_vcf).map(BufWriter::new)?);
            run_with_writer(&mut writer, args)?;
        }
    } else {
        let path_output_tsv = args
            .output
            .path_output_tsv
            .as_ref()
            .expect("tsv path must be set; vcf and tsv are mutually exclusive, vcf unset");
        let mut writer = VarFishTsvWriter::with_path(path_output_tsv);
        run_with_writer(&mut writer, args)?;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use clap_verbosity_flag::Verbosity;
    use pretty_assertions::assert_eq;
    use temp_testdir::TempDir;

    use crate::annotate::seqvars::bin_from_range;

    use super::{run, Args, PathOutput};

    #[test]
    fn smoke_test_output_vcf() -> Result<(), anyhow::Error> {
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
            output: PathOutput {
                path_output_vcf: Some(path_out.into_os_string().into_string().unwrap()),
                path_output_tsv: None,
            },
            max_var_count: None,
            max_fb_tables: 5_000_000,
            path_input_ped: String::from(
                "tests/data/db/create/seqvar_freqs/db-rs1263393206/input.ped",
            ),
        };

        run(&args_common, &args)?;

        let actual = std::fs::read_to_string(args.output.path_output_vcf.unwrap())?;
        let expected = std::fs::read_to_string(
            "tests/data/db/create/seqvar_freqs/db-rs1263393206/output.vcf",
        )?;
        assert_eq!(&expected, &actual);

        Ok(())
    }

    #[test]
    fn smoke_test_output_tsv() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.tsv");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            genome_release: None,
            path_db: String::from("tests/data/annotate/db"),
            path_input_vcf: String::from(
                "tests/data/db/create/seqvar_freqs/db-rs1263393206/input.vcf",
            ),
            output: PathOutput {
                path_output_vcf: None,
                path_output_tsv: Some(path_out.into_os_string().into_string().unwrap()),
            },
            max_var_count: None,
            max_fb_tables: 5_000_000,
            path_input_ped: String::from(
                "tests/data/db/create/seqvar_freqs/db-rs1263393206/input.ped",
            ),
        };

        run(&args_common, &args)?;

        let actual = std::fs::read_to_string(args.output.path_output_tsv.unwrap())?;
        let expected = std::fs::read_to_string(
            "tests/data/db/create/seqvar_freqs/db-rs1263393206/output.tsv",
        )?;
        assert_eq!(&expected, &actual);

        Ok(())
    }

    #[test]
    fn test_bin_from_range() -> Result<(), anyhow::Error> {
        assert_eq!(bin_from_range(0, 0)?, 585);
        assert_eq!(bin_from_range(1, 0)?, 585);
        assert_eq!(bin_from_range(0, 1)?, 585);
        assert_eq!(bin_from_range(0, 42)?, 585);
        assert_eq!(bin_from_range(42_424_242, 42_424_243)?, 908);
        assert_eq!(bin_from_range(0, 42_424_243)?, 1);

        Ok(())
    }
}
