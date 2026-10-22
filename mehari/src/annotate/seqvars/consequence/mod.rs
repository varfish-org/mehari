pub mod logic;
mod pick;
pub mod plugin;
pub mod terms;

use crate::annotate::cli::{HgvsParentheses, PredictorSettings, TranscriptSettings};
use crate::annotate::seqvars;
use crate::annotate::seqvars::consequence::terms::Consequence;
use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::db::keys::Var;
use crate::pbs::txs::TxSeqDatabase;
use anyhow::{Error, anyhow};
use enumflags2::BitFlags;
use hgvs::data::interface::Provider;
use hgvs::parser::{HgvsVariant, NoRef};
use itertools::Itertools;
pub(crate) use logic::ConsequencePredictor;
use plugin::{CompoundPlugin, RunInfo, TranscriptPlugin};
use prost::Message;
use std::fmt;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read};
use std::path::Path;
use std::sync::Arc;
use terms::{ANN_PICK, AnnField};

/// Load protobuf transcripts.
pub fn load_tx_db(tx_path: impl AsRef<Path> + Display) -> anyhow::Result<TxSeqDatabase> {
    // Open file and if necessary, wrap in a decompressor.
    let file = File::open(tx_path.as_ref())
        .map_err(|e| anyhow!("failed to open file {}: {}", tx_path, e))?;
    // The decompressed size, if the zstd frame header records it.
    let mut content_size = None;
    let mut reader: Box<dyn Read> = match tx_path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gz") => Box::new(flate2::read::MultiGzDecoder::new(file)),
        Some("zst") => {
            let mut file = BufReader::with_capacity(zstd::zstd_safe::DCtx::in_size(), file);
            let header = file
                .fill_buf()
                .map_err(|e| anyhow!("failed to read file {}: {}", tx_path, e))?;
            content_size = zstd::zstd_safe::get_frame_content_size(header)
                .ok()
                .flatten();
            Box::new(
                zstd::Decoder::with_buffer(file)
                    .map_err(|e| anyhow!("failed to open zstd decoder for {}: {}", tx_path, e))?,
            )
        }
        _ => Box::new(file),
    };

    // Now read the whole file into a byte buffer.  Reserving the decompressed size up front
    // avoids growing, and thus copying, the buffer while reading.
    let mut buffer = Vec::new();
    if let Some(size) = content_size.and_then(|size| usize::try_from(size).ok()) {
        // If this fails, e.g., for a corrupt header, the buffer grows while reading instead.
        let _ = buffer.try_reserve_exact(size);
    }
    reader
        .read_to_end(&mut buffer)
        .map_err(|e| anyhow!("failed to read file {}: {}", tx_path, e))?;

    // Deserialize the buffer with prost.
    let mut db = TxSeqDatabase::decode(&mut Cursor::new(buffer))
        .map_err(|e| anyhow!("failed to decode protobuf file {}: {}", tx_path, e));

    // back-fill deprecated fields into new fields
    if let Ok(ref mut db) = db {
        for sv in &mut db.source_version {
            #[allow(deprecated)]
            if sv.assembly.is_empty() {
                sv.assembly = seqvars::_assembly_from_assembly_enum(sv.assembly_enum);
            }
            #[allow(deprecated)]
            if sv.source_name.is_empty() {
                sv.source_name = seqvars::_source_from_source_enum(sv.source_name_enum);
            }
        }

        if let Some(tx_db) = &mut db.tx_db {
            for tx in &mut tx_db.transcripts {
                for aln in &mut tx.genome_alignments {
                    #[allow(deprecated)]
                    if aln.genome_build.is_empty() {
                        aln.genome_build = seqvars::_genome_build_from_enum(aln.genome_build_enum);
                    }
                }
            }
        }
    }

    let db = db?;
    refuse_without_gene_ids(&db, &tx_path)?;
    Ok(db)
}

/// Fail if the database predates annotation gene IDs, i.e., if it has transcripts and none of
/// them has a `gene_id`.
fn refuse_without_gene_ids(db: &TxSeqDatabase, tx_path: impl Display) -> anyhow::Result<()> {
    let predates_gene_ids = db.tx_db.as_ref().is_some_and(|tx_db| {
        !tx_db.transcripts.is_empty() && tx_db.transcripts.iter().all(|tx| tx.gene_id.is_empty())
    });
    if !predates_gene_ids {
        return Ok(());
    }

    let mehari_versions = db
        .source_version
        .iter()
        .map(|sv| sv.mehari_version.as_str())
        .filter(|version| !version.is_empty())
        .unique()
        .join(", ");
    let built_by = if mehari_versions.is_empty() {
        "an unknown mehari version".to_string()
    } else {
        format!("mehari {mehari_versions}")
    };
    Err(anyhow!(
        "transcript database {} predates annotation gene IDs: it was built by {}. Rebuild it \
        with mehari {} or newer.",
        tx_path,
        built_by,
        crate::common::version()
    ))
}

pub(crate) fn load_transcript_dbs_for_assembly(
    tx_sources: &[String],
    assembly: &str,
) -> Result<Vec<TxSeqDatabase>, Error> {
    let req_assembly = assembly;

    // Filter out any transcript databases that do not match the requested assembly.
    let check_assembly = |db: &TxSeqDatabase, required: &str| {
        db.source_version
            .iter()
            .flat_map(|s| {
                vec![
                    s.assembly.clone(),
                    #[allow(deprecated)]
                    seqvars::_assembly_from_assembly_enum(s.assembly_enum),
                ]
            })
            .any(|a| a.eq_ignore_ascii_case(required))
    };
    let databases = tx_sources
        .iter()
        .enumerate()
        .map(|(i, path)| (i, load_tx_db(path)))
        .filter_map(|(i, txdb)| match txdb {
            Ok(db) => {
                if check_assembly(&db, req_assembly) {
                    Some(Ok(db))
                } else {
                    tracing::info!("Skipping transcript database {} as its version {:?} does not support the requested assembly ({:?})", &tx_sources[i], &db.source_version, &assembly);
                    None
                }
            },
            Err(_) => Some(txdb),
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    Ok(databases)
}

pub(crate) struct ConsequenceAnnotator {
    pub(crate) predictor: ConsequencePredictor,
}

impl ConsequenceAnnotator {
    fn new(predictor: ConsequencePredictor) -> Self {
        Self { predictor }
    }

    pub(crate) fn from_db_and_settings(
        tx_db: TxSeqDatabase,
        reference: Option<impl AsRef<Path>>,
        in_memory_reference: bool,
        predictor_settings: &PredictorSettings,
    ) -> anyhow::Result<Self> {
        let args = predictor_settings;
        let s = &args.reporting_settings;

        let run = RunInfo::new(&tx_db, args.vep_compatibility_mode);
        let mut plugins = s
            .plugins
            .iter()
            .map(|plugin| plugin.build(&run))
            .collect::<Result<Vec<_>, _>>()?;
        if args.compound_settings.enable_compound_variants {
            plugins.push(Arc::new(CompoundPlugin::default()));
        }

        let provider = Arc::new(MehariProvider::new(tx_db, reference, in_memory_reference));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .transcript_settings(args.transcript_settings.clone())
                .keep_intergenic(args.reporting_settings.keep_intergenic)
                .discard_utr_splice_variants(args.reporting_settings.discard_utr_splice_variants)
                .normalize(!args.do_not_normalize_variants())
                .renormalize_g(!args.do_not_renormalize_g())
                .vep_consequence_terms(args.vep_consequence_terms())
                .vep_compatibility_mode(args.vep_compatibility_mode)
                .hgvs_parentheses(args.hgvs_parentheses())
                .plugins(plugins)
                .build()?,
        );
        Ok(Self::new(predictor))
    }

    pub(crate) fn annotate(&self, vcf_var: &Var) -> anyhow::Result<Vec<AnnField>> {
        let Var {
            chrom,
            pos,
            reference,
            alternative,
        } = vcf_var.clone();

        // Annotate with variant effect.
        Ok(self
            .predictor
            .predict(&VcfVariant {
                chromosome: chrom,
                position: pos,
                reference,
                alternative,
            })?
            .unwrap_or_default())
    }
}

/// Helper to find the maximum genomic bounds of all transcripts overlapping a variant.
pub fn get_transcript_boundaries(
    predictor: &ConsequencePredictor,
    chrom: &str,
    pos: i32,
    ref_len: usize,
) -> (std::collections::HashSet<String>, i32, i32) {
    let mut accessions = std::collections::HashSet::new();
    let mut min_start = i32::MAX;
    let mut max_end = -1;

    let txs = predictor
        .provider
        .get_tx_for_region(
            chrom,
            logic::ALT_ALN_METHOD,
            pos - logic::PADDING,
            pos + ref_len as i32 + logic::PADDING,
        )
        .unwrap_or_default();

    for tx_record in txs {
        accessions.insert(tx_record.tx_ac.clone());
        if let Some(tx_details) = predictor.provider.get_tx(&tx_record.tx_ac)
            && let Some(aln) = tx_details.genome_alignments.first()
        {
            for exon in &aln.exons {
                min_start = std::cmp::min(min_start, exon.alt_start_i);
                max_end = std::cmp::max(max_end, exon.alt_end_i);
            }
        }
    }

    (accessions, min_start, max_end)
}

/// A variant description how VCF would do it.
#[derive(Debug, PartialEq, Eq, Clone, Default)]
pub struct VcfVariant {
    /// Chromosome name.
    pub chromosome: String,
    /// 1-based position on the chromosome of first base of `reference`.
    pub position: i32,
    /// Reference bases.
    pub reference: String,
    /// Alternative bases.
    pub alternative: String,
}

/// `chromosome:position:reference:alternative`, e.g. `17:41215920:G:T`.
impl Display for VcfVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{}:{}:{}",
            self.chromosome, self.position, self.reference, self.alternative
        )
    }
}

impl std::str::FromStr for VcfVariant {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let [chromosome, position, reference, alternative] = s.split(':').collect::<Vec<_>>()[..]
        else {
            anyhow::bail!("expected chromosome:position:reference:alternative, got {s:?}");
        };
        Ok(Self {
            chromosome: chromosome.into(),
            position: position.parse()?,
            reference: reference.into(),
            alternative: alternative.into(),
        })
    }
}

/// Configuration for consequence prediction.
#[derive(Debug, Clone, derive_builder::Builder)]
#[builder(pattern = "immutable", build_fn(validate = "Self::validate"))]
pub struct Config {
    /// Which transcripts to report.
    #[builder(default)]
    pub transcript_settings: TranscriptSettings,

    /// Whether to keep intergenic variants.
    #[builder(default = "false")]
    pub keep_intergenic: bool,

    /// Whether to report splice variants in UTRs.
    #[builder(default = "false")]
    pub discard_utr_splice_variants: bool,

    /// Whether to normalize HGVS variants.
    #[builder(default = "true")]
    pub normalize: bool,

    /// Whether re-normalize genomic variants.
    #[builder(default = "true")]
    pub renormalize_g: bool,

    /// Whether to report VEP consequence terms.
    #[builder(default = "false")]
    pub vep_consequence_terms: bool,

    /// Whether to reproduce VEP's behaviour where it differs from HGVS or SO, e.g. its stop terms
    /// of indels. `--vep-compatibility-mode` also sets `vep_consequence_terms`.
    #[builder(default = "false")]
    pub vep_compatibility_mode: bool,

    /// Where to write parentheses, e.g. `p.(Gly2Arg)` for a predicted change.
    #[builder(default)]
    pub hgvs_parentheses: HgvsParentheses,

    /// Plugins filling the custom ANN columns, in column order.
    #[builder(default)]
    pub plugins: Vec<Arc<dyn TranscriptPlugin>>,
}

impl ConfigBuilder {
    fn validate(&self) -> Result<(), String> {
        let flag_pick = self
            .transcript_settings
            .as_ref()
            .is_some_and(|settings| settings.flag_pick);
        let columns = self
            .plugins
            .iter()
            .flatten()
            .flat_map(|plugin| plugin.columns())
            .chain(flag_pick.then(|| ANN_PICK.to_string()));
        match columns.duplicates().next() {
            Some(column) => Err(format!(
                "column {column:?} is filled by more than one plugin, or by a plugin and \
                 --flag-pick"
            )),
            None => Ok(()),
        }
    }
}

impl Default for Config {
    fn default() -> Self {
        ConfigBuilder::default().build().unwrap()
    }
}

pub type Consequences = BitFlags<Consequence>;

struct FormattedLoc<'a>(&'a HgvsVariant);

impl<'a> fmt::Display for FormattedLoc<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            HgvsVariant::CdsVariant { loc_edit, .. } => write!(f, "c.{}", NoRef(loc_edit)),
            HgvsVariant::GenomeVariant { loc_edit, .. } => write!(f, "g.{}", NoRef(loc_edit)),
            HgvsVariant::MtVariant { loc_edit, .. } => write!(f, "m.{}", NoRef(loc_edit)),
            HgvsVariant::TxVariant { loc_edit, .. } => write!(f, "n.{}", NoRef(loc_edit)),
            HgvsVariant::ProtVariant { loc_edit, .. } => write!(f, "p.{}", NoRef(loc_edit)),
            HgvsVariant::RnaVariant { loc_edit, .. } => write!(f, "r.{}", NoRef(loc_edit)),
        }
    }
}

#[cfg(test)]
mod test {
    use super::load_tx_db;
    use crate::pbs::txs::{SourceVersion, Transcript, TranscriptDb, TxSeqDatabase};
    use prost::Message;
    use temp_testdir::TempDir;

    /// Write a database with transcripts of the given `gene_id`s, built by the given mehari
    /// versions, and load it.
    fn write_and_load(
        gene_ids: &[&str],
        mehari_versions: &[&str],
    ) -> anyhow::Result<TxSeqDatabase> {
        let db = TxSeqDatabase {
            tx_db: Some(TranscriptDb {
                transcripts: gene_ids
                    .iter()
                    .map(|&gene_id| Transcript {
                        gene_id: gene_id.into(),
                        ..Default::default()
                    })
                    .collect(),
                ..Default::default()
            }),
            source_version: mehari_versions
                .iter()
                .map(|&mehari_version| SourceVersion {
                    mehari_version: mehari_version.into(),
                    ..Default::default()
                })
                .collect(),
            ..Default::default()
        };
        let temp = TempDir::default();
        let path = temp.join("txs.bin");
        std::fs::write(&path, db.encode_to_vec())?;
        load_tx_db(path.display().to_string())
    }

    /// A database whose transcripts have no `gene_id` predates annotation gene IDs. Loading it
    /// fails and names each mehari version that built it once.
    #[rstest::rstest]
    #[case::merged_database(&["0.38.0", "0.37.1", "0.38.0"], "built by mehari 0.38.0, 0.37.1.")]
    #[case::no_mehari_version(&[], "built by an unknown mehari version.")]
    fn old_database_is_refused(#[case] mehari_versions: &[&str], #[case] built_by: &str) {
        let err = write_and_load(&["", ""], mehari_versions).unwrap_err();
        let msg = err.to_string();

        assert!(
            msg.contains("txs.bin predates annotation gene IDs"),
            "{msg}"
        );
        assert!(msg.contains(built_by), "{msg}");
        assert!(
            msg.ends_with("Rebuild it with mehari 0.0.0 or newer."),
            "{msg}"
        );
    }

    /// A database loads if a transcript has a `gene_id`, or if it has no transcripts.
    #[rstest::rstest]
    #[case::new_database(&["672", ""])]
    #[case::no_transcripts(&[])]
    fn database_loads(#[case] gene_ids: &[&str]) {
        assert!(write_and_load(gene_ids, &["0.38.0"]).is_ok());
    }
}
