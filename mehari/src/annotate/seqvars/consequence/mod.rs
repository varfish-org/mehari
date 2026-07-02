pub mod logic;
pub mod terms;

use crate::annotate::cli::{ConsequenceBy, PredictorSettings};
use crate::annotate::seqvars;
use crate::annotate::seqvars::consequence::terms::Consequence;
use crate::annotate::seqvars::provider::{
    ConfigBuilder as MehariProviderConfigBuilder, Provider as MehariProvider,
};
use crate::db::keys::Var;
use crate::pbs::txs::TxSeqDatabase;
use anyhow::{Error, anyhow};
use enumflags2::BitFlags;
use hgvs::data::interface::Provider;
use hgvs::parser::{HgvsVariant, NoRef};
pub(crate) use logic::ConsequencePredictor;
use prost::Message;
use std::fmt;
use std::fmt::Display;
use std::fs::File;
use std::io::{Cursor, Read};
use std::path::Path;
use std::sync::Arc;
use terms::{
    ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_COMPOUND_IDS, ANN_COMPOUND_VARIANTS, ANN_TX_SEQ_ALT,
    ANN_TX_SEQ_REF, AnnField,
};

/// Load protobuf transcripts.
pub fn load_tx_db(tx_path: impl AsRef<Path> + Display) -> anyhow::Result<TxSeqDatabase> {
    // Open file and if necessary, wrap in a decompressor.
    let file = File::open(tx_path.as_ref())
        .map_err(|e| anyhow!("failed to open file {}: {}", tx_path, e))?;
    let mut reader: Box<dyn Read> = match tx_path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gz") => Box::new(flate2::read::MultiGzDecoder::new(file)),
        Some("zst") => Box::new(
            zstd::Decoder::new(file)
                .map_err(|e| anyhow!("failed to open zstd decoder for {}: {}", tx_path, e))?,
        ),
        _ => Box::new(file),
    };

    // Now read the whole file into a byte buffer.
    let mut buffer = Vec::new();
    reader
        .read_to_end(&mut buffer)
        .map_err(|e| anyhow!("failed to read file {}: {}", tx_path, e))?;

    // Deserialize the buffer with prost.
    let mut db = TxSeqDatabase::decode(&mut Cursor::new(buffer))
        .map_err(|e| anyhow!("failed to decode protobuf file {}: {}", tx_path, e));

    // back-fill deprecated fields into new fields
    if let Ok(ref mut db) = db {
        for sv in &mut db.source_version {
            if sv.assembly.is_empty() {
                sv.assembly = seqvars::_assembly_from_assembly_enum(sv.assembly_enum);
            }
            if sv.source_name.is_empty() {
                sv.source_name = seqvars::_source_from_source_enum(sv.source_name_enum);
            }
        }

        if let Some(tx_db) = &mut db.tx_db {
            for tx in &mut tx_db.transcripts {
                for aln in &mut tx.genome_alignments {
                    if aln.genome_build.is_empty() {
                        aln.genome_build = seqvars::_genome_build_from_enum(aln.genome_build_enum);
                    }
                }
            }
        }
    }

    db
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

        let mut custom_columns = Vec::new();
        if s.report_cdna_sequence.includes_ref() {
            custom_columns.push(ANN_TX_SEQ_REF.to_string());
        }
        if s.report_cdna_sequence.includes_alt() {
            custom_columns.push(ANN_TX_SEQ_ALT.to_string());
        }
        if s.report_protein_sequence.includes_ref() {
            custom_columns.push(ANN_AA_SEQ_REF.to_string());
        }
        if s.report_protein_sequence.includes_alt() {
            custom_columns.push(ANN_AA_SEQ_ALT.to_string());
        }
        if args.compound_settings.enable_compound_variants {
            custom_columns.push(ANN_COMPOUND_IDS.to_string());
            custom_columns.push(ANN_COMPOUND_VARIANTS.to_string());
        }

        let provider = Arc::new(MehariProvider::new(
            tx_db,
            reference,
            in_memory_reference,
            MehariProviderConfigBuilder::default()
                .pick_transcript(args.transcript_settings.pick_transcript.clone())
                .pick_transcript_mode(args.transcript_settings.pick_transcript_mode)
                .build()?,
        ));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .report_most_severe_consequence_by(
                    args.transcript_settings.report_most_severe_consequence_by,
                )
                .keep_intergenic(args.reporting_settings.keep_intergenic)
                .discard_utr_splice_variants(args.reporting_settings.discard_utr_splice_variants)
                .normalize(!args.do_not_normalize_variants())
                .renormalize_g(!args.do_not_renormalize_g())
                .vep_consequence_terms(args.vep_consequence_terms())
                .report_cdna_sequence(s.report_cdna_sequence)
                .report_protein_sequence(s.report_protein_sequence)
                .custom_columns(custom_columns)
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
        if let Some(mut ann_fields) = self.predictor.predict(&VcfVariant {
            chromosome: chrom.clone(),
            position: pos,
            reference: reference.clone(),
            alternative: alternative.clone(),
        })? && !ann_fields.is_empty()
        {
            let has_group_id = self
                .predictor
                .config
                .custom_columns
                .contains(&ANN_COMPOUND_IDS.to_string());
            let has_group_vars = self
                .predictor
                .config
                .custom_columns
                .contains(&ANN_COMPOUND_VARIANTS.to_string());

            if has_group_id || has_group_vars {
                let var_id = format!("{}:{}:{}:{}", chrom, pos, reference, alternative);
                for ann in &mut ann_fields {
                    if has_group_id {
                        ann.custom_fields
                            .insert(ANN_COMPOUND_IDS.to_string(), Some("Single".to_string()));
                    }
                    if has_group_vars {
                        ann.custom_fields
                            .insert(ANN_COMPOUND_VARIANTS.to_string(), Some(var_id.clone()));
                    }
                }
            }

            return Ok(ann_fields);
        }

        Ok(vec![])
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

/// Configuration for consequence prediction.
#[derive(Debug, Clone, derive_builder::Builder)]
#[builder(pattern = "immutable")]
pub struct Config {
    /// Whether to report only the worst consequence for each picked transcript.
    #[builder(default)]
    pub report_most_severe_consequence_by: Option<ConsequenceBy>,

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

    /// Whether to report cDNA sequence.
    #[builder(default)]
    pub report_cdna_sequence: SequenceReporting,

    /// Whether to report protein sequence.
    #[builder(default)]
    pub report_protein_sequence: SequenceReporting,

    /// Ordered list of extra columns registered by plugins/features.
    #[builder(default)]
    pub custom_columns: Vec<String>,
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

#[derive(
    Debug,
    Default,
    Clone,
    Copy,
    PartialEq,
    Eq,
    clap::ValueEnum,
    parse_display::FromStr,
    parse_display::Display,
)]
#[display(style = "kebab-case")]
pub enum SequenceReporting {
    #[default]
    None,
    Reference,
    Alternative,
    Both,
}

impl SequenceReporting {
    pub fn includes_ref(&self) -> bool {
        matches!(self, Self::Reference | Self::Both)
    }
    pub fn includes_alt(&self) -> bool {
        matches!(self, Self::Alternative | Self::Both)
    }
}
