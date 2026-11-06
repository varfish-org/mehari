pub mod logic;
mod pick;
pub mod plugin;
pub mod terms;

use crate::annotate::cli::{HgvsParentheses, PredictorSettings, TranscriptSettings};
use crate::annotate::seqvars;
use crate::annotate::seqvars::consequence::terms::Consequence;
use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::db::keys::Var;
use crate::errors::SeqvarsError;
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
            (pos.saturating_add(ref_len as i32)).saturating_add(logic::PADDING),
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

/// A structural variant from a symbolic allele such as `<DEL>`, which has no bases. The
/// `VcfVariant` of its record holds the padding base at `POS`, the base in front of the change.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct StructuralVariant {
    pub kind: SvKind,
    /// 1-based position of the last changed base; `POS` for an insertion.
    pub end: i32,
    /// Whether the breakpoints are uncertain: `IMPRECISE`, unless both confidence intervals are
    /// 0, or a claim of abundance only (`SVCLAIM=D`).
    pub imprecise: bool,
    /// The confidence intervals around `POS` and `end`, if given. On an imprecise allele, a
    /// missing one is unknown.
    pub intervals: (Option<ConfidenceInterval>, Option<ConfidenceInterval>),
    /// Whether the call claims that its breakpoints are joined: an `SVCLAIM` with `J`, or no
    /// `SVCLAIM` on an allele other than `<CNV>` and `<CN#>`, which claim abundance only.
    pub junction: bool,
}

/// A confidence interval around a position as offsets `(min, max)`.
pub type ConfidenceInterval = (i32, i32);

/// The type of a `StructuralVariant`. A subtype, e.g. `<DEL:ME:ALU>`, does not change it.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum SvKind {
    Deletion,
    /// `tandem` if the copy follows the original, i.e. with an adjacency claim (`SVCLAIM` with
    /// `J`), or `<DUP:TANDEM>` without `SVCLAIM`.
    Duplication {
        tandem: bool,
    },
    Inversion,
    /// Of `length` bases, if `SVLEN` is given.
    Insertion {
        length: Option<i32>,
    },
    /// `<CNV>` without a known direction.
    CopyNumberChange,
}

/// `SVCLAIM` of VCF 4.5: what the call of a structural variant claims.
#[derive(
    Debug,
    PartialEq,
    Eq,
    Clone,
    Copy,
    strum::EnumString,
    strum::Display,
    serde::Serialize,
    serde::Deserialize,
)]
#[cfg_attr(feature = "server", derive(utoipa::ToSchema))]
pub enum SvClaim {
    /// `D`: a change of the copy number.
    #[strum(serialize = "D")]
    #[serde(rename = "D")]
    Abundance,
    /// `J`: breakpoints that are joined.
    #[strum(serialize = "J")]
    #[serde(rename = "J")]
    Adjacency,
    /// `DJ`: both.
    #[strum(serialize = "DJ")]
    #[serde(rename = "DJ")]
    Both,
}

/// `SVTYPE`, which VCF 4.4 deprecates. `Other` keeps a type that callers define, e.g. `TRA` or
/// `CPX`.
#[derive(Debug, PartialEq, Eq, Clone, strum::EnumString, strum::Display)]
#[strum(serialize_all = "UPPERCASE")]
pub enum SvType {
    Del,
    Ins,
    Dup,
    Inv,
    Cnv,
    Bnd,
    #[strum(default)]
    Other(String),
}

/// The INFO fields of the record of a symbolic allele that `StructuralVariant::from_info` reads,
/// and of a breakend that `Breakend::from_info` reads.
#[derive(Debug, Default, Clone)]
pub struct SvInfo {
    /// `END`: the last base of the longest allele of the record.
    pub end: Option<i32>,
    /// `SVLEN` of the allele.
    pub sv_length: Option<i32>,
    /// `SVTYPE`, which a custom allele counts as.
    pub sv_type: Option<SvType>,
    /// `SVCLAIM`.
    pub claim: Option<SvClaim>,
    /// `CN` of a `<CNV>`.
    pub copy_number: Option<i32>,
    /// `CIPOS` as offsets `(min, max)`.
    pub cipos: Option<(i32, i32)>,
    /// `CIEND` as offsets `(min, max)`; `CIPOS` counts if it is unset.
    pub ciend: Option<(i32, i32)>,
    /// `IMPRECISE`.
    pub imprecise: bool,
}

/// The copy number `copies` of an allele, e.g. from `CN`, rounded to whole copies; saturates at
/// `i32::MAX`. A negative or non-finite number gives `SeqvarsError::InvalidVariant`.
pub fn round_copy_number(copies: f64) -> Result<i32, SeqvarsError> {
    if !copies.is_finite() || copies < 0.0 {
        return Err(SeqvarsError::InvalidVariant(format!(
            "CN {copies} is not a copy number of 0 or more"
        )));
    }
    Ok(copies.round() as i32)
}

impl StructuralVariant {
    /// The structural variant of the symbolic allele `alt` at `pos` with the INFO fields `info`,
    /// if the allele is supported, with its last base from `SVLEN` or `END`, see `span_end`. As in
    /// VCF 4.5, a `<DEL>` or `<DUP>` with a claim of abundance only (`SVCLAIM=D`) is a copy number
    /// call with uncertain breakpoints, and so is a `<CNV>` with `CN`, the copy number of the
    /// allele: with 0 copies one of a deletion, with 2 or more one of a duplication, and with 1
    /// no variant (`None`, as VCF 4.5 treats it like `GT=0`). So are `<CN0>`, `<CN1>` and
    /// `<CN2>`, as older files write them. `<CNV>` without `CN` has no direction. `<CNV:TR>`, a
    /// copy number of tandem repeats, is not supported. An `<INS>` of no bases has an unknown
    /// length. Without `CIPOS`, `POS` is precise, and `CIEND` defaults to `CIPOS`; a precise
    /// allele of an `IMPRECISE` record has intervals of 0.
    pub fn from_info(pos: i32, alt: &str, info: &SvInfo) -> Option<Self> {
        let (base, subtype) = symbolic_type(alt, info.sv_type.as_ref())?;
        let abundance_only = info.claim == Some(SvClaim::Abundance);
        let copy_number = |copies: Option<i32>| match copies {
            Some(0) => SvKind::Deletion,
            Some(copies) if copies >= 2 => SvKind::Duplication { tandem: false },
            _ => SvKind::CopyNumberChange,
        };
        let (kind, abundance_only) = match base {
            "CNV" if subtype == Some("TR") => return None,
            "CNV" if info.copy_number == Some(1) => return None,
            "DEL" => (SvKind::Deletion, abundance_only),
            "DUP" => {
                let tandem = match info.claim {
                    Some(SvClaim::Adjacency | SvClaim::Both) => true,
                    Some(SvClaim::Abundance) => false,
                    None => subtype == Some("TANDEM"),
                };
                (SvKind::Duplication { tandem }, abundance_only)
            }
            "INV" => (SvKind::Inversion, false),
            "INS" => (
                SvKind::Insertion {
                    length: (info.sv_length.map(i32::saturating_abs)).filter(|&length| length > 0),
                },
                false,
            ),
            "CNV" => (copy_number(info.copy_number), true),
            _ => {
                let copies = base.strip_prefix("CN")?.parse().ok()?;
                if copies == 1 {
                    return None;
                }
                (copy_number(Some(copies)), true)
            }
        };
        let end = span_end(pos, alt, info.sv_type.as_ref(), info.sv_length, info.end);
        if end == pos && !matches!(kind, SvKind::Insertion { .. }) {
            return None;
        }
        // `CIEND` defaults to `CIPOS`. As in VCF 4.5, `IMPRECISE` makes the breakpoints
        // uncertain, unless both intervals are 0, as for a precise allele of the record, and a
        // missing `CIPOS` is unknown. Without `IMPRECISE`, the intervals give the micro-homology
        // of a precise call.
        let intervals = (info.cipos, info.ciend.or(info.cipos));
        let imprecise =
            abundance_only || info.imprecise && intervals != (Some((0, 0)), Some((0, 0)));
        let junction = match info.claim {
            Some(SvClaim::Adjacency | SvClaim::Both) => true,
            Some(SvClaim::Abundance) => false,
            None => matches!(base, "DEL" | "DUP" | "INV" | "INS"),
        };
        Some(Self {
            kind,
            end,
            imprecise,
            intervals,
            junction,
        })
    }

    /// Whether the symbolic allele `alt` with the INFO fields `info` keeps the copy number of the
    /// reference, so that `from_info` gives `None` as for no variant: a `<CNV>` with a `CN` of 1,
    /// or a `<CN1>`.
    pub fn keeps_copy_number(alt: &str, info: &SvInfo) -> bool {
        match symbolic_type(alt, info.sv_type.as_ref()) {
            Some(("CNV", subtype)) => subtype != Some("TR") && info.copy_number == Some(1),
            Some((base, _)) => base.strip_prefix("CN").and_then(|n| n.parse().ok()) == Some(1),
            None => false,
        }
    }
}

/// The type of the symbolic allele `alt` and its subtype, e.g. `DEL` and `ME` for
/// `<DEL:ME:ALU>`. A custom allele, which `##ALT` declares, counts as its `SVTYPE`, `sv_type`, if
/// that is the type of a symbolic allele.
fn symbolic_type<'a>(alt: &'a str, sv_type: Option<&SvType>) -> Option<(&'a str, Option<&'a str>)> {
    let mut parts = alt.strip_prefix('<')?.strip_suffix('>')?.split(':');
    let (base, subtype) = (parts.next()?, parts.next());
    let base = match base {
        "DEL" | "DUP" | "INV" | "INS" | "CNV" => base,
        _ if base
            .strip_prefix("CN")
            .is_some_and(|n| !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit())) =>
        {
            base
        }
        _ => match sv_type? {
            SvType::Del => "DEL",
            SvType::Ins => "INS",
            SvType::Dup => "DUP",
            SvType::Inv => "INV",
            SvType::Cnv => "CNV",
            SvType::Bnd | SvType::Other(_) => return None,
        },
    };
    Some((base, subtype))
}

/// The last base of the symbolic allele `alt` at `pos` with the INFO fields `SVTYPE`, `SVLEN`
/// and `END`, see `compound::event_end`. An insertion lies at `pos`. Else the own `SVLEN` of the
/// allele decides, as in VCF 4.5, then `END`. An `SVLEN` of 0 is the length difference of VCF
/// 4.3 and earlier, e.g. of an `<INV>`, so `END` decides.
pub fn span_end(
    pos: i32,
    alt: &str,
    sv_type: Option<&SvType>,
    sv_length: Option<i32>,
    end: Option<i32>,
) -> i32 {
    if symbolic_type(alt, sv_type).is_some_and(|(base, _)| base == "INS") {
        return pos;
    }
    match (sv_length, end) {
        (Some(len), _) if len != 0 => pos.saturating_add(len.saturating_abs()),
        (_, Some(end)) => pos.max(end),
        (_, None) => pos,
    }
}

/// A breakend from its ALT, e.g. `G]17:198982]` (VCF 4.5): the base at `POS` joins a piece of
/// the mate contig, with the bases `inserted` between them.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Breakend {
    /// Whether `POS` and the bases in front of it stay, as in `t[p[` and `t]p]`, else `POS` and
    /// the bases behind it, as in `]p]t` and `[p[t`.
    pub keeps_left: bool,
    /// The bases of t besides the one at `POS`, in the order of the joined chromosome.
    pub inserted: String,
    /// `None` for a single breakend, e.g. `G.`.
    pub mate: Option<Mate>,
    /// Whether `POS` is uncertain: with `IMPRECISE` and a `CIPOS` other than `0,0` or none. The
    /// INFO of the record tells, not the ALT, so `from_str` gives `false`, and `from_info`
    /// reads it.
    pub imprecise: bool,
    /// The confidence interval around `POS` as offsets `(min, max)`, if given.
    pub cipos: Option<(i32, i32)>,
}

impl Breakend {
    /// The breakend of `alt` with the uncertainty of `POS` from the INFO fields `info`. As in
    /// VCF 4.5, `IMPRECISE` makes `POS` uncertain unless `CIPOS` is `0,0`, and a missing `CIPOS`
    /// leaves the interval open; `CIPOS` without `IMPRECISE` is the micro-homology of a precise
    /// call.
    pub fn from_info(alt: &str, info: &SvInfo) -> Result<Self, Error> {
        let mut breakend: Self = alt.parse()?;
        breakend.cipos = info.cipos;
        breakend.imprecise = info.imprecise && info.cipos != Some((0, 0));
        Ok(breakend)
    }
}

/// The mate of a `Breakend`.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Mate {
    pub chromosome: String,
    /// 1-based; 0 and the contig length + 1 are the telomeres.
    pub position: i32,
    /// Whether the joined piece runs from `position` to the end of the contig, as `[` says, else
    /// from the start of the contig to `position`.
    pub extends_right: bool,
}

/// Two transcripts of different genes that a junction of a structural variant joins, see
/// `ConsequencePredictor::sv_fusions`.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Fusion {
    /// The transcript whose 5' part the junction keeps and the one whose 3' part follows it, or
    /// the transcripts before and behind a head-to-head or tail-to-tail junction on the joined
    /// chromosome, as a bidirectional fusion has no 5' gene.
    pub transcripts: (String, String),
    /// The genes in the same order, as `GENE1::GENE2`, by symbol or else by ID.
    pub label: String,
    /// Whether the genes point into or away from each other on the joined chromosome, i.e. lie
    /// on alternate strands of it, else they read in one direction.
    pub bidirectional: bool,
    /// Whether a deletion makes the junction.
    pub deletion: bool,
    /// Whether a unidirectional fusion keeps the reading frame of the 5' gene in the 3' gene,
    /// see `ConsequencePredictor::fusion_frame`; `None` where that is unknown.
    pub in_frame: Option<bool>,
    /// Whether `in_frame` assumes normal splicing, as a junction lies in an intron.
    pub assumes_splicing: bool,
}

/// The annotations of an allele, see `ConsequencePredictor::predict_variant`.
#[derive(Debug, Default, Clone)]
pub struct Annotations {
    /// One per reported transcript.
    pub anns: Vec<terms::AnnField>,
    /// The gene fusions of a structural variant or a breakend, sorted, see `Fusion::label`.
    pub fusions: Vec<String>,
}

/// The bracket forms `t[p[`, `t]p]`, `]p]t` and `[p[t`, with `.` for t at a telomere, and the
/// single breakends `t.` and `.t`.
impl std::str::FromStr for Breakend {
    type Err = anyhow::Error;

    fn from_str(alt: &str) -> Result<Self, Self::Err> {
        let (t, keeps_left, mate) = match alt.find(['[', ']']) {
            Some(open) => {
                let bracket = &alt[open..=open];
                let (mate, tail) = (alt[open + 1..].split_once(bracket))
                    .ok_or_else(|| anyhow!("unclosed bracket in breakend {alt:?}"))?;
                let (t, keeps_left) = match (&alt[..open], tail) {
                    (head, "") => (head, true),
                    ("", tail) => (tail, false),
                    _ => anyhow::bail!("bases on both sides of breakend {alt:?}"),
                };
                let (chromosome, position) = (mate.rsplit_once(':'))
                    .filter(|(chromosome, _)| !chromosome.is_empty())
                    .ok_or_else(|| anyhow!("no mate contig in breakend {alt:?}"))?;
                let position = (position.parse::<i32>().ok())
                    .filter(|position| *position >= 0)
                    .ok_or_else(|| anyhow!("invalid mate position in breakend {alt:?}"))?;
                let mate = Mate {
                    chromosome: chromosome.into(),
                    position,
                    extends_right: bracket == "[",
                };
                (t, keeps_left, Some(mate))
            }
            None => match (alt.strip_suffix('.'), alt.strip_prefix('.')) {
                (Some(t), _) => (t, true, None),
                (_, Some(t)) => (t, false, None),
                _ => anyhow::bail!("not a breakend: {alt:?}"),
            },
        };
        // The base at `POS` ends t on the side away from the mate; t is `.` at a telomere.
        let inserted = if t == "." && mate.is_some() {
            ""
        } else {
            anyhow::ensure!(
                !t.is_empty() && t.bytes().all(|base| base.is_ascii_alphabetic()),
                "invalid bases in breakend {alt:?}"
            );
            if keeps_left {
                &t[1..]
            } else {
                &t[..t.len() - 1]
            }
        };
        Ok(Self {
            keeps_left,
            inserted: inserted.to_ascii_uppercase(),
            mate,
            imprecise: false,
            cipos: None,
        })
    }
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

/// A variant with its allele classified once, see `Variant::from_info`.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Variant {
    pub var: VcfVariant,
    pub kind: VariantKind,
}

/// The allele of a `Variant`.
#[derive(Debug, PartialEq, Eq, Clone)]
pub enum VariantKind {
    /// Bases, or none for a deletion.
    Sequence,
    /// A symbolic allele such as `<DEL>`.
    Structural(StructuralVariant),
    /// A breakend such as `G]17:198982]`.
    Breakend(Breakend),
}

impl Variant {
    /// `var` with its allele classified, with `info` from the INFO fields of its record: a
    /// symbolic allele as a structural variant, see `StructuralVariant::from_info`, a breakend
    /// as such, see `Breakend::from_info`, and bases as a sequence allele. `*`, `<*>`,
    /// `<NON_REF>` and `ALT` `.` are no variants of their own (`None`). A symbolic allele that
    /// keeps the copy number of the reference is no variant either, see
    /// `StructuralVariant::keeps_copy_number`. A symbolic allele that
    /// `StructuralVariant::from_info` rejects otherwise, a malformed breakend, a type other than
    /// those of `SVTYPE` for a custom allele, which alone takes it, see `symbolic_type`, and a
    /// confidence interval with its minimum above its maximum give `SeqvarsError::InvalidVariant`.
    pub fn from_info(var: VcfVariant, info: &SvInfo) -> Result<Option<Self>, SeqvarsError> {
        let invalid = |reason: String| Err(SeqvarsError::InvalidVariant(reason));
        let alt = &var.alternative;
        if let Some(SvType::Other(sv_type)) = &info.sv_type
            && alt.starts_with('<')
            && symbolic_type(alt, None).is_none()
        {
            return invalid(format!(
                "SVTYPE {sv_type:?} is none of DEL, INS, DUP, INV, CNV and BND"
            ));
        }
        for (key, interval) in [("CIPOS", info.cipos), ("CIEND", info.ciend)] {
            if let Some((min, max)) = interval
                && min > max
            {
                return invalid(format!(
                    "{key} {min},{max} has its minimum above its maximum"
                ));
            }
        }
        let kind = if logic::is_sequence_allele(alt) {
            VariantKind::Sequence
        } else if matches!(alt.as_str(), "*" | "<*>" | "<NON_REF>" | ".") {
            return Ok(None);
        } else if alt.starts_with('<') {
            if let Some(sv) = StructuralVariant::from_info(var.position, alt, info) {
                VariantKind::Structural(sv)
            } else if StructuralVariant::keeps_copy_number(alt, info) {
                return Ok(None);
            } else {
                // With a length, the allele passes, so its end is what is missing.
                let with_length = SvInfo {
                    sv_length: Some(1),
                    ..info.clone()
                };
                return invalid(
                    if StructuralVariant::from_info(var.position, alt, &with_length).is_some() {
                        format!("{alt} has no end past its position (END or SVLEN)")
                    } else {
                        format!("{alt} is not supported")
                    },
                );
            }
        } else {
            match Breakend::from_info(alt, info) {
                Ok(breakend) => VariantKind::Breakend(breakend),
                Err(e) => return invalid(e.to_string()),
            }
        };
        Ok(Some(Self { var, kind }))
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

    /// Drop the splice terms of variants in UTRs.
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
    /// of indels. Implies `vep_consequence_terms`.
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
    use super::{Breakend, Mate, SvClaim, load_tx_db};
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

    /// The breakend of an ALT: the side that stays, the bases of t besides the one at `POS`, and
    /// the mate as `(chromosome, position, extends_right)`. The examples are those of VCF 4.5.
    #[rstest::rstest]
    #[case::joins_right("G[13:123457[", true, "", Some(("13", 123457, true)))]
    #[case::joins_left_reversed("G]13:123457]", true, "", Some(("13", 123457, false)))]
    #[case::left_joins("]13:123457]G", false, "", Some(("13", 123457, false)))]
    #[case::right_joins_reversed("[13:123457[G", false, "", Some(("13", 123457, true)))]
    #[case::inserted_left("CAGTNNNNNCA[2:321682[", true, "AGTNNNNNCA", Some(("2", 321682, true)))]
    #[case::inserted_right("]13:123456]AGTNNNNNCAT", false, "AGTNNNNNCA", Some(("13", 123456, false)))]
    #[case::telomere(".[13:123457[", true, "", Some(("13", 123457, true)))]
    #[case::mate_at_telomere("]1:0]A", false, "", Some(("1", 0, false)))]
    #[case::mate_contig_with_colon("A[HLA-A*01:01:01:01:100[", true, "", Some(("HLA-A*01:01:01:01", 100, true)))]
    #[case::single_left("G.", true, "", None)]
    #[case::single_right(".A", false, "", None)]
    #[case::single_inserted("GAC.", true, "AC", None)]
    #[case::lower_case("gac[2:321682[", true, "AC", Some(("2", 321682, true)))]
    #[case::lower_case_right("]13:123456]acG", false, "AC", Some(("13", 123456, false)))]
    fn breakend_from_alt(
        #[case] alt: &str,
        #[case] keeps_left: bool,
        #[case] inserted: &str,
        #[case] mate: Option<(&str, i32, bool)>,
    ) -> Result<(), anyhow::Error> {
        let mate = mate.map(|(chromosome, position, extends_right)| Mate {
            chromosome: chromosome.into(),
            position,
            extends_right,
        });
        let expected = Breakend {
            keeps_left,
            inserted: inserted.into(),
            mate,
            imprecise: false,
            cipos: None,
        };
        assert_eq!(alt.parse::<Breakend>()?, expected);
        Ok(())
    }

    /// A malformed breakend, or another allele, is an error, not a panic.
    #[rstest::rstest]
    #[case::bases_on_both_sides("A]1:100]C")]
    #[case::mixed_brackets("A[1:100]")]
    #[case::unclosed_bracket("A[1:100")]
    #[case::no_bases("[1:100[")]
    #[case::no_position("A[1[")]
    #[case::no_contig("A[:100[")]
    #[case::negative_position("A[1:-1[")]
    #[case::invalid_bases("A-[1:100[")]
    #[case::telomere_without_mate(".")]
    #[case::sequence("ACGT")]
    #[case::symbolic("<DEL>")]
    #[case::star("*")]
    fn malformed_breakend(#[case] alt: &str) {
        assert!(alt.parse::<Breakend>().is_err());
    }

    /// A call claims that its breakpoints are joined with `J` in `SVCLAIM`, or without
    /// `SVCLAIM` unless it is a copy number allele.
    #[rstest::rstest]
    #[case::deletion("<DEL>", None, None, true)]
    #[case::abundance_claim("<DEL>", Some(SvClaim::Abundance), None, false)]
    #[case::both_claims("<DEL>", Some(SvClaim::Both), None, true)]
    #[case::inversion("<INV>", None, None, true)]
    #[case::cn0("<CN0>", None, None, false)]
    #[case::cnv_lost("<CNV>", None, Some(0), false)]
    fn junction_claim(
        #[case] alt: &str,
        #[case] claim: Option<SvClaim>,
        #[case] copy_number: Option<i32>,
        #[case] expected: bool,
    ) {
        let info = super::SvInfo {
            end: Some(2000),
            claim,
            copy_number,
            ..Default::default()
        };
        let sv = super::StructuralVariant::from_info(1000, alt, &info);
        assert_eq!(sv.map(|sv| sv.junction), Some(expected));
    }
}
