//! Compute molecular consequence of variants.
use super::pick::{self, Candidate};
use super::plugin::{TxAnnotation, TxContext};
use super::terms::{ANN_PICK, FeatureTag, GroupedAlleles};
use super::terms::{
    Allele, AnnField, Consequence, FeatureBiotype, FeatureType, Pos, Rank, SoFeature,
};
use crate::annotate::seqvars::consequence::{Config, Consequences, FormattedLoc, VcfVariant};
use crate::annotate::seqvars::provider::PbsTranscriptExt;
use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::errors::{GroupValidationError, SeqvarsError};
use crate::pbs::txs::{GenomeAlignment, Strand, Transcript, TranscriptBiotype, TranscriptTag};
use hgvs::mapper::altseq::{AltSeqBuilder, RefTranscriptData, VariantLocation};
use hgvs::parser::{NoRef, ProteinEdit};
use hgvs::sequences::translate_cds_with_exceptions;
use hgvs::{
    data::interface::{Provider, TxForRegionRecord},
    mapper::{Error, assembly},
    parser::{
        Accession, CdsFrom, CdsPos, GenomeInterval, GenomeLocEdit, HgvsVariant, MtInterval,
        MtLocEdit, Mu, NaEdit, ProtLocEdit,
    },
};
use itertools::Itertools;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

/// Wrap mapper, provider, and map for consequence prediction.
#[derive(derivative::Derivative)]
#[derivative(Debug)]
pub struct ConsequencePredictor {
    /// The internal transcript provider for locating transcripts.
    #[derivative(Debug = "ignore")]
    pub(crate) provider: Arc<MehariProvider>,

    /// Assembly mapper for variant consequence prediction.
    #[derivative(Debug = "ignore")]
    mapper: assembly::Mapper,

    /// Configuration for the predictor.
    #[derivative(Debug = "ignore")]
    pub(crate) config: Config,

    /// The columns of `config.plugins`, in order, then `pick` for `--flag-pick`.
    custom_columns: Vec<String>,
    /// Gene annotations dropped by `--restrict-transcripts`, for the summary at the end.
    restrict_dropped_genes: AtomicUsize,

    /// Variants with gene annotations dropped by `--restrict-transcripts`.
    restrict_dropped_variants: AtomicUsize,
}

/// Padding to look for genes upstream/downstream.
pub const PADDING: i32 = 5_000;
/// Generally used alternative alignment method.
pub const ALT_ALN_METHOD: &str = "splign";

#[derive(Debug, Clone)]
struct HgvsProjectionContext {
    n: Option<HgvsVariant>,
    c: Option<HgvsVariant>,
    p: Option<HgvsVariant>,
    /// The region of the transcript that `c` changes, as hgvs-rs uses it to project `c` to `p`.
    region: Option<VariantLocation>,
    /// With VEP terms, the c. variant at the input position, before the 3' shift. VEP checks
    /// with it whether an indel overlaps the stop codon.
    unshifted_c: Option<HgvsVariant>,
}

impl HgvsProjectionContext {
    /// Check if the variant falls within the boundaries of the CDS (start to stop),
    /// including intronic parts.
    fn is_within_cds_bounds(&self) -> bool {
        if let Some(HgvsVariant::CdsVariant { loc_edit, .. }) = &self.c {
            let loc = loc_edit.loc.inner();
            let start = &loc.start;
            let end = &loc.end;

            let is_5_prime = start.cds_from == CdsFrom::Start
                && start.base < 0
                && end.cds_from == CdsFrom::Start
                && end.base < 0;

            let is_3_prime = start.cds_from == CdsFrom::End && end.cds_from == CdsFrom::End;

            !is_5_prime && !is_3_prime
        } else {
            false
        }
    }

    #[allow(dead_code)]
    /// Check if the variant is strictly within the coding sequence, i.e., not in UTRs or intronic.
    fn is_within_coding_sequence(&self) -> bool {
        if let Some(HgvsVariant::CdsVariant { loc_edit, .. }) = &self.c {
            let loc = loc_edit.loc.inner();

            if !self.is_within_cds_bounds() {
                return false;
            }

            let start_offset = loc.start.offset.unwrap_or(0);
            let end_offset = loc.end.offset.unwrap_or(0);

            start_offset == 0 && end_offset == 0
        } else {
            false
        }
    }
}

/// Extension trait for checks on `HgvsVariant`s that consequence prediction needs.
pub(crate) trait HgvsVariantExt {
    /// Check if the alternative transcript of the n. variant depends on an unknown splice outcome.
    ///
    /// This is the case if the variant has an intronic offset at either end. It is also the case
    /// if the variant starts and ends in different exons of `alignment`: then it covers a whole
    /// intron with both splice sites, although its positions have no intronic offset. Insertions
    /// have no reference bases, so they never cover an intron.
    fn alt_depends_on_splicing(&self, alignment: &GenomeAlignment) -> bool;

    /// VEP's `stop_gained`, `stop_lost`, `stop_retained_variant` and `synonymous_variant` for an
    /// indel in the CDS, as in VEP 116.3.
    ///
    /// VEP translates only the codons that the variant changes
    /// (`TranscriptVariationAllele::codon`): the reference codons, and the same stretch of the
    /// altered CDS and 3' UTR. A partial codon at the end of this stretch gives an `X`. A stop
    /// that a new frame reaches behind these codons thus gives no `stop_gained`. After an `X`,
    /// VEP reads the bases at the old stop codon instead (`_ins_del_stop_altered`). For
    /// `stop_retained_variant`, VEP then checks the overlap with the stop codon at `unshifted`,
    /// the input position of the variant before the 3' shift (`_overlaps_stop_codon_cil`). There,
    /// an insertion overlaps the stop codon only if its bases fit in front of the stop codon end.
    ///
    /// Returns `None` for an edit other than an insertion, duplication, deletion or delins, for
    /// an edit that keeps the length, and for an edit with an intronic offset.
    fn vep_indel_stop_terms(
        &self,
        unshifted: &HgvsVariant,
        ref_data: &RefTranscriptData,
        tx: &Transcript,
    ) -> Option<Consequences>;
}

impl HgvsVariantExt for HgvsVariant {
    fn alt_depends_on_splicing(&self, alignment: &GenomeAlignment) -> bool {
        if self.spans_intron() {
            return true;
        }
        let HgvsVariant::TxVariant { loc_edit, .. } = self else {
            return false;
        };
        if loc_edit.edit.inner().is_ins() {
            return false;
        }
        let loc = loc_edit.loc.inner();

        // An intron follows the end of each exon except the last one.
        let exon_ends = alignment.exons.iter().filter_map(|exon| exon.alt_cds_end_i);
        let last_exon_end = exon_ends.clone().max();
        exon_ends
            .filter(|end| Some(*end) != last_exon_end)
            .any(|end| (loc.start.base..loc.end.base).contains(&end))
    }

    fn vep_indel_stop_terms(
        &self,
        unshifted: &HgvsVariant,
        ref_data: &RefTranscriptData,
        tx: &Transcript,
    ) -> Option<Consequences> {
        let codon = |cds_pos: i32| (cds_pos + 2) / 3;

        // The CDS with its stop codon, followed by the 3' UTR. The stored sequence can end with
        // `A` bases that are not part of the transcript.
        let tx_seq = &ref_data.transcript_sequence;
        let tx_end = usize::try_from(tx.tx_len()).ok()?.min(tx_seq.len());
        let cds_and_utr = tx_seq.get(usize::try_from(ref_data.cds_start - 1).ok()?..tx_end)?;
        let cds_len = ref_data.cds_stop - ref_data.cds_start + 1;
        let cds_pos = |pos: &CdsPos| match pos.cds_from {
            CdsFrom::Start => pos.base,
            CdsFrom::End => cds_len + pos.base,
        };
        let index = |cds_pos: i32| usize::try_from(cds_pos).ok();
        // The first and last base that a variant changes, and the bases that replace them. As in
        // VEP, an insertion lies between `end` and `start = end + 1`.
        let changed = |variant: &HgvsVariant| -> Option<(i32, i32, String)> {
            let HgvsVariant::CdsVariant { loc_edit, .. } = variant else {
                return None;
            };
            let loc = loc_edit.loc.inner();
            if loc.start.offset.unwrap_or(0) != 0 || loc.end.offset.unwrap_or(0) != 0 {
                return None;
            }
            let (start, end) = (cds_pos(&loc.start), cds_pos(&loc.end));
            Some(match loc_edit.edit.inner() {
                NaEdit::Ins { alternative } => (start + 1, start, alternative.clone()),
                NaEdit::Dup { .. } => (
                    end + 1,
                    end,
                    cds_and_utr.get(index(start - 1)?..index(end)?)?.to_string(),
                ),
                NaEdit::DelRef { .. } | NaEdit::DelNum { .. } => (start, end, String::new()),
                NaEdit::RefAlt { alternative, .. } => (start, end, alternative.clone()),
                _ => return None,
            })
        };
        let (mut start, mut end, mut alt) = changed(self)?;
        // VEP shifts an insertion or deletion 3' on the genome, also across the CDS end, where
        // hgvs-rs stops. On the genome, an intron follows the end of an exon.
        let base = |cds_pos: i32| cds_and_utr.as_bytes().get(index(cds_pos - 1)?).copied();
        let exon_ends = tx
            .genome_alignments
            .first()?
            .exons
            .iter()
            .filter_map(|exon| Some(exon.alt_cds_end_i? - ref_data.cds_start + 1))
            .collect::<Vec<_>>();
        if start > end {
            while !alt.is_empty()
                && !exon_ends.contains(&end)
                && base(start) == alt.as_bytes().first().copied()
            {
                let first = alt.remove(0);
                alt.push(first);
                (start, end) = (start + 1, end + 1);
            }
        } else if alt.is_empty() {
            while !exon_ends.contains(&end) && base(start).is_some() && base(start) == base(end + 1)
            {
                (start, end) = (start + 1, end + 1);
            }
        }
        let alt = alt.as_str();
        let len_change = alt.len() as i32 - (end - start + 1);
        if len_change == 0 {
            return None;
        }
        // VEP checks the stop codon only for an indel that overlaps the CDS after the 3' shift.
        if start.min(end) > cds_len {
            return Some(Consequences::empty());
        }
        // The altered CDS and 3' UTR. VEP builds none for an indel beyond the transcript end.
        let before = cds_and_utr.get(..index(start - 1)?)?;
        let alt_seq = cds_and_utr
            .get(index(end)?..)
            .map(|behind| [before, alt, behind].concat());

        // The changed codons `first..=last`. An insertion between two codons changes none.
        let (first, last) = (codon(start), codon(end));
        let from = index(3 * (first - 1))?;
        let translate = |codons: &str, exceptions: &[(usize, char)]| {
            translate_cds_with_exceptions(codons, true, "*", ref_data.translation_table, exceptions)
                .ok()
        };
        // VEP applies the translation exceptions, e.g., selenocysteine, to the reference codons.
        let ref_exceptions = ref_data
            .translation_exceptions
            .iter()
            .filter_map(|exception| {
                let codon = usize::try_from(exception.position).ok()?;
                Some((codon.checked_sub(index(first)?)?, exception.amino_acid))
            })
            .collect::<Vec<_>>();
        // The translation of the complete codons, and an `X` for a partial codon behind them
        // unless they translate to a single stop.
        let peptide = |codons: &str, exceptions: &[(usize, char)]| -> Option<String> {
            let complete = codons.len() / 3 * 3;
            let mut peptide = translate(codons.get(..complete)?, exceptions)?;
            if complete < codons.len() && peptide != "*" {
                peptide.push('X');
            }
            Some(peptide)
        };
        let in_cds = (1..=cds_len).contains(&start) && (1..=cds_len).contains(&end);
        let peptides = if in_cds {
            let alt_len = index(3 * (last - first + 1) + len_change)?;
            let ref_codons = cds_and_utr.get(from..index(3 * last)?)?;
            let alt_seq = alt_seq.as_ref()?;
            let alt_codons = alt_seq.get(from..(from + alt_len).min(alt_seq.len()))?;
            Some((
                peptide(ref_codons, &ref_exceptions)?,
                peptide(alt_codons, &[])?,
            ))
        } else {
            None
        };
        // VEP's `partial_codon`: the changed codons start at an incomplete last codon.
        let partial_codon = in_cds && (1..3).contains(&(cds_len - 3 * (first - 1)));

        let incomplete_3p = tx.is_incomplete_3p();
        let overlaps_stop = !incomplete_3p && end >= cds_len - 2 && start <= cds_len;
        // VEP's `_overlaps_stop_codon_cil` works on the genome, at the input position. An
        // insertion overlaps the stop codon only if its bases fit between its position and the
        // stop codon end. An input position in an intron overlaps no stop codon.
        let strand = Strand::try_from(tx.genome_alignments.first()?.strand).ok()?;
        let overlaps_stop_cil = !incomplete_3p
            && changed(unshifted).is_some_and(|(start, end, alt)| {
                let ins_len = alt.len() as i32;
                if start <= end {
                    end >= cds_len - 2 && start <= cds_len
                } else if strand == Strand::Minus {
                    (cds_len - 2 - ins_len..=cds_len - 1).contains(&end)
                } else {
                    (cds_len - 2..=cds_len - 1 - ins_len).contains(&end)
                }
            });
        // VEP's `_ins_del_stop_altered`: the altered CDS and 3' UTR read no stop codon at the
        // position of the old one. Without an altered sequence, VEP counts the stop codon as
        // kept.
        let stop_altered = |overlaps: bool| {
            let stop_codon = index(cds_len - 3)?..index(cds_len)?;
            Some(
                overlaps
                    && len_change != 0
                    && (1..=cds_len).contains(&start)
                    && alt_seq.as_ref().is_some_and(|alt_seq| {
                        alt_seq
                            .get(stop_codon)
                            .and_then(|codon| translate(codon, &[]))
                            .as_deref()
                            != Some("*")
                    }),
            )
        };
        // VEP's `ref_eq_alt_sequence`: the protein with the changed amino acids replaced
        // continues with a stop behind its last amino acid, or the stop stays in place.
        let protein = ref_data
            .aa_sequence
            .strip_suffix('*')
            .unwrap_or(&ref_data.aa_sequence);
        let ref_eq_alt_sequence = |ref_pep: &str, alt_pep: &str| -> Option<bool> {
            let (from, to) = (index(first - 1)?, index(last)?.min(protein.len()));
            if from >= protein.len() {
                return Some(alt_pep.starts_with('*'));
            }
            let mutated = [protein.get(..from)?, alt_pep, protein.get(to..)?].concat();
            Some(
                (mutated.starts_with(protein)
                    && mutated
                        .get(protein.len()..)
                        .is_some_and(|rest| rest.starts_with('*')))
                    || (ref_pep.contains('*') && ref_pep.find('*') == alt_pep.find('*')),
            )
        };

        let stop_lost = !partial_codon
            && match &peptides {
                Some((ref_pep, alt_pep)) if !alt_pep.contains('X') => {
                    ref_pep.contains('*') && !alt_pep.contains('*')
                }
                _ => stop_altered(overlaps_stop)?,
            };
        let stop_retained = !partial_codon
            && !stop_lost
            && match &peptides {
                Some((ref_pep, alt_pep)) if !alt_pep.is_empty() && !alt_pep.contains('X') => {
                    ref_eq_alt_sequence(ref_pep, alt_pep)?
                }
                _ => len_change != 0 && overlaps_stop_cil && !stop_altered(overlaps_stop_cil)?,
            };
        let stop_gained = !stop_lost
            && !stop_retained
            && peptides
                .as_ref()
                .is_some_and(|(ref_pep, alt_pep)| alt_pep.contains('*') && !ref_pep.contains('*'));
        // VEP's `synonymous_variant`: the changed codons translate to the same peptide without
        // an `X`, and VEP reports no `stop_retained_variant` or `start_retained_variant`. The
        // latter holds if the changed codons keep the peptide and include the start codon.
        let synonymous = !stop_retained
            && first > 1
            && peptides.as_ref().is_some_and(|(ref_pep, alt_pep)| {
                !ref_pep.is_empty() && ref_pep == alt_pep && !alt_pep.contains('X')
            });

        Some(if stop_retained {
            Consequence::StopRetainedVariant.into()
        } else if stop_lost {
            Consequence::StopLost.into()
        } else if stop_gained {
            Consequence::StopGained.into()
        } else if synonymous {
            Consequence::SynonymousVariant.into()
        } else {
            Consequences::empty()
        })
    }
}

/// Applies non-overlapping n. variants to the transcript sequence `tx_seq`.
///
/// Returns `None` if a variant is not an n. variant, has an intronic offset, uses an
/// unsupported edit (inversion), or lies outside `tx_seq`.
pub(crate) fn apply_n_edits(tx_seq: &str, vars_n: &[&HgvsVariant]) -> Option<String> {
    struct NEdit {
        replace_start: usize,
        replace_end: usize,
        alt: String,
    }
    let mut n_edits = Vec::with_capacity(vars_n.len());

    for var_n in vars_n {
        let HgvsVariant::TxVariant { loc_edit, .. } = var_n else {
            return None;
        };
        if var_n.spans_intron() {
            return None;
        }
        let loc = loc_edit.loc.inner();
        let edit = loc_edit.edit.inner();

        let n_loc_start = loc.start.base;
        let n_loc_end = loc.end.base;

        if n_loc_start > n_loc_end {
            tracing::warn!(
                "Invalid transcript coordinates ({} > {}) after HGVS projection on transcript {}.",
                n_loc_start,
                n_loc_end,
                var_n.accession().value
            );
            return None;
        }

        let (replace_start, replace_end, alt) = match edit {
            NaEdit::RefAlt { alternative, .. } | NaEdit::NumAlt { alternative, .. } => (
                (n_loc_start - 1) as usize,
                n_loc_end as usize,
                alternative.clone(),
            ),
            NaEdit::DelRef { .. } | NaEdit::DelNum { .. } => (
                (n_loc_start - 1) as usize,
                n_loc_end as usize,
                "".to_string(),
            ),
            NaEdit::Dup { reference } => (
                (n_loc_start - 1) as usize,
                n_loc_end as usize,
                format!("{}{}", reference, reference),
            ),
            NaEdit::Ins { alternative } => (
                n_loc_start as usize,
                n_loc_start as usize,
                alternative.clone(),
            ),
            _ => {
                tracing::warn!(
                    "Unsupported NaEdit type {:?} for transcript sequence assembly. Skipping.",
                    edit
                );
                return None;
            }
        };

        n_edits.push(NEdit {
            replace_start,
            replace_end,
            alt,
        });
    }

    // Apply the edits from right to left, so that an edit does not shift the indices of the
    // edits still to apply. An insertion at index `k` lies left of an edit that replaces bases
    // from `k` on, so of two edits with the same start, the one with the larger end goes first.
    n_edits.sort_by_key(|edit| std::cmp::Reverse((edit.replace_start, edit.replace_end)));

    let mut alt_seq = tx_seq.to_string();
    for edit in &n_edits {
        if edit.replace_start > alt_seq.len() || edit.replace_end > alt_seq.len() {
            tracing::warn!(
                "Edit range out of bounds: {}..{} exceeds sequence length {}. Cannot assemble variant.",
                edit.replace_start,
                edit.replace_end,
                alt_seq.len()
            );
            return None;
        }

        if edit.replace_start > edit.replace_end {
            tracing::warn!(
                "Invalid edit range: start {} > end {}. Cannot assemble variant.",
                edit.replace_start,
                edit.replace_end
            );
            return None;
        }

        alt_seq.replace_range(edit.replace_start..edit.replace_end, &edit.alt);
    }

    Some(alt_seq)
}

#[derive(Debug)]
struct TranscriptLocationContext {
    rank: Rank,
    distance: Option<i32>,
    is_exonic: bool,
    is_intronic: bool,
    is_upstream: bool,
    is_downstream: bool,
}

#[derive(Debug)]
struct ConsequenceContext {
    cds_consequences: Consequences,
    protein_consequences: Consequences,
    cdna_pos: Option<Pos>,
    cds_pos: Option<Pos>,
    protein_pos: Option<Pos>,
}

/// Consequences of one placement of a variant on a transcript, with the location and the HGVS
/// projection they come from.
#[derive(Debug)]
struct PlacementConsequences {
    consequences: Consequences,
    location: TranscriptLocationContext,
    /// `None` up- or downstream of the transcript.
    projection: Option<(HgvsProjectionContext, ConsequenceContext)>,
}

/// One placement of an indel, in 0-based genome coordinates. In a repeat, an indel has several
/// placements with the same alternate sequence.
#[derive(Debug)]
enum Placement {
    /// `bases` deleted from `start` on.
    Del { start: i32, bases: String },
    /// `bases` inserted in front of `pos`.
    Ins { pos: i32, bases: String },
}

impl Placement {
    /// The genome variant of this placement on the contig `accession`.
    fn to_var_g(&self, accession: &Accession) -> HgvsVariant {
        let (start, end, edit) = match self {
            Placement::Del { start, bases } => (
                start + 1,
                start + bases.len() as i32,
                NaEdit::DelRef {
                    reference: bases.clone(),
                },
            ),
            Placement::Ins { pos, bases } => (
                *pos,
                pos + 1,
                NaEdit::Ins {
                    alternative: bases.clone(),
                },
            ),
        };
        HgvsVariant::GenomeVariant {
            accession: accession.clone(),
            gene_symbol: None,
            loc_edit: GenomeLocEdit {
                loc: Mu::Certain(GenomeInterval {
                    start: Some(start),
                    end: Some(end),
                }),
                edit: Mu::Certain(edit),
            },
        }
    }

    /// Whether this placement lies outside the transcript, i.e., changes none of its bases and
    /// no intron either.
    fn is_outside(&self, alignment: &GenomeAlignment) -> bool {
        let (Some(first), Some(last)) = (alignment.exons.first(), alignment.exons.last()) else {
            return true;
        };
        let (tx_start, tx_end) = (first.alt_start_i, last.alt_end_i);
        match self {
            Placement::Del { start, bases } => {
                !overlaps(*start, start + bases.len() as i32, tx_start, tx_end)
            }
            Placement::Ins { pos, .. } => *pos <= tx_start || *pos >= tx_end,
        }
    }

    /// How much of the transcript this placement changes: 0 nothing (intron only),
    /// 1 exon bases outside the CDS, 2 CDS bases, 3 an essential splice site.
    ///
    /// An insertion changes a splice site only between its two bases. An insertion at an exon
    /// edge adds its bases to the exon.
    fn rank(&self, alignment: &GenomeAlignment) -> u8 {
        let exons = &alignment.exons;
        // The two bases at each end of each intron.
        let mut sites = exons.windows(2).flat_map(|pair| {
            let (intron_start, intron_end) = (pair[0].alt_end_i, pair[1].alt_start_i);
            [
                (intron_start, intron_start + 2),
                (intron_end - 2, intron_end),
            ]
        });
        let cds_start = alignment.cds_start.unwrap_or(-1);
        let cds_end = alignment.cds_end.unwrap_or(-1);
        let in_cds = |pos: i32| cds_start <= pos && pos < cds_end;

        match self {
            Placement::Del { start, bases } => {
                let end = start + bases.len() as i32;
                if sites.any(|(site_start, site_end)| overlaps(*start, end, site_start, site_end)) {
                    3
                } else if exons.iter().any(|exon| {
                    overlaps(
                        *start,
                        end,
                        exon.alt_start_i.max(cds_start),
                        exon.alt_end_i.min(cds_end),
                    )
                }) {
                    2
                } else if exons
                    .iter()
                    .any(|exon| overlaps(*start, end, exon.alt_start_i, exon.alt_end_i))
                {
                    1
                } else {
                    0
                }
            }
            Placement::Ins { pos, .. } => {
                if sites.any(|(site_start, _)| *pos == site_start + 1) {
                    return 3;
                }
                // The exon bases in front of and behind the inserted bases, if these join an exon.
                let flanks = exons.iter().enumerate().find_map(|(i, exon)| {
                    if exon.alt_start_i < *pos && *pos < exon.alt_end_i {
                        Some((pos - 1, *pos))
                    } else if *pos == exon.alt_end_i {
                        exons.get(i + 1).map(|next| (pos - 1, next.alt_start_i))
                    } else if *pos == exon.alt_start_i {
                        i.checked_sub(1)
                            .and_then(|j| exons.get(j))
                            .map(|prev| (prev.alt_end_i - 1, *pos))
                    } else {
                        None
                    }
                });
                match flanks {
                    Some((before, behind)) if in_cds(before) && in_cds(behind) => 2,
                    Some(_) => 1,
                    None => 0,
                }
            }
        }
    }
}

impl ConsequencePredictor {
    pub fn new(provider: Arc<MehariProvider>, config: Config) -> Self {
        let mut custom_columns: Vec<_> = config
            .plugins
            .iter()
            .flat_map(|plugin| plugin.columns())
            .collect();
        if config.transcript_settings.flag_pick {
            custom_columns.push(ANN_PICK.to_string());
        }

        tracing::info!("Building transcript interval trees ...");

        let reference_available = provider.reference_available();

        let mapper_config = assembly::Config {
            // TODO: add ability to construct assembly mapper/config with custom mappings (for contig <-> accession lookups) in hgvs-rs to avoid lock-in to bioutils assemblies.
            assembly: provider.assembly(),
            replace_reference: reference_available,
            strict_bounds: false,
            renormalize_g: reference_available && config.renormalize_g,
            genome_seq_available: reference_available,
            normalize: config.normalize,
            // An insertion between an exon and an intron keeps the splice site, so its bases
            // extend the exon, also in the protein. hgvs-rs reads a dup of the intron bases
            // next to an exon as such an insertion.
            ins_at_boundary_is_intronic: false,
            ..Default::default()
        };
        let mapper = assembly::Mapper::new(mapper_config, provider.clone());
        tracing::info!("... done building transcript interval trees");

        ConsequencePredictor {
            provider,
            mapper,
            config,
            custom_columns,
            restrict_dropped_genes: AtomicUsize::new(0),
            restrict_dropped_variants: AtomicUsize::new(0),
        }
    }

    /// The custom ANN columns, in output order.
    pub fn custom_columns(&self) -> &[String] {
        &self.custom_columns
    }

    /// The context of a transcript annotation; empty if no plugin reads it.
    fn tx_context(
        &self,
        tx: Option<&Transcript>,
        var_c: Option<HgvsVariant>,
        var_n: Option<HgvsVariant>,
    ) -> TxContext {
        if self.config.plugins.is_empty() {
            TxContext::empty()
        } else {
            TxContext::new(self.provider.clone(), tx, var_c, var_n)
        }
    }

    /// Run the plugins on the reported annotations of one prediction.
    fn run_plugins(
        &self,
        variants: &[VcfVariant],
        mut txs: Vec<TxAnnotation>,
    ) -> Result<Vec<AnnField>, SeqvarsError> {
        for plugin in &self.config.plugins {
            plugin.annotate(variants, &mut txs)?;
        }
        Ok(txs.into_iter().map(|tx| tx.ann).collect())
    }

    /// Predict the consequences of a variant.
    ///
    /// Note that the predictions will be affected by whether transcript picking has been
    /// enabled in the data provider and the configuration of the predictor, in particular
    /// `Config::report_all_transcripts`.
    ///
    /// # Args
    ///
    /// * `var`: The variant to predict consequences for.
    ///
    /// # Returns
    ///
    /// A list of `AnnField` records, one for each transcript affected by the variant
    /// sorted lexicographically by transcript accession.
    ///
    /// If the accessio is not valid, then `None` will be returned.
    ///
    /// # Errors
    ///
    /// If there was any error during the prediction.
    pub fn predict(&self, var: &VcfVariant) -> Result<Option<Vec<AnnField>>, SeqvarsError> {
        // Normalize variant by stripping common prefix and suffix.
        let mut norm_var = self.normalize_variant(var);

        // TODO check for VCF specification version.
        // According to VCF specification (>=4.1), an alternative of "N" means REF=ALT
        // Prior to 4.1, it indicated a deletion.
        if norm_var.alternative == "N" {
            norm_var.alternative = norm_var.reference.clone();
        }

        // Obtain accession from chromosome name.
        let chrom_acc = self
            .provider
            .contig_manager
            .get_accession(&norm_var.chromosome);
        let chrom_acc = if let Some(chrom_acc) = chrom_acc {
            chrom_acc
        } else {
            tracing::warn!(
                "Could not determine chromosome accession for {:?}; giving up on annotation",
                &norm_var
            );
            return Ok(None);
        };

        // We follow hgvs conventions and therefore normalize input variants
        let var_g = Self::get_var_g(&norm_var, chrom_acc);
        let (var_g_fwd, var_g_rev) = if self.mapper.config.renormalize_g {
            // On circular contigs, normalize one lap further so that shifting 5' can cross the
            // origin, then shift back.  Positions may then lie just beyond either contig end.
            let lap = self.circular_length(chrom_acc).unwrap_or(0);
            let var_g_lap = Self::shift_var_g(&var_g, lap);
            let right = self
                .mapper
                .variant_mapper()
                .right_normalizer()?
                .normalize(&var_g_lap)?;
            let left = self
                .mapper
                .variant_mapper()
                .left_normalizer()?
                .normalize(&var_g_lap)?;
            (
                Self::shift_var_g(&right, -lap),
                Self::shift_var_g(&left, -lap),
            )
        } else {
            (var_g.clone(), var_g.clone())
        };
        // The consequence terms come from one of these placements, see `terms_var_g`.
        let placements = if self.mapper.config.renormalize_g {
            self.equivalent_placements(&var_g_rev, &var_g_fwd)
        } else {
            None
        };

        // Get all affected transcripts.
        let (var_start_fwd, var_end_fwd) = Self::get_var_start_end(&var_g_fwd);
        let (var_start_rev, var_end_rev) = Self::get_var_start_end(&var_g_rev);

        let qry_start = var_start_fwd.min(var_start_rev) - PADDING;
        let qry_end = var_end_fwd.max(var_end_rev) + PADDING;

        let mut txs = self
            .provider
            .get_tx_for_region(chrom_acc, ALT_ALN_METHOD, qry_start, qry_end)
            .map_err(|e| SeqvarsError::Provider(e.to_string()))?;
        txs.sort_by(|a, b| a.tx_ac.cmp(&b.tx_ac));

        // Handle case of no overlapping transcripts -> intergenic.
        if txs.is_empty() {
            if !self.config.keep_intergenic {
                // Plugins may want predictions without annotations, too.
                return Ok(Some(
                    self.run_plugins(std::slice::from_ref(var), Vec::new())?,
                ));
            }
            let hgvs_g = Some(self.format_hgvs_g(&var_g));

            let ann = AnnField {
                allele: Allele::Alt {
                    alternative: var.alternative.clone(),
                },
                gene_id: "".to_string(),
                consequences: vec![Consequence::IntergenicVariant],
                putative_impact:
                    crate::annotate::seqvars::consequence::terms::PutativeImpact::Modifier,
                feature_type: FeatureType::Custom {
                    value: "Intergenic".to_string(),
                },
                feature_id: "".to_string(),
                feature_biotype: vec![FeatureBiotype::Noncoding],
                feature_tags: vec![],
                rank: None,
                distance: None,
                strand: 0,
                hgvs_g,
                hgvs_n: None,
                hgvs_c: None,
                hgvs_p: None,
                cdna_pos: None,
                cds_pos: None,
                protein_pos: None,
                gene_symbol: "".to_string(),
                messages: None,
                custom_fields: BTreeMap::new(),
            };
            let intergenic = TxAnnotation {
                context: self.tx_context(None, None, None),
                ann,
            };
            return Ok(Some(
                self.run_plugins(std::slice::from_ref(var), vec![intergenic])?,
            ));
        }

        let anns = self.select_annotations(txs, |tx| {
            let tx_var_g = if tx.alt_strand == -1 {
                &var_g_rev
            } else {
                &var_g_fwd
            };
            self.build_ann_field(var, tx_var_g, &var_g, placements.as_deref(), tx.clone())
        })?;
        Ok(Some(self.run_plugins(std::slice::from_ref(var), anns)?))
    }

    /// Annotate the transcripts near the variant that `Config::transcript_settings` selects.
    fn select_annotations(
        &self,
        txs: Vec<TxForRegionRecord>,
        annotate: impl FnMut(&TxForRegionRecord) -> Result<Option<TxAnnotation>, SeqvarsError>,
    ) -> Result<Vec<TxAnnotation>, SeqvarsError> {
        let cands = txs
            .into_iter()
            .filter_map(|record| match self.provider.get_tx(&record.tx_ac) {
                Some(tx) => Some(Candidate::new(record, tx)),
                None => {
                    tracing::warn!(
                        "Requested transcript accession {}, got None (potentially filtered)",
                        &record.tx_ac
                    );
                    None
                }
            })
            .collect();
        let selection = pick::select(&self.config.transcript_settings, cands, annotate)?;
        if selection.dropped_genes > 0 {
            self.restrict_dropped_genes
                .fetch_add(selection.dropped_genes, AtomicOrdering::Relaxed);
            self.restrict_dropped_variants
                .fetch_add(1, AtomicOrdering::Relaxed);
        }
        Ok(selection.anns)
    }

    /// Log how many gene annotations `--restrict-transcripts` dropped so far.
    pub(crate) fn log_restrict_summary(&self) {
        let genes = self.restrict_dropped_genes.load(AtomicOrdering::Relaxed);
        if genes > 0 {
            tracing::warn!(
                "--restrict-transcripts: no annotation for {} gene(s) on {} variant(s), as none \
                 of their transcripts has the requested tags",
                genes,
                self.restrict_dropped_variants.load(AtomicOrdering::Relaxed)
            );
        }
    }

    fn get_var_start_end(var_g: &HgvsVariant) -> (i32, i32) {
        match &var_g {
            HgvsVariant::GenomeVariant { loc_edit, .. } => {
                let loc = loc_edit.loc.inner();
                let edit = loc_edit.edit.inner();
                let start = loc
                    .start
                    .map(|s| s - 1)
                    .expect("Failed to get start position");
                let end = loc.end.expect("Failed to get end position");
                // In insertion / duplication cases, range end is exclusive.
                // See https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/
                let end = if edit.is_ins() || edit.is_dup() {
                    end - 1
                } else {
                    end
                };
                (start, end)
            }
            _ => unreachable!(),
        }
    }

    /// All placements of a deletion, insertion or duplication, in genome order, from its left-
    /// and right-normalized form.
    ///
    /// Returns `None` for other variants, or if the genome sequence is not available.
    fn equivalent_placements(
        &self,
        var_g_left: &HgvsVariant,
        var_g_right: &HgvsVariant,
    ) -> Option<Vec<Placement>> {
        let (
            HgvsVariant::GenomeVariant {
                accession,
                loc_edit: left,
                ..
            },
            HgvsVariant::GenomeVariant {
                loc_edit: right, ..
            },
        ) = (var_g_left, var_g_right)
        else {
            return None;
        };
        let left_start = left.loc.inner().start?;
        let (right_start, right_end) = (right.loc.inner().start?, right.loc.inner().end?);
        // Genome sequence of the 0-based range `begin..end`.
        let fetch = |begin: i32, end: i32| {
            self.provider
                .get_seq_part(
                    &accession.value,
                    Some(usize::try_from(begin).ok()?),
                    Some(usize::try_from(end).ok()?),
                )
                .ok()
        };

        match (left.edit.inner(), right.edit.inner()) {
            (
                NaEdit::DelRef { .. } | NaEdit::DelNum { .. },
                NaEdit::DelRef { .. } | NaEdit::DelNum { .. },
            ) => {
                let len = usize::try_from(right_end - right_start + 1).ok()?;
                let first = left_start - 1;
                let seq = fetch(first, right_end)?;
                (first..right_start)
                    .map(|start| {
                        let i = usize::try_from(start - first).ok()?;
                        Some(Placement::Del {
                            start,
                            bases: seq.get(i..i + len)?.to_string(),
                        })
                    })
                    .collect()
            }
            (
                left_edit @ (NaEdit::Ins { .. } | NaEdit::Dup { .. }),
                right_edit @ (NaEdit::Ins { .. } | NaEdit::Dup { .. }),
            ) => {
                // An insertion lies behind its start base. The copy of a duplication can go in
                // front of or behind the duplicated bases.
                let first = if matches!(left_edit, NaEdit::Dup { .. }) {
                    left_start - 1
                } else {
                    left_start
                };
                let (last, bases) = match right_edit {
                    NaEdit::Dup { .. } => (right_end, fetch(right_start - 1, right_end)?),
                    NaEdit::Ins { alternative } => (right_start, alternative.clone()),
                    _ => return None,
                };
                // The alternate sequence from `first` to the end of the inserted bases.
                let alt = fetch(first, last)? + &bases;
                (first..=last)
                    .map(|pos| {
                        let i = usize::try_from(pos - first).ok()?;
                        Some(Placement::Ins {
                            pos,
                            bases: alt.get(i..i + bases.len())?.to_string(),
                        })
                    })
                    .collect()
            }
            _ => None,
        }
    }

    /// The genome variant to compute the consequence terms from, if it differs from `var_g`.
    ///
    /// All `placements` give the same alternate sequence. If one of them keeps the essential
    /// splice site or the stop codon, the alternate sequence keeps it. So the placement that
    /// changes the least of the transcript shows what the variant does. Ties go to a placement
    /// that spares the first and last three bases of each exon: then the alternate sequence
    /// keeps them. The intronic splice region windows are not compared. mehari locates an
    /// insertion by the base behind it for a window that counts from its left end (the donor
    /// region from +3, the fifth base, the polypyrimidine tract on the minus strand, the exonic
    /// splice region at an exon start), and by the base in front of it for the donor and
    /// acceptor sites themselves, which test whether the insertion joins the exon. Remaining
    /// ties go to the placement nearest to `var_g`, the 3'-most one.
    ///
    /// The ends of a transcript do not follow from its sequence. So the terms come from a
    /// placement outside the transcript only if `var_g` lies outside, and then from `var_g`.
    fn terms_var_g(
        &self,
        placements: &[Placement],
        alignment: &GenomeAlignment,
        strand: Strand,
        var_g: &HgvsVariant,
    ) -> Option<HgvsVariant> {
        let HgvsVariant::GenomeVariant { accession, .. } = var_g else {
            return None;
        };
        // The placements from the 3' end on; the first one is `var_g`.
        let mut from_3p = if strand == Strand::Minus {
            itertools::Either::Left(placements.iter())
        } else {
            itertools::Either::Right(placements.iter().rev())
        }
        .peekable();
        if from_3p.peek()?.is_outside(alignment) {
            return None;
        }
        let placement = from_3p
            .filter(|placement| !placement.is_outside(alignment))
            .min_by_key(|placement| {
                (
                    placement.rank(alignment),
                    self.in_exonic_splice_region(placement, accession, alignment, strand),
                )
            })?;
        let terms_var_g = placement.to_var_g(accession);
        (terms_var_g != *var_g).then_some(terms_var_g)
    }

    /// Whether `placement` changes the first or last three bases of an exon, i.e., hits an
    /// exonic splice region window.
    fn in_exonic_splice_region(
        &self,
        placement: &Placement,
        accession: &Accession,
        alignment: &GenomeAlignment,
        strand: Strand,
    ) -> bool {
        let var_g = placement.to_var_g(accession);
        let (var_start, var_end) = Self::get_var_start_end(&var_g);
        let (_, consequences) =
            self.determine_transcript_context(alignment, strand, &var_g, var_start, var_end);
        consequences.contains(Consequence::ExonicSpliceRegionVariant)
    }

    fn determine_transcript_context(
        &self,
        alignment: &GenomeAlignment,
        strand: Strand,
        var_g: &HgvsVariant,
        var_start: i32,
        var_end: i32,
    ) -> (TranscriptLocationContext, Consequences) {
        let mut consequences = Consequences::empty();
        let mut rank = Rank::default();
        let mut is_exonic = false;
        let mut is_intronic = false;
        let mut distance: Option<i32> = None;

        // The range of an insertion is the base in front of it. An insertion at the start of
        // an exon adds its bases to the exon, so its range is the first exon base instead.
        let (var_start, var_end) = match var_g {
            HgvsVariant::GenomeVariant { loc_edit, .. }
                if matches!(loc_edit.edit.inner(), NaEdit::Ins { .. })
                    && (alignment.exons.iter().skip(1)).any(|exon| exon.alt_start_i == var_end) =>
            {
                (var_start + 1, var_end + 1)
            }
            _ => (var_start, var_end),
        };

        let var_overlaps =
            |start: i32, end: i32| -> bool { overlaps(var_start, var_end, start, end) };

        let cds_start = alignment.cds_start.unwrap_or(-1);
        let cds_end = alignment.cds_end.unwrap_or(-1);

        let var_overlaps_cds = var_overlaps(cds_start, cds_end);

        // Insertion ranges are end-exclusive, so subtract 1 from the end, where applicable.
        let ins_shift = Self::ins_shift(var_g);

        // Find first exon that overlaps with variant or intron that contains the variant.
        //
        // Note that exons are stored in genome position order.
        let mut prev_end = None;
        let mut min_start = None;
        let mut max_end = None;

        for exon_alignment in &alignment.exons {
            let exon_start = exon_alignment.alt_start_i;
            let exon_end = exon_alignment.alt_end_i;

            let is_utr = !(var_overlaps_cds || exon_start >= cds_start && exon_end <= cds_end);

            // Check the cases where the variant overlaps with the exon or is contained within an
            // intron.
            if var_overlaps(exon_start, exon_end) {
                rank = Rank {
                    ord: exon_alignment.ord + 1,
                    total: alignment.exons.len() as i32,
                };
                is_exonic = true;
                distance = Some(0);
                consequences |= Self::analyze_exonic_variant(
                    ins_shift, strand, var_start, var_end, exon_start, exon_end, &rank, is_utr,
                );
            } else if let Some(intron_start) = prev_end
                && var_start >= intron_start
                && var_end <= exon_end
                && !is_exonic
            {
                // The intron takes the number of the exon before it in transcript order, which
                // is this exon only on the minus strand.
                rank = Rank {
                    ord: if strand == Strand::Minus {
                        exon_alignment.ord + 1
                    } else {
                        exon_alignment.ord
                    },
                    total: alignment.exons.len() as i32 - 1,
                };
                is_intronic = true;

                // We compute the "distance" with "+1", the first base of the
                // intron is "+1", the last one is "-1".
                let dist_start: i32 = var_start + 1 - intron_start;
                let dist_end: i32 = -(exon_start + 1 - var_end);
                let dist_start_end = if dist_start.abs() <= dist_end.abs() {
                    dist_start
                } else {
                    dist_end
                };
                if distance.is_none()
                    || dist_start_end.abs() <= distance.expect("cannot be None").abs()
                {
                    distance = Some(dist_start_end);
                }
            }

            if let Some(intron_start) = prev_end {
                consequences |= Self::analyze_intronic_variant(
                    ins_shift,
                    alignment,
                    strand,
                    var_start,
                    var_end,
                    intron_start,
                    exon_start,
                    is_utr,
                );
            }

            min_start = Some(std::cmp::min(min_start.unwrap_or(exon_start), exon_start));
            max_end = Some(std::cmp::max(max_end.unwrap_or(exon_end), exon_end));
            prev_end = Some(exon_end);
        }

        let min_start = min_start.expect("must have seen exon");
        let max_end = max_end.expect("must have seen exon");
        let is_upstream = var_end <= min_start;
        let is_downstream = var_start >= max_end;

        if !is_exonic && !is_intronic {
            if is_upstream {
                let val = -(min_start + 1 - var_end);
                if val.abs() <= PADDING {
                    consequences |= match strand {
                        Strand::Plus => Consequence::UpstreamGeneVariant,
                        Strand::Minus => Consequence::DownstreamGeneVariant,
                        _ => unreachable!("invalid strand: {}", alignment.strand),
                    };
                }
                if distance.is_none() {
                    distance = Some(val);
                }
            } else if is_downstream {
                let val = var_start + 1 - max_end;
                if val.abs() <= PADDING {
                    consequences |= match strand {
                        Strand::Plus => Consequence::DownstreamGeneVariant,
                        Strand::Minus => Consequence::UpstreamGeneVariant,
                        _ => unreachable!("invalid strand: {}", alignment.strand),
                    };
                }
                if distance.is_none() {
                    distance = Some(val);
                }
            }
        }

        (
            TranscriptLocationContext {
                rank,
                distance,
                is_exonic,
                is_intronic,
                is_upstream,
                is_downstream,
            },
            consequences,
        )
    }

    /// The HGVS projection of `var_g` onto `tx`. `unshifted_g` is `var_g` at the input position,
    /// before the 3' shift.
    fn project_hgvs(
        &self,
        var_g: &HgvsVariant,
        unshifted_g: &HgvsVariant,
        tx: &Transcript,
        transcript_biotype: TranscriptBiotype,
    ) -> Result<HgvsProjectionContext, SeqvarsError> {
        let mut projection = HgvsProjectionContext {
            n: None,
            c: None,
            p: None,
            region: None,
            unshifted_c: None,
        };

        projection.n = self.mapper.g_to_n(var_g, &tx.id).map_or_else(
            |e| match e {
                Error::NonAdjacentExons(_, _, _, _) => {
                    tracing::warn!("{}, {}: NonAdjacentExons, skipping", &tx.id, var_g);
                    Ok(None)
                }
                _ => Err(SeqvarsError::from(e)),
            },
            |v| Ok(Some(v)),
        )?;

        if let Some(var_n) = &projection.n {
            projection.c = match transcript_biotype {
                TranscriptBiotype::Coding => {
                    self.mapper
                        .variant_mapper()
                        .n_to_c(var_n)
                        .map(Some)
                        .map_err(|e| SeqvarsError::HgvsProjection(e.to_string()))?
                }
                TranscriptBiotype::NonCoding => Some(var_n.clone()),
                _ => None,
            };

            if let Some(var_c) = &projection.c
                && transcript_biotype == TranscriptBiotype::Coding
            {
                projection.region = self.variant_location(var_c);

                // if the variant is purely intronic, we can skip safe_project_c_to_p,
                // and simply inject p.?
                //
                // Except for a dup of the one intron base next to an exon (`c.10+1dup`,
                // `c.11-1dup`): it equals an insertion of that base at the exon edge, which
                // keeps the splice site and extends the exon. So hgvs-rs projects it. A longer
                // dup there has a copy in the intron that keeps the site, see `terms_var_g`.
                let is_purely_intronic = match var_c {
                    HgvsVariant::CdsVariant { loc_edit, .. } => {
                        let loc = loc_edit.loc.inner();
                        let start_offset = loc.start.offset.unwrap_or(0);
                        let end_offset = loc.end.offset.unwrap_or(0);
                        let dup_next_to_exon = matches!(loc_edit.edit.inner(), NaEdit::Dup { .. })
                            && start_offset == end_offset
                            && start_offset.abs() == 1;
                        start_offset != 0
                            && end_offset != 0
                            && loc.start.base == loc.end.base
                            && loc.start.cds_from == loc.end.cds_from
                            && !dup_next_to_exon
                    }
                    _ => false,
                };

                if is_purely_intronic {
                    projection.p = Some(HgvsVariant::ProtVariant {
                        accession: var_c.accession().clone(),
                        gene_symbol: var_c.gene_symbol().clone(),
                        loc_edit: ProtLocEdit::Unknown,
                    });
                } else {
                    projection.p = self.safe_project_c_to_p(var_c)?;
                }
            }
        }

        if self.config.vep_consequence_terms && transcript_biotype == TranscriptBiotype::Coding {
            projection.unshifted_c = self
                .mapper
                .variant_mapper()
                .g_to_c(unshifted_g, &tx.id, ALT_ALN_METHOD)
                .inspect_err(|e| tracing::debug!("{}: cannot map {}: {}", &tx.id, unshifted_g, e))
                .ok();
        }

        Ok(projection)
    }

    fn analyze_transcript_consequences(
        &self,
        projection: &HgvsProjectionContext,
        tx: &Transcript,
        tx_record: &TxForRegionRecord,
        transcript_location: &TranscriptLocationContext,
        tx_len: i32,
        transcript_biotype: TranscriptBiotype,
    ) -> Result<ConsequenceContext, SeqvarsError> {
        let mut context = ConsequenceContext {
            cds_consequences: Consequences::empty(),
            protein_consequences: Consequences::empty(),
            cdna_pos: None,
            cds_pos: None,
            protein_pos: None,
        };

        if let Some(var_n) = &projection.n {
            context.cdna_pos = transcript_location.is_exonic.then_some(match var_n {
                HgvsVariant::TxVariant { loc_edit, .. } => Pos {
                    ord: loc_edit.loc.inner().start.base,
                    total: Some(tx_len),
                },
                _ => panic!("Invalid tx position: {:?}", var_n),
            });
        }

        if let Some(var_c) = &projection.c
            && transcript_biotype == TranscriptBiotype::Coding
        {
            // If there's no stop codon, we can't compute the cds_len.
            // We can, however, still analyze any consequences before the (missing) stop codon.
            let cds_len = tx.stop_codon.map(|stop| stop - tx.start_codon.unwrap());
            let cds_pos = match var_c {
                HgvsVariant::CdsVariant { loc_edit, .. } => {
                    let start = &loc_edit.loc.inner().start;
                    // c.-N (5' UTR) and c.*N (3' UTR) are not positions in the CDS.
                    (start.cds_from == CdsFrom::Start && start.base >= 1).then_some(Pos {
                        ord: start.base,
                        total: cds_len.map(|len| len - tx.stop_codon_padding(tx_len)),
                    })
                }
                _ => panic!("Invalid CDS position: {:?}", var_c),
            };
            context.cds_pos = cds_pos.filter(|_| transcript_location.is_exonic);

            let conservative = is_conservative_cds_variant(var_c);
            let incomplete_3p = tx.is_incomplete_3p();
            let available_cds_len = tx.available_cds_len(tx_len);
            // VEP reports CDS terms for an indel that starts in the CDS, also if the indel can go
            // before c.1, e.g., `c.1del` or `c.1_3dup`. With VEP terms, the position of such an
            // indel counts, not its 5' UTR region.
            let region = projection.region.filter(|&region| {
                !(self.config.vep_consequence_terms && region == VariantLocation::FivePrimeUtr)
            });
            context.cds_consequences = Self::analyze_cds_variant(
                var_c,
                region,
                transcript_location.is_exonic,
                conservative,
                incomplete_3p,
                available_cds_len,
            );

            if let Some(var_p) = &projection.p {
                if matches!(
                    var_p,
                    HgvsVariant::ProtVariant {
                        loc_edit: ProtLocEdit::Unknown,
                        ..
                    }
                ) {
                    // protein_pos and protein_consequences remain intentionally empty (or rather None)
                } else {
                    // Like VEP, do not count the stop codon (if the transcript has one).
                    let prot_len = cds_len
                        .expect("cds_len cannot be None if hgvs.p projection has been successful")
                        / 3
                        - i32::from(!incomplete_3p);
                    context.protein_pos = match var_p {
                        HgvsVariant::ProtVariant { loc_edit, .. } => match loc_edit {
                            ProtLocEdit::Ordinary { loc, .. } => Some(Pos {
                                ord: loc.inner().start.number,
                                total: Some(prot_len),
                            }),
                            _ => None,
                        },
                        _ => panic!("Not a protein position: {:?}", var_p),
                    };

                    context.protein_consequences = self.analyze_protein_variant(
                        var_c,
                        var_p,
                        &context.protein_pos,
                        conservative,
                        &tx_record.tx_ac,
                        incomplete_3p,
                    );

                    let indel_terms = Consequence::FrameshiftVariant
                        | Consequence::ConservativeInframeInsertion
                        | Consequence::DisruptiveInframeInsertion
                        | Consequence::ConservativeInframeDeletion
                        | Consequence::DisruptiveInframeDeletion;
                    let indel_in_cds = context.cds_consequences.intersects(indel_terms);
                    let protein_unchanged = matches!(
                        var_p,
                        HgvsVariant::ProtVariant {
                            loc_edit: ProtLocEdit::NoChange | ProtLocEdit::NoChangeUncertain,
                            ..
                        }
                    );
                    // A silent change, e.g., `p.Phe703=` or `p.Ter704=`.
                    let silent = matches!(
                        var_p,
                        HgvsVariant::ProtVariant {
                            loc_edit: ProtLocEdit::Ordinary { edit, .. },
                            ..
                        } if match edit.inner() {
                            ProteinEdit::Ident => true,
                            ProteinEdit::Subst { alternative } => alternative.is_empty(),
                            _ => false,
                        }
                    );
                    if indel_in_cds
                        && (protein_unchanged || silent)
                        && !self.config.vep_consequence_terms
                    {
                        // An indel that keeps the protein (`p.=`) is synonymous. An indel with a
                        // silent change keeps only the term of that change. VEP reports the indel
                        // terms instead.
                        context.cds_consequences.remove(indel_terms);
                        if protein_unchanged {
                            context.protein_consequences |= Consequence::SynonymousVariant;
                        }
                    } else if indel_in_cds
                        && let Ok(ref_data) = hgvs::mapper::altseq::ref_transcript_data_cached(
                            self.provider.clone(),
                            &tx_record.tx_ac,
                            None,
                        )
                    {
                        // An indel that keeps the protein (`p.=`) keeps the stop codon, if the
                        // CDS has one.
                        if ref_data.aa_sequence.ends_with('*') && protein_unchanged {
                            context.protein_consequences |= Consequence::StopRetainedVariant;
                        }
                    }
                }
            }

            // VEP takes the stop terms and `synonymous_variant` of an indel from the changed
            // codons, whatever `p.` says.
            if self.config.vep_consequence_terms
                && let Ok(ref_data) = hgvs::mapper::altseq::ref_transcript_data_cached(
                    self.provider.clone(),
                    &tx_record.tx_ac,
                    None,
                )
                && let Some(stop_terms) = var_c.vep_indel_stop_terms(
                    projection.unshifted_c.as_ref().unwrap_or(var_c),
                    &ref_data,
                    tx,
                )
            {
                context.protein_consequences.remove(
                    Consequence::StopGained
                        | Consequence::StopLost
                        | Consequence::StopRetainedVariant
                        | Consequence::SynonymousVariant,
                );
                context.protein_consequences |= stop_terms;
            }
        }

        Ok(context)
    }

    /// Consequences of `var_g` on `tx`, with the location and the HGVS projection they come
    /// from. `unshifted_g` is the variant at the input position, before the 3' shift. Returns
    /// `None` if the projection onto the transcript fails.
    fn placement_consequences(
        &self,
        var_g: &HgvsVariant,
        unshifted_g: &HgvsVariant,
        tx: &Transcript,
        tx_record: &TxForRegionRecord,
        alignment: &GenomeAlignment,
        strand: Strand,
        transcript_biotype: TranscriptBiotype,
    ) -> Result<Option<PlacementConsequences>, SeqvarsError> {
        let (var_start, var_end) = Self::get_var_start_end(var_g);
        let (location, mut consequences) =
            self.determine_transcript_context(alignment, strand, var_g, var_start, var_end);

        if location.is_exonic {
            if transcript_biotype == TranscriptBiotype::NonCoding {
                consequences |= Consequence::NonCodingTranscriptExonVariant;
            }
        } else if location.is_intronic {
            if transcript_biotype == TranscriptBiotype::NonCoding {
                consequences |= Consequence::NonCodingTranscriptIntronVariant;
            } else {
                consequences |= Consequence::CodingTranscriptIntronVariant;
            }
        }

        let projection = if !location.is_upstream && !location.is_downstream {
            let projection = self.project_hgvs(var_g, unshifted_g, tx, transcript_biotype)?;
            if projection.n.is_none() {
                return Ok(None);
            }

            let consequence_ctx = self.analyze_transcript_consequences(
                &projection,
                tx,
                tx_record,
                &location,
                tx.tx_len(),
                transcript_biotype,
            )?;

            consequences |= consequence_ctx.cds_consequences | consequence_ctx.protein_consequences;

            self.consequences_fix_special_cases(
                &mut consequences,
                // exon_alignment_consequences, // TODO include these as well
                consequence_ctx.cds_consequences,
                consequence_ctx.protein_consequences,
                &projection,
                tx.is_incomplete_5p(),
            );

            Some((projection, consequence_ctx))
        } else {
            None
        };

        Ok(Some(PlacementConsequences {
            consequences,
            location,
            projection,
        }))
    }

    /// Annotation of `var_g` on the transcript of `tx_record`. `unshifted_g` is `var_g` at the
    /// input position, before the 3' shift.
    ///
    /// The HGVS descriptions and positions come from `var_g`, the consequence terms from the
    /// placement that `terms_var_g` picks from `placements`.
    fn build_ann_field(
        &self,
        orig_var: &VcfVariant,
        var_g: &HgvsVariant,
        unshifted_g: &HgvsVariant,
        placements: Option<&[Placement]>,
        tx_record: TxForRegionRecord,
    ) -> Result<Option<TxAnnotation>, SeqvarsError> {
        let tx = match self.provider.get_tx(&tx_record.tx_ac) {
            Some(tx) => {
                if TranscriptBiotype::try_from(tx.biotype).expect("invalid tx biotype")
                    == TranscriptBiotype::Coding
                    && tx.start_codon.is_none()
                {
                    tracing::debug!(
                        "Skipping transcript {} because it is coding but has no known start codon",
                        &tx_record.tx_ac
                    );
                    return Ok(None);
                }
                tx
            }
            None => {
                tracing::warn!(
                    "Requested transcript accession {}, got None (potentially filtered)",
                    &tx_record.tx_ac
                );
                return Ok(None);
            }
        };

        assert_eq!(
            tx.genome_alignments.len(),
            1,
            "At this point, only one genome alignment is expected"
        );

        let alignment = tx.genome_alignments.first().unwrap();
        let strand = Strand::try_from(alignment.strand).expect("invalid strand");
        let transcript_biotype =
            TranscriptBiotype::try_from(tx.biotype).expect("invalid transcript biotype");

        let Some(at_var_g) = self.placement_consequences(
            var_g,
            unshifted_g,
            tx,
            &tx_record,
            alignment,
            strand,
            transcript_biotype,
        )?
        else {
            return Ok(None); // Early exit if g->n projection failed.
        };
        let at_terms_var_g = placements
            .and_then(|placements| self.terms_var_g(placements, alignment, strand, var_g))
            .and_then(|terms_var_g| {
                self.placement_consequences(
                    &terms_var_g,
                    unshifted_g,
                    tx,
                    &tx_record,
                    alignment,
                    strand,
                    transcript_biotype,
                )
                .inspect_err(|e| {
                    tracing::debug!("{}: keeping the terms of {}: {}", &tx.id, var_g, e);
                })
                .ok()
                .flatten()
            });
        let at_terms_var_g = at_terms_var_g.as_ref().unwrap_or(&at_var_g);
        let mut consequences = at_terms_var_g.consequences;
        let terms_projection = at_terms_var_g.projection.as_ref().map(|(p, _)| p);

        let transcript_location = &at_var_g.location;
        let projection = at_var_g.projection.as_ref().map(|(p, _)| p);
        let (rank, cdna_pos, cds_pos, protein_pos) = match &at_var_g.projection {
            Some((_, ctx)) => (
                Some(transcript_location.rank.clone()),
                ctx.cdna_pos.clone(),
                ctx.cds_pos.clone(),
                ctx.protein_pos.clone(),
            ),
            None => (None, None, None, None),
        };

        let hgvs_g = Some(self.format_hgvs_g(var_g));
        let hgvs_n = projection
            .and_then(|p| p.n.as_ref())
            .map(|var| FormattedLoc(var).to_string());
        let hgvs_c = projection
            .and_then(|p| p.c.as_ref())
            .map(|var| FormattedLoc(var).to_string());
        let hgvs_p = projection
            .and_then(|p| p.p.as_ref())
            .map(|var| FormattedLoc(var).to_string());

        let feature_biotype = vec![match transcript_biotype {
            TranscriptBiotype::Coding => FeatureBiotype::Coding,
            TranscriptBiotype::NonCoding => FeatureBiotype::Noncoding,
            _ => unreachable!("invalid biotype: {:?}", transcript_biotype),
        }];
        let feature_tags = tx
            .tags
            .iter()
            .map(|tag| TranscriptTag::try_from(*tag).expect("invalid transcript tag"))
            .filter(|tag| !matches!(tag, TranscriptTag::EnsemblGraft))
            .filter_map(|transcript_tag| {
                let tag = FeatureTag::from(transcript_tag);
                if !matches!(tag, FeatureTag::Other(_)) {
                    Some(tag)
                } else {
                    None
                }
            })
            .collect_vec();

        if consequences.is_empty() {
            tracing::error!(
                "No consequences for {:?} on {} (hgvs_n={}, hgvs_c={}, hgvs_p={}) - adding `gene_variant`; \
                most likely the transcript has multiple stop codons and the variant \
                lies behind the first.",
                orig_var,
                &tx_record.tx_ac,
                hgvs_n.as_deref().unwrap_or("None"),
                hgvs_c.as_deref().unwrap_or("None"),
                hgvs_p.as_deref().unwrap_or("None")
            );
            consequences |= Consequence::GeneVariant;
        }

        if self.config.vep_consequence_terms {
            self.adjust_vep_terms(&mut consequences, terms_projection);
        }

        let consequences = consequences.iter().collect_vec();
        let putative_impact = (*consequences.first().unwrap()).into();

        let strand = match strand {
            Strand::Unknown => 0,
            Strand::Plus => 1,
            Strand::Minus => -1,
        };

        let ann = AnnField {
            allele: Allele::Alt {
                alternative: orig_var.alternative.clone(),
            },
            consequences,
            putative_impact,
            gene_symbol: tx.gene_symbol.clone(),
            gene_id: tx.gene_id.clone(),
            feature_type: FeatureType::SoTerm {
                term: SoFeature::Transcript,
            },
            feature_id: tx.id.clone(),
            feature_biotype,
            feature_tags,
            rank,
            hgvs_g,
            hgvs_n,
            hgvs_c,
            hgvs_p,
            cdna_pos,
            cds_pos,
            protein_pos,
            strand,
            distance: transcript_location.distance,
            messages: None,
            custom_fields: BTreeMap::new(),
        };
        let region = projection.as_ref().and_then(|p| p.region);
        let alt_depends_on_splicing = projection
            .as_ref()
            .and_then(|p| p.n.as_ref())
            .is_some_and(|var_n| var_n.alt_depends_on_splicing(alignment));
        let (var_c, var_n) = projection.map_or((None, None), |p| (p.c, p.n));

        Ok(Some(TxAnnotation {
            context: self
                .tx_context(Some(tx), var_c, var_n)
                .with_alt_inputs(region, alt_depends_on_splicing),
            ann,
        }))
    }

    fn adjust_vep_terms(
        &self,
        consequences: &mut Consequences,
        projection_context: Option<&HgvsProjectionContext>,
    ) {
        use super::terms::Consequence::*;

        // vep reports the umbrella intron variant term.
        if consequences.contains(CodingTranscriptIntronVariant) {
            consequences.remove(CodingTranscriptIntronVariant);
            consequences.insert(IntronVariant);
        }
        if consequences.contains(NonCodingTranscriptIntronVariant) {
            consequences.remove(NonCodingTranscriptIntronVariant);
            consequences.insert(IntronVariant);
            consequences.insert(NonCodingTranscriptVariant);
        }
        if consequences.contains(FivePrimeUtrIntronVariant) {
            consequences.remove(FivePrimeUtrIntronVariant);
            consequences.insert(IntronVariant);
        }
        if consequences.contains(ThreePrimeUtrIntronVariant) {
            consequences.remove(ThreePrimeUtrIntronVariant);
            consequences.insert(IntronVariant);
        }

        if consequences.contains(FivePrimeUtrExonVariant) {
            consequences.remove(FivePrimeUtrExonVariant);
            consequences.insert(FivePrimeUtrVariant);
        }
        if consequences.contains(ThreePrimeUtrExonVariant) {
            consequences.remove(ThreePrimeUtrExonVariant);
            consequences.insert(ThreePrimeUtrVariant);
        }

        if consequences.contains(SelenocysteineGain | SelenocysteineLoss) {
            consequences.remove(SelenocysteineGain | SelenocysteineLoss);
            consequences.insert(MissenseVariant);
        }

        // VEP marks some intronic variants as coding_sequence_variant,
        // which is wrong, but we will mimic here for compatibility
        let within_cds_bounds = projection_context
            .map(|c| c.is_within_cds_bounds())
            .unwrap_or(false);
        if consequences.contains(ExonicSpliceRegionVariant) {
            consequences.remove(ExonicSpliceRegionVariant);
            consequences.insert(SpliceRegionVariant);
            if within_cds_bounds {
                consequences.insert(CodingSequenceVariant);
            }
        }

        if consequences.intersects(FrameshiftElongation | FrameshiftTruncation) {
            consequences.remove(FrameshiftElongation | FrameshiftTruncation);
            consequences.insert(FrameshiftVariant);
        }

        let inframe_specifics = DisruptiveInframeDeletion
            | DisruptiveInframeInsertion
            | ConservativeInframeDeletion
            | ConservativeInframeInsertion;

        if consequences.intersects(inframe_specifics) {
            if consequences.intersects(DisruptiveInframeDeletion | ConservativeInframeDeletion) {
                consequences.insert(InframeDeletion);
            }
            if consequences.intersects(DisruptiveInframeInsertion | ConservativeInframeInsertion) {
                consequences.insert(InframeInsertion);
            }
            consequences.remove(inframe_specifics);
        }

        let is_inframe = consequences.intersects(InframeDeletion | InframeInsertion);
        let essential_splice = SpliceDonorVariant | SpliceAcceptorVariant;

        if is_inframe {
            if consequences.intersects(essential_splice) {
                consequences.remove(InframeDeletion | InframeInsertion);
                if within_cds_bounds {
                    consequences.insert(CodingSequenceVariant);
                }
            } else if consequences.intersects(
                SpliceRegionVariant
                    | SpliceDonorFifthBaseVariant
                    | SpliceDonorRegionVariant
                    | SplicePolypyrimidineTractVariant,
            ) && within_cds_bounds
            {
                consequences.insert(CodingSequenceVariant);
            }
        }

        if consequences.contains(ExonLossVariant)
            && *consequences != Into::<Consequences>::into(ExonLossVariant)
        {
            consequences.remove(ExonLossVariant);
        }

        let suppress_splice_region = SpliceDonorVariant
            | SpliceAcceptorVariant
            | SpliceDonorFifthBaseVariant
            | SpliceDonorRegionVariant;

        if consequences.intersects(suppress_splice_region) {
            consequences.remove(SpliceRegionVariant);
        }

        if consequences.contains(SpliceDonorFifthBaseVariant) {
            consequences.remove(SpliceDonorRegionVariant);
        }

        if consequences.contains(SpliceAcceptorVariant) {
            consequences.remove(SplicePolypyrimidineTractVariant);
        }

        let suppress_intron =
            SpliceDonorVariant | SpliceAcceptorVariant | SpliceDonorFifthBaseVariant;

        if consequences.intersects(suppress_intron) {
            consequences.remove(IntronVariant);
        }
    }

    #[allow(clippy::too_many_arguments, unused_variables)]
    fn consequences_fix_special_cases(
        &self,
        consequences: &mut Consequences,
        consequences_cds: Consequences,
        consequences_protein: Consequences,
        projection: &HgvsProjectionContext,
        incomplete_5p: bool,
    ) {
        // If we have a transcript_ablation, we can remove all other consequences
        if consequences.contains(Consequence::TranscriptAblation) {
            *consequences = Consequence::TranscriptAblation.into();
            return;
        }

        if let Some(var_c) = projection.c.as_ref() {
            // Do not report splice variants in UTRs.
            if self.config.discard_utr_splice_variants {
                let splice_variants = Consequence::ExonicSpliceRegionVariant
                    | Consequence::SpliceDonorVariant
                    | Consequence::SpliceAcceptorVariant
                    | Consequence::SpliceRegionVariant
                    | Consequence::SpliceDonorFifthBaseVariant
                    | Consequence::SpliceDonorRegionVariant
                    | Consequence::SplicePolypyrimidineTractVariant;
                let utr_intron_variants = Consequence::FivePrimeUtrIntronVariant
                    | Consequence::ThreePrimeUtrIntronVariant;
                let utr_exon_variants =
                    Consequence::FivePrimeUtrExonVariant | Consequence::ThreePrimeUtrExonVariant;
                let is_utr = match var_c {
                    HgvsVariant::CdsVariant { loc_edit, .. } => {
                        let loc = loc_edit.loc.inner();
                        loc.start.base < 0
                            && loc.end.base < 0
                            && loc.start.cds_from == CdsFrom::Start
                            && loc.end.cds_from == CdsFrom::Start
                            || loc.start.base > 0
                                && loc.end.base > 0
                                && loc.start.cds_from == CdsFrom::End
                                && loc.end.cds_from == CdsFrom::End
                    }
                    _ => false, // Not a CDS variant, can't be UTR in this context
                };
                if is_utr
                    && consequences.intersects(splice_variants)
                    && consequences.intersects(utr_intron_variants | utr_exon_variants)
                {
                    *consequences &= !splice_variants;
                }
            }
        }

        // vep simply reports the frameshift on the protein level, irrespective of the
        // actual outcome
        // i.e., this depends on if you want to have the mechanism or the outcome described
        if !self.config.vep_consequence_terms {
            // If a frameshift/ins/del was predicted on the CDS level,
            // but any relevant consequence (i.e. not just GeneVariant) was produced on the protein level,
            // then it is likely that the frameshift induced a more specific consequence.
            let check_cds_csqs: Consequences = Consequence::FrameshiftVariant.into();
            let checked = consequences_cds & check_cds_csqs;
            if checked != Consequences::empty()
                // if the protein consequence is not effectively empty, we remove the CDS frameshift consequence
                && !(consequences_protein.eq(&Consequence::GeneVariant) || consequences_protein.is_empty())
                // if the protein consequence also includes a frameshift, then we keep it,
                // and also next to a stop as the first changed residue (`p.Val170Ter`)
                && !consequences_protein.intersects(
                Consequence::FrameshiftElongation
                    | Consequence::FrameshiftTruncation
                    | Consequence::FrameshiftVariant
                    | Consequence::StopGained,
            ) {
                *consequences &= !checked;
            }
        }

        // In some cases, we predict a stop lost based on the cds variant
        // but the protein translation does not confirm this.
        //
        // e.g.:
        // 20:35511609:CAAGCCGCCTCCAGGTAGCAGCCACAGCCAGGAGCACACAGACAGAAGACTGTGTCATGGGTCATGGCCCCTCCGCACACCTACAGGTTTGCCAAAGGAA:C
        if consequences_cds.contains(Consequence::StopLost)
            && !consequences_protein
                .intersects(Consequence::StopLost | Consequence::ProteinAlteringVariant)
            && projection.p.as_ref().is_some_and(|p| {
                !matches!(
                    p,
                    HgvsVariant::ProtVariant {
                        loc_edit: ProtLocEdit::NoProteinUncertain,
                        ..
                    }
                )
            })
        {
            *consequences &= !Consequence::StopLost;
        }

        // Similarly, for the start lost case
        //
        // e.g.:
        // 13:32316456:TA:T
        // (This case just shortens a poly-A from which the start codon starts)
        if consequences_cds.contains(Consequence::StartLost)
            && !consequences_protein.contains(Consequence::StartLost)
            && !consequences_protein.is_empty()
        {
            *consequences &= !Consequence::StartLost;
        }

        // VEP reports `start_retained_variant` for a deletion in the start codon that leaves the
        // CDS intact. With default terms, hgvs-rs puts such a deletion in the 5' UTR.
        if self.config.vep_consequence_terms
            && let (
                Some(HgvsVariant::TxVariant {
                    loc_edit: n_loc_edit,
                    accession,
                    ..
                }),
                Some(HgvsVariant::CdsVariant {
                    loc_edit: c_loc_edit,
                    ..
                }),
            ) = (&projection.n, &projection.c)
        {
            let n_loc = n_loc_edit.loc.inner();
            let c_edit = c_loc_edit.edit.inner();
            let c_loc = c_loc_edit.loc.inner();

            // If edit occurs within the first 3 bases of the CDS,
            let (start, end) = (c_loc.start.base, c_loc.end.base);
            if within_first_codon(c_loc.start.clone(), c_loc.end.clone()) {
                // … then we need to check whether this is a start lost or a start retained.
                // To that end, check the deletion against the transcript sequence (0-based).
                let start_retained = match c_edit {
                    NaEdit::DelRef { .. } => match (
                        self.provider.get_seq_part(&accession.value, None, None),
                        usize::try_from(n_loc.start.base - start),
                        usize::try_from(n_loc.start.base - 1),
                        usize::try_from(n_loc.end.base),
                    ) {
                        (Ok(tx_seq), Ok(cds_start), Ok(del_start), Ok(del_end)) => {
                            deletion_keeps_cds(&tx_seq, cds_start, del_start..del_end)
                        }
                        _ => false,
                    },
                    // TODO: handle other cases
                    _ => false,
                };
                if start_retained {
                    tracing::trace!("Fixing StartLost → StartRetained for {:?}", &projection,);
                    *consequences &= !Consequence::StartLost;
                    *consequences |= Consequence::StartRetainedVariant;
                }

                // VEP puts the copy of a dup behind the dup. If the inserted bases lie inside
                // the start codon and can go before c.1, VEP reports `start_lost` and
                // `start_retained_variant`, but no in-frame insertion.
                let ins_after = match c_edit {
                    NaEdit::Ins { .. } => Some(start),
                    NaEdit::Dup { .. } => Some(end),
                    _ => None,
                };
                if projection.region == Some(VariantLocation::FivePrimeUtr)
                    && ins_after.is_some_and(|base| base <= 2)
                {
                    *consequences |= Consequence::StartLost | Consequence::StartRetainedVariant;
                    *consequences &= !(Consequence::ConservativeInframeInsertion
                        | Consequence::DisruptiveInframeInsertion);
                }
            }
        }

        if let Some(HgvsVariant::CdsVariant { loc_edit, .. }) = projection.c.as_ref() {
            let loc = loc_edit.loc.inner();
            let start_base = loc.start.base;
            let start_cds_from = loc.start.cds_from;
            let end_base = loc.end.base;
            let end_cds_from = loc.end.cds_from;

            let starts_left_of_start = start_cds_from == CdsFrom::Start && start_base < 0;
            let ends_left_of_start = end_cds_from == CdsFrom::Start && end_base < 0;

            if consequences.contains(Consequence::ExonLossVariant)
                && starts_left_of_start
                && ends_left_of_start
            {
                *consequences &= !Consequence::ExonLossVariant;
            }
        }

        if let Some(HgvsVariant::ProtVariant {
            loc_edit: ProtLocEdit::Unknown,
            ..
        }) = projection.p.as_ref()
            && consequences.is_empty()
            && projection.is_within_coding_sequence()
        {
            *consequences |= Consequence::CodingSequenceVariant;
        }

        // A CDS that is incomplete at its 5' end, e.g. `cds_start_NF`, has no start codon, so a
        // change of its first codon is just a change of the CDS. VEP gives no start terms there
        // either.
        if incomplete_5p {
            *consequences &= !(Consequence::StartLost | Consequence::StartRetainedVariant);
            if consequences.is_empty() {
                *consequences |= Consequence::CodingSequenceVariant;
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn analyze_exonic_variant(
        ins_shift: i32,
        strand: Strand,
        var_start: i32,
        var_end: i32,
        exon_start: i32,
        exon_end: i32,
        rank: &Rank,
        _is_utr: bool,
    ) -> Consequences {
        let mut consequences: Consequences = Consequences::empty();

        let var_overlaps =
            |start: i32, end: i32| -> bool { overlaps(var_start, var_end, start, end) };
        // The range of an insertion is the base in front of it. The insertion moves the bases
        // behind it, so a window that counts from the exon start contains it if it contains
        // the base behind it.
        let var_overlaps_behind = |start: i32, end: i32| -> bool {
            overlaps(var_start + ins_shift, var_end + ins_shift, start, end)
        };

        // Check the cases where the variant overlaps with whole exon.
        if var_start <= exon_start && var_end >= exon_end {
            // FIXME: this is not true if the var_c variant is effectively pre start completely
            // we address that in fix_special_cases
            consequences |= Consequence::ExonLossVariant;
            if var_start < exon_start {
                if strand == Strand::Plus && !rank.is_first() {
                    consequences |= Consequence::SpliceAcceptorVariant;
                } else if strand == Strand::Minus && !rank.is_last() {
                    consequences |= Consequence::SpliceDonorVariant;
                }
            }
            if var_end > exon_end {
                if strand == Strand::Plus && !rank.is_last() {
                    consequences |= Consequence::SpliceDonorVariant;
                } else if strand == Strand::Minus && !rank.is_last() {
                    consequences |= Consequence::SpliceAcceptorVariant;
                }
            }
        }

        // Check splice region variants
        if var_overlaps(exon_end - 3, exon_end) {
            if strand == Strand::Plus {
                if !rank.is_last() {
                    consequences |= Consequence::ExonicSpliceRegionVariant;
                }
            } else {
                // alignment.strand == Strand::Minus
                if !rank.is_first() {
                    consequences |= Consequence::ExonicSpliceRegionVariant;
                }
            }
        }
        if var_overlaps_behind(exon_start, exon_start + 3) {
            if strand == Strand::Plus {
                if !rank.is_first() {
                    consequences |= Consequence::ExonicSpliceRegionVariant;
                }
            } else {
                // alignment.strand == Strand::Minus
                if !rank.is_last() {
                    consequences |= Consequence::ExonicSpliceRegionVariant;
                }
            }
        }
        consequences
    }

    #[allow(clippy::too_many_arguments)]
    fn analyze_intronic_variant(
        ins_shift: i32,
        alignment: &GenomeAlignment,
        strand: Strand,
        var_start: i32,
        var_end: i32,
        intron_start: i32,
        intron_end: i32,
        _is_utr: bool,
    ) -> Consequences {
        let mut consequences: Consequences = Consequences::empty();

        let var_overlaps =
            |start: i32, end: i32| -> bool { overlaps(var_start, var_end, start, end) };
        // The range of an insertion is the base in front of it. The insertion moves the bases
        // behind it, so a window that counts from the intron start, e.g. donor +3 to +8 on the
        // plus strand, contains it if it contains the base behind it.
        let var_overlaps_behind = |start: i32, end: i32| -> bool {
            overlaps(var_start + ins_shift, var_end + ins_shift, start, end)
        };

        // Check the cases where the variant overlaps with the splice acceptor/donor site.
        if var_overlaps(intron_start, intron_start + 2 - ins_shift) {
            // Left side, is acceptor/donor depending on transcript's strand.
            match strand {
                Strand::Plus => {
                    consequences |= Consequence::SpliceDonorVariant;
                }
                Strand::Minus => {
                    consequences |= Consequence::SpliceAcceptorVariant;
                }
                _ => unreachable!("invalid strand: {}", alignment.strand),
            }
        }

        // Check the case where the variant overlaps with the splice donor site.
        if var_overlaps(intron_end - 2, intron_end - ins_shift) {
            // Left side, is acceptor/donor depending on transcript's strand.
            match strand {
                Strand::Plus => {
                    consequences |= Consequence::SpliceAcceptorVariant;
                }
                Strand::Minus => {
                    consequences |= Consequence::SpliceDonorVariant;
                }
                _ => unreachable!("invalid strand: {}", alignment.strand),
            }
        }

        // Check the case where the variant overlaps with the splice region (1-3 bases in exon
        // or 3-8 bases in intron).
        // n.b. the 1-3 bases in exon check is already done within `analyze_exonic_variant`.
        // We have to check all cases independently and not with `else`
        // because the variant may be larger.
        if var_overlaps_behind(intron_start + 2, intron_start + 8)
            || var_overlaps(intron_end - 8, intron_end - 2)
        {
            consequences |= Consequence::SpliceRegionVariant;
        }

        // Check the case where the variant overlaps with the polypyrimidine tract.
        // (A sequence variant that falls in the polypyrimidine tract at 3' end of intron between 17 and 3 bases from the end (acceptor -3 to acceptor -17))
        if strand == Strand::Plus && var_overlaps(intron_end - 17, intron_end - 2) {
            consequences |= Consequence::SplicePolypyrimidineTractVariant;
        }
        if strand == Strand::Minus && var_overlaps_behind(intron_start + 2, intron_start + 17) {
            consequences |= Consequence::SplicePolypyrimidineTractVariant;
        }

        // Check conditions for splice_donor_region_variant
        // (A sequence variant that falls in the region between the 3rd and 6th base after splice junction (5' end of intron))
        // Note that this is two bases short of the intronic part of splice_region_variant
        if strand == Strand::Plus && var_overlaps_behind(intron_start + 2, intron_start + 6) {
            consequences |= Consequence::SpliceDonorRegionVariant;
        }
        if strand == Strand::Minus && var_overlaps(intron_end - 6, intron_end - 2) {
            consequences |= Consequence::SpliceDonorRegionVariant;
        }

        // Check conditions for splice_donor_5th_base_variant
        // (A sequence variant that causes a change at the 5th base pair after the start of the intron in the orientation of the transcript.)
        if strand == Strand::Plus && var_overlaps_behind(intron_start + 4, intron_start + 5) {
            consequences |= Consequence::SpliceDonorFifthBaseVariant;
        }
        if strand == Strand::Minus && var_overlaps(intron_end - 5, intron_end - 4) {
            consequences |= Consequence::SpliceDonorFifthBaseVariant;
        }

        consequences
    }

    fn ins_shift(var_g: &HgvsVariant) -> i32 {
        // The range of an insertion is the base in front of it. An insertion changes a
        // splice site only between its two bases, i.e., behind the first one. We can express
        // this with a shift of 1 that drops the second site base.

        match var_g {
            HgvsVariant::GenomeVariant {
                loc_edit: GenomeLocEdit { edit, .. },
                ..
            } => {
                let edit = edit.inner();
                if edit.is_ins() || edit.is_dup() { 1 } else { 0 }
            }
            _ => unreachable!(),
        }
    }

    fn analyze_cds_variant(
        var_c: &HgvsVariant,
        region: Option<VariantLocation>,
        is_exonic: bool,
        conservative: bool,
        incomplete_3p: bool,
        available_cds_len: Option<i32>,
    ) -> Consequences {
        let mut consequences: Consequences = Consequences::empty();

        if let HgvsVariant::CdsVariant { loc_edit, .. } = &var_c {
            // Handle the cases where the variant touches the start or stop codon based on `var_c`
            // coordinates.  The cases where the start/stop codon is touched by the variant
            // directly is handled above based on the `var_p` prediction.
            let loc = loc_edit.loc.inner();
            let edit = loc_edit.edit.inner();
            let start_base = loc.start.base;
            let start_cds_from = loc.start.cds_from;
            let end_base = loc.end.base;
            let end_cds_from = loc.end.cds_from;
            let loc_start_offset = loc.start.offset.unwrap_or(0);
            let loc_end_offset = loc.end.offset.unwrap_or(0);

            // Update is_intronic flag with information from var_c.
            // From hgvs spec:
            // > Base-Offset coordinates use a base position,
            //   which is an index in the specified sequence,
            //   and an optional offset from that base position.
            //   Non-zero offsets refer to non-coding sequence,
            //   such as 5’ UTR, 3’ UTR, or intronic position.
            let is_intronic_or_utr = loc_start_offset != 0 && loc_end_offset != 0;

            // An insertion at an exon edge (`c.10_10+1ins`, `c.11-1_11ins`) keeps the splice
            // site and adds its bases to the exon. In front of the first CDS base (`c.1-1_1ins`),
            // they join the 5' UTR. Behind the last one (`c.N_N+1ins`), they join the 3' UTR.
            let at_exon_edge = matches!(edit, NaEdit::Ins { .. })
                && start_base == end_base
                && start_cds_from == end_cds_from
                && matches!((loc_start_offset, loc_end_offset), (0, 1) | (-1, 0));
            let ins_in_front_of_start = at_exon_edge
                && start_cds_from == CdsFrom::Start
                && start_base == 1
                && loc_start_offset == -1;
            let ins_behind_stop = at_exon_edge
                && end_cds_from == CdsFrom::Start
                && loc_end_offset == 1
                && !incomplete_3p
                && available_cds_len == Some(end_base);

            // The variables below mean "VARIANT_{starts,stops}_{left,right}_OF_{start,stop}_CODON".
            //
            // start codon
            let starts_left_of_start = start_cds_from == CdsFrom::Start && start_base < 0;
            let ends_right_of_start = start_cds_from != CdsFrom::Start || start_base > 0;
            if starts_left_of_start && ends_right_of_start {
                consequences |= Consequence::StartLost;
            }
            // hgvs-rs places an insertion or duplication in a UTR if it leaves the CDS intact,
            // e.g., `c.1dup` in front of the start codon or `c.Ndup` behind the stop codon.
            // Without a stop codon, the CDS goes on behind `c.N`, so hgvs-rs keeps a `c.Ndup` in
            // the CDS.
            let in_5_prime_utr = region == Some(VariantLocation::FivePrimeUtr);
            let in_3_prime_utr = region == Some(VariantLocation::ThreePrimeUtr);
            // stop codon
            let starts_left_of_stop = start_cds_from == CdsFrom::Start;
            let ends_right_of_stop =
                end_cds_from == CdsFrom::End || in_3_prime_utr || ins_behind_stop;
            if starts_left_of_stop
                && ends_right_of_stop
                && !incomplete_3p
                && !in_3_prime_utr
                && !ins_behind_stop
            {
                consequences |= Consequence::StopLost;
            }

            // The last codon starts at the last codon boundary. It is partial if the CDS
            // length is not a multiple of 3.
            if incomplete_3p
                && let Some(cds_len) = available_cds_len
                && start_base <= cds_len
                && end_base >= cds_len - (cds_len - 1) % 3
            {
                consequences |= Consequence::IncompleteTerminalCodonVariant;
            }

            if (start_cds_from == CdsFrom::Start && start_base <= 0)
                && ends_right_of_stop
                && matches!(edit, NaEdit::DelNum { .. } | NaEdit::DelRef { .. })
            {
                consequences |= Consequence::TranscriptAblation;
            }

            // Detect variants affecting the 5'/3' UTRs.
            if (starts_left_of_start && start_base < 0) || in_5_prime_utr || ins_in_front_of_start {
                if is_intronic_or_utr {
                    consequences |= Consequence::FivePrimeUtrIntronVariant;
                } else if is_exonic {
                    consequences |= Consequence::FivePrimeUtrExonVariant;
                }
            }
            if ends_right_of_stop {
                if is_intronic_or_utr {
                    consequences |= Consequence::ThreePrimeUtrIntronVariant;
                } else if is_exonic {
                    consequences |= Consequence::ThreePrimeUtrExonVariant;
                }
            }

            if matches!(edit, NaEdit::DelNum { .. } | NaEdit::DelRef { .. })
                && end_base > start_base
                && end_base > 0
            {
                // A deletion that can go before c.1 leaves the start codon intact.
                if start_base <= 3 && start_cds_from == CdsFrom::Start && !in_5_prime_utr {
                    consequences |= Consequence::StartLost;
                }

                if loc_start_offset < 0 && loc_end_offset > 0 {
                    consequences |= Consequence::ExonLossVariant;
                }
            }

            // Make sure not to report frameshift variants
            // that occur completely within intronic sequence
            // i.e. not within the CDS, as the definition is
            // "A sequence variant which causes a disruption of the translational reading frame,
            // because the number of nucleotides inserted or deleted is not a multiple of three."
            let within_exonic_sequence =
                loc_start_offset == 0 && loc_end_offset == 0 || at_exon_edge;
            let _crosses_boundary = (loc_start_offset != 0) ^ (loc_end_offset != 0);

            if !(ends_right_of_stop
                || starts_left_of_start
                || in_5_prime_utr
                || ins_in_front_of_start
                || is_intronic_or_utr)
            {
                match edit {
                    NaEdit::RefAlt {
                        reference,
                        alternative,
                    } => {
                        if reference.len().abs_diff(alternative.len()) % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else {
                            // Check for inframe insertions/deletions (that are not delins)
                            match (
                                reference.len().cmp(&alternative.len()),
                                alternative.is_empty() ^ reference.is_empty(),
                            ) {
                                (Ordering::Less, true) => {
                                    if conservative {
                                        consequences |= Consequence::ConservativeInframeInsertion;
                                    } else {
                                        consequences |= Consequence::DisruptiveInframeInsertion;
                                    }
                                }
                                (Ordering::Greater, true) => {
                                    if conservative {
                                        consequences |= Consequence::ConservativeInframeDeletion;
                                    } else {
                                        consequences |= Consequence::DisruptiveInframeDeletion;
                                    }
                                }
                                _ => {}
                            }
                        }
                    }
                    NaEdit::DelRef { reference } => {
                        if reference.len() % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else if conservative {
                            consequences |= Consequence::ConservativeInframeDeletion;
                        } else {
                            consequences |= Consequence::DisruptiveInframeDeletion;
                        }
                    }
                    NaEdit::DelNum { count } => {
                        if count % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else if conservative {
                            consequences |= Consequence::ConservativeInframeDeletion;
                        } else {
                            consequences |= Consequence::DisruptiveInframeDeletion;
                        }
                    }
                    NaEdit::Ins { alternative } => {
                        if alternative.len() % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else if conservative {
                            consequences |= Consequence::ConservativeInframeInsertion;
                        } else {
                            consequences |= Consequence::DisruptiveInframeInsertion;
                        }
                    }
                    NaEdit::Dup { reference } => {
                        if reference.len() % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else if conservative {
                            consequences |= Consequence::ConservativeInframeInsertion;
                        } else {
                            consequences |= Consequence::DisruptiveInframeInsertion;
                        }
                    }
                    _ => {}
                }
            }
        } else {
            panic!("Must be CDS variant: {}", var_c)
        };
        consequences
    }

    fn analyze_protein_variant(
        &self,
        var_c: &HgvsVariant,
        var_p: &HgvsVariant,
        protein_pos: &Option<Pos>,
        conservative: bool,
        tx_accession: &str,
        incomplete_3p: bool,
    ) -> Consequences {
        let mut consequences: Consequences = Consequences::empty();

        // TODO move to hgvs-rs library as method of `ProtPos` or similar
        fn is_stop(s: &str) -> bool {
            s == "X" || s == "Ter" || s == "*"
        }

        // hgvs-rs gives `X` (`Xaa`) for an incomplete codon at the transcript end, not for a
        // stop codon.
        fn has_stop(s: &str) -> bool {
            s.contains('*') || s.contains("Ter")
        }

        match var_p {
            HgvsVariant::ProtVariant { loc_edit, .. } => match loc_edit {
                ProtLocEdit::Ordinary { loc, edit } => {
                    let loc = loc.inner();
                    match edit.inner() {
                        ProteinEdit::Fs { .. } => {
                            consequences |= Consequence::FrameshiftVariant;

                            // in the case of frameshifts, we will get the altered protein sequence
                            // in order to compare it with the unaltered one

                            if let Ok(reference_data) =
                                hgvs::mapper::altseq::ref_transcript_data_cached(
                                    self.provider.clone(),
                                    tx_accession,
                                    None,
                                )
                            {
                                let original_sequence_len = reference_data.aa_sequence.len();
                                if let Ok(alt_data) =
                                    self.alt_seq_builder(var_c, &reference_data).build_altseq()
                                    && let Some(alt_data) = alt_data.first()
                                {
                                    let altered_sequence = &alt_data.aa_sequence;

                                    // Compare the lengths up to the first stop of the new frame.
                                    // A new frame without a stop (`fsTer?`) runs to the
                                    // transcript end. It is an elongation only if it reads past
                                    // the reference stop codon. Without a reference stop codon,
                                    // hgvs cuts it to the reference length.
                                    //
                                    // do not use the 'X' fallback here,
                                    // as that is _usually_ only added
                                    // when the number of bases is not divisible by 3.
                                    // We only want to identify cases where a new/later
                                    // stop codon is encountered
                                    // .or_else(|| altered_sequence.find('X'))
                                    if let Some(pos) = altered_sequence.find('*') {
                                        // Count the amino acids before the stops. A CDS with an
                                        // incomplete end has no stop codon.
                                        let original_amino_acids =
                                            reference_data.aa_sequence.trim_end_matches('*').len();
                                        match pos.cmp(&original_amino_acids) {
                                            Ordering::Less => {
                                                consequences |= Consequence::FrameshiftTruncation;
                                            }
                                            Ordering::Equal => {
                                                if !self.config.vep_consequence_terms {
                                                    consequences |= Consequence::MissenseVariant;
                                                    // TODO: discuss stop_retained
                                                    // consequences |= Consequence::StopRetainedVariant;
                                                    consequences &= !Consequence::FrameshiftVariant;
                                                }
                                            }
                                            Ordering::Greater => {
                                                consequences |= Consequence::FrameshiftElongation;
                                            }
                                        }
                                    } else if reference_data.aa_sequence.ends_with('*')
                                        && altered_sequence.len() > original_sequence_len
                                    {
                                        consequences |= Consequence::FrameshiftElongation;
                                    }
                                }
                            }
                        }
                        ProteinEdit::Ext { .. } => {
                            if !incomplete_3p {
                                consequences |= Consequence::StopLost;
                            }
                            consequences |= Consequence::FeatureElongation;
                        }
                        ProteinEdit::Subst { alternative } => {
                            if alternative.is_empty() {
                                consequences |= Consequence::SynonymousVariant;
                            } else if is_stop(alternative) {
                                if loc.start == loc.end && is_stop(&loc.start.aa) {
                                    consequences |= Consequence::StopRetainedVariant;
                                } else {
                                    consequences |= Consequence::StopGained;
                                    // if the substitution happens right before the stop codon
                                    // and if it is a conservative change
                                    // then it is not a stop gained
                                    // cf. 1:43450470:GCCT:G, ENST00000634258.3:c.10294_10296del/p.Leu3432Ter
                                    if let Some(p) = protein_pos
                                        && p.total.is_some_and(|t| p.ord == t)
                                        && conservative
                                    {
                                        consequences &= !Consequence::StopGained;
                                        consequences |= Consequence::ConservativeInframeDeletion;
                                    }
                                }
                            } else {
                                consequences |= Consequence::MissenseVariant;
                                // Missense variants that affect selenocysteine are marked
                                // as rare amino acid variants / selenocysteine gain/loss variants.
                                let alt_has_selenocysteine = alternative.contains('U');
                                let ref_has_selenocysteine =
                                    (loc.start == loc.end) && loc.start.aa == "U";
                                match (ref_has_selenocysteine, alt_has_selenocysteine) {
                                    (true, false) => {
                                        consequences |= Consequence::SelenocysteineLoss;
                                    }
                                    (false, true) => {
                                        consequences |= Consequence::SelenocysteineGain;
                                    }
                                    (true, true) => {
                                        consequences |= Consequence::RareAminoAcidVariant;
                                    }
                                    _ => {}
                                }
                            }
                        }
                        ProteinEdit::DelIns { alternative } => {
                            consequences |= Consequence::ProteinAlteringVariant;
                            if alternative
                                .len()
                                .cmp(&(loc.start.number.abs_diff(loc.end.number) as usize + 1))
                                == Ordering::Equal
                            {
                                // When the delins does not change the CDS length,
                                // it is a missense variant, not an inframe deletion
                                // cf https://github.com/Ensembl/ensembl-vep/issues/1388
                                consequences |= Consequence::MissenseVariant;
                                consequences &= !Consequence::ProteinAlteringVariant;
                            }

                            if (is_stop(&loc.start.aa) || is_stop(&loc.end.aa))
                                && !has_stop(alternative)
                                && !incomplete_3p
                            {
                                consequences |= Consequence::StopLost;
                            }

                            if has_stop(alternative) {
                                consequences |= Consequence::StopGained;
                            }
                        }
                        // Inserted bases that encode a stop codon give an insertion up to the
                        // stop, e.g., `p.Lys168_Ile169insThrValSerTer`, also if their number is
                        // not a multiple of 3, as HGVS prescribes. The frame still shifts then.
                        ProteinEdit::Ins { alternative } if has_stop(alternative) => {
                            consequences |= Consequence::StopGained;
                            consequences |= if inserts_frameshift(var_c) {
                                Consequence::FrameshiftVariant
                            } else if conservative {
                                Consequence::ConservativeInframeInsertion
                            } else {
                                Consequence::DisruptiveInframeInsertion
                            };
                        }
                        ProteinEdit::Ins { .. } | ProteinEdit::Dup => {
                            if conservative {
                                consequences |= Consequence::ConservativeInframeInsertion;
                            } else {
                                consequences |= Consequence::DisruptiveInframeInsertion;
                            }
                        }
                        ProteinEdit::Del => {
                            if conservative {
                                consequences |= Consequence::ConservativeInframeDeletion;
                            } else {
                                consequences |= Consequence::DisruptiveInframeDeletion;
                            }
                        }
                        ProteinEdit::Ident => {
                            if loc.start == loc.end && is_stop(&loc.start.aa) {
                                consequences |= Consequence::StopRetainedVariant;
                            } else {
                                consequences |= Consequence::SynonymousVariant;
                            }
                        }
                    };
                }
                // The protein does not change (`p.=`), e.g., for an insertion that keeps the stop
                // codon.
                ProtLocEdit::NoChange | ProtLocEdit::NoChangeUncertain => {}
                ProtLocEdit::InitiationUncertain => {
                    consequences |= Consequence::StartLost;
                }
                ProtLocEdit::NoProtein | ProtLocEdit::NoProteinUncertain | ProtLocEdit::Unknown => {
                }
            },
            _ => panic!("Must be protein variant: {}", var_p),
        }
        consequences
    }

    /// Length of `chrom_acc` if it is a circular contig.
    fn circular_length(&self, chrom_acc: &str) -> Option<i32> {
        self.provider
            .contig_manager
            .get_circular_length(chrom_acc)
            .and_then(|len| i32::try_from(len).ok())
    }

    /// Shift the positions of the genome variant `var_g` by `offset`.
    fn shift_var_g(var_g: &HgvsVariant, offset: i32) -> HgvsVariant {
        let mut var_g = var_g.clone();
        if let HgvsVariant::GenomeVariant { loc_edit, .. } = &mut var_g {
            let loc = loc_edit.loc.inner_mut();
            loc.start = loc.start.map(|pos| pos + offset);
            loc.end = loc.end.map(|pos| pos + offset);
        }
        var_g
    }

    /// Format `var_g` for the `hgvs_g` field.  On circular contigs, this uses the `m.` prefix
    /// and folds positions into `1..=len`, so a start after the end spans the origin.
    fn format_hgvs_g(&self, var_g: &HgvsVariant) -> String {
        if let HgvsVariant::GenomeVariant {
            accession,
            gene_symbol,
            loc_edit,
        } = var_g
            && let Some(len) = self.circular_length(&accession.value)
        {
            let fold = |pos: Option<i32>| pos.map(|pos| (pos - 1).rem_euclid(len) + 1);
            let loc = loc_edit.loc.inner();
            let var_m = HgvsVariant::MtVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: MtLocEdit {
                    loc: Mu::from(
                        MtInterval {
                            start: fold(loc.start),
                            end: fold(loc.end),
                        },
                        loc_edit.loc.is_certain(),
                    ),
                    edit: loc_edit.edit.clone(),
                },
            };
            FormattedLoc(&var_m).to_string()
        } else {
            FormattedLoc(var_g).to_string()
        }
    }

    fn get_var_g(var: &VcfVariant, chrom_acc: &str) -> HgvsVariant {
        let chrom_acc = chrom_acc.to_string();
        HgvsVariant::GenomeVariant {
            accession: Accession { value: chrom_acc },
            gene_symbol: None,
            loc_edit: if var.reference.is_empty() {
                // insertion
                GenomeLocEdit {
                    loc: Mu::Certain(GenomeInterval {
                        start: Some(var.position - 1),
                        end: Some(var.position),
                    }),
                    edit: Mu::Certain(NaEdit::Ins {
                        alternative: var.alternative.clone(),
                    }),
                }
            } else if var.alternative.is_empty() {
                // deletion
                GenomeLocEdit {
                    loc: Mu::Certain(GenomeInterval {
                        start: Some(var.position),
                        end: Some(var.position + var.reference.len() as i32 - 1),
                    }),
                    edit: Mu::Certain(NaEdit::DelRef {
                        reference: var.reference.clone(),
                    }),
                }
            } else {
                // substitution
                GenomeLocEdit {
                    loc: Mu::Certain(GenomeInterval {
                        start: Some(var.position),
                        end: Some(var.position + var.reference.len() as i32 - 1),
                    }),
                    edit: Mu::Certain(NaEdit::RefAlt {
                        reference: var.reference.clone(),
                        alternative: var.alternative.clone(),
                    }),
                }
            },
        }
    }

    // Normalize variant by stripping common suffixes and prefixes.
    fn normalize_variant(&self, var: &VcfVariant) -> VcfVariant {
        let mut result = var.clone();

        // Strip common suffixes.
        while result.reference.len() > 1 && result.alternative.len() > 1 {
            if result.reference.chars().last().unwrap()
                == result.alternative.chars().last().unwrap()
            {
                result.reference.pop();
                result.alternative.pop();
            } else {
                break;
            }
        }

        // Strip common suffixes.
        while !result.reference.is_empty() && !result.alternative.is_empty() {
            if result.reference.chars().next().unwrap()
                == result.alternative.chars().next().unwrap()
            {
                result.position += 1;
                result.reference.remove(0);
                result.alternative.remove(0);
            } else {
                break;
            }
        }

        result
    }

    /// Predict the combined consequence of multiple (phased) variants.
    /// This is an _experimental_ feature.
    pub fn predict_multiple(
        &self,
        variants: &[VcfVariant],
    ) -> Result<Option<Vec<AnnField>>, SeqvarsError> {
        // Run each input variant through single-variant normalization first.
        let normalized_variants: Vec<VcfVariant> = variants
            .iter()
            .map(|v| {
                let (_, r1, a1) =
                    hgvs::sequences::trim_common_suffixes_slice(&v.reference, &v.alternative);
                let (prefix_trim, r2, a2) = hgvs::sequences::trim_common_prefixes_slice(r1, a1);
                // Detect when a2 == "N" and rewrite ALT to the normalized REF
                let alternative = if a2 == "N" {
                    r2.to_string()
                } else {
                    a2.to_string()
                };
                VcfVariant {
                    chromosome: v.chromosome.clone(),
                    position: v.position + prefix_trim as i32,
                    reference: r2.to_string(),
                    alternative,
                }
            })
            .collect();

        let mut paired_vars: Vec<_> = variants.iter().cloned().zip(normalized_variants).collect();
        paired_vars.sort_by_key(|(_, norm)| norm.position);

        let sorted_originals: Vec<VcfVariant> =
            paired_vars.iter().map(|(orig, _)| orig.clone()).collect();
        let sorted_normalized: Vec<VcfVariant> =
            paired_vars.into_iter().map(|(_, norm)| norm).collect();

        let sorted_vars = Self::validate_and_sort_variant_group(&sorted_normalized)?;
        if sorted_vars.is_empty() {
            return Ok(None);
        }

        let chrom_acc = self
            .provider
            .contig_manager
            .get_accession(&sorted_vars[0].chromosome)
            .ok_or_else(|| SeqvarsError::UnknownChromosomeAccession)?;

        // build a pseudo hgvs.g description
        // Treat insertions (reference.len()==0) as interbase events
        let min_var = sorted_vars.first().unwrap();
        let min_pos = if min_var.reference.is_empty() {
            min_var.position - 1
        } else {
            min_var.position
        };

        let max_var = sorted_vars.last().unwrap();
        let max_pos = if max_var.reference.is_empty() {
            max_var.position
        } else {
            max_var.position + max_var.reference.len() as i32 - 1
        };

        // Reject two or more insertion-only variants at the same normalized position
        let insertion_positions: Vec<i32> = sorted_vars
            .iter()
            .filter(|v| v.reference.is_empty())
            .map(|v| v.position)
            .collect();
        if insertion_positions.len() >= 2 {
            let mut unique_positions = insertion_positions.clone();
            unique_positions.sort_unstable();
            unique_positions.dedup();
            if unique_positions.len() < insertion_positions.len() {
                return Err(SeqvarsError::GroupValidation(
                    GroupValidationError::MultipleInsertionsSamePosition,
                ));
            }
        }

        let mut hgvs_g = None;
        if let Ok(ref_seq_g) = self.provider.get_seq_part(
            chrom_acc,
            Some((min_pos as usize).saturating_sub(1)),
            Some(max_pos as usize),
        ) {
            let mut alt_seq_g = ref_seq_g.clone();

            let mut g_edits: Vec<_> = sorted_vars
                .iter()
                .map(|var| {
                    let start = (var.position - min_pos) as usize;
                    let end = start + var.reference.len();
                    (start, end, var.alternative.clone())
                })
                .collect();

            // Apply the edits from right to left, like `apply_n_edits`.
            g_edits.sort_by_key(|(start, end, _)| std::cmp::Reverse((*start, *end)));
            for (start, end, alt) in g_edits {
                alt_seq_g.replace_range(start..end, &alt);
            }

            let compound_var_g = HgvsVariant::GenomeVariant {
                accession: Accession::new(chrom_acc),
                gene_symbol: None,
                loc_edit: GenomeLocEdit {
                    loc: Mu::Certain(GenomeInterval {
                        start: Some(min_pos),
                        end: Some(max_pos),
                    }),
                    edit: Mu::Certain(NaEdit::RefAlt {
                        reference: ref_seq_g,
                        alternative: alt_seq_g,
                    }),
                },
            };
            hgvs_g = Some(self.format_hgvs_g(&compound_var_g));
        }

        let txs = self.get_transcripts_for_variant_group(&sorted_vars, chrom_acc)?;
        if txs.is_empty() {
            return Ok(None);
        }

        let multi_anns = self.select_annotations(txs, |tx_record| {
            self.predict_multiple_for_transcript(
                &sorted_vars,
                &sorted_originals,
                tx_record,
                chrom_acc,
                hgvs_g.clone(),
            )
        })?;

        Ok(Some(self.run_plugins(variants, multi_anns)?))
    }

    /// Ensures variants don't overlap and are in the correct order.
    fn validate_and_sort_variant_group(
        variants: &[VcfVariant],
    ) -> Result<Vec<VcfVariant>, SeqvarsError> {
        if variants.is_empty() {
            return Ok(Vec::new());
        }

        let chrom = &variants[0].chromosome;
        if !variants.iter().all(|v| v.chromosome == *chrom) {
            return Err(SeqvarsError::GroupValidation(
                GroupValidationError::DifferentChromosomes,
            ));
        }

        let mut sorted_vars = variants.to_vec();
        sorted_vars.sort_by_key(|v| v.position);

        for window in sorted_vars.windows(2) {
            let v1 = &window[0];
            let v2 = &window[1];
            let v1_end = v1.position + v1.reference.len() as i32 - 1;

            if v1_end >= v2.position {
                return Err(SeqvarsError::GroupValidation(
                    GroupValidationError::OverlappingVariants(v1.position, v2.position),
                ));
            }
        }

        Ok(sorted_vars)
    }

    /// Fetches all relevant transcripts for the group of variants.
    fn get_transcripts_for_variant_group(
        &self,
        sorted_vars: &[VcfVariant],
        chrom_acc: &str,
    ) -> Result<Vec<TxForRegionRecord>, SeqvarsError> {
        let min_pos = sorted_vars.first().unwrap().position;
        let max_pos = sorted_vars.last().unwrap().position
            + sorted_vars.last().unwrap().reference.len() as i32
            - 1;

        let mut txs = self
            .provider
            .get_tx_for_region(
                chrom_acc,
                ALT_ALN_METHOD,
                min_pos - PADDING,
                max_pos + PADDING,
            )
            .map_err(|e| SeqvarsError::Provider(e.to_string()))?;
        txs.sort_by(|a, b| a.tx_ac.cmp(&b.tx_ac));
        Ok(txs)
    }

    /// Evaluates a group of variants for a specific transcript.
    /// Returns a tuple of projections for exonic variants (used for cDNA sequence assembly)
    /// and the combined baseline consequences of all variants in the group.
    fn get_group_projections(
        &self,
        sorted_vars: &[VcfVariant],
        tx: &Transcript,
        chrom_acc: &str,
        transcript_biotype: TranscriptBiotype,
    ) -> Result<Option<(Vec<HgvsProjectionContext>, Consequences)>, SeqvarsError> {
        let alignment = tx.genome_alignments.first().unwrap();
        let strand = Strand::try_from(alignment.strand).expect("invalid strand");

        let splice_csqs = Consequence::SpliceAcceptorVariant
            | Consequence::SpliceDonorVariant
            | Consequence::SpliceRegionVariant
            | Consequence::SplicePolypyrimidineTractVariant
            | Consequence::SpliceDonorRegionVariant
            | Consequence::SpliceDonorFifthBaseVariant
            | Consequence::ExonicSpliceRegionVariant;

        let mut projections = Vec::new();
        let mut group_consequences = Consequences::empty();

        for var in sorted_vars {
            let var_g = Self::get_var_g(var, chrom_acc);
            let (var_start, var_end) = Self::get_var_start_end(&var_g);
            let (tx_loc, tx_csqs) =
                self.determine_transcript_context(alignment, strand, &var_g, var_start, var_end);

            group_consequences |= tx_csqs;

            if tx_csqs.intersects(splice_csqs) {
                tracing::warn!(
                    "Phased variant {:?} affects splicing on {}. Skipping multi-variant assembly.",
                    var,
                    tx.id
                );
                return Ok(None);
            }

            if tx_loc.is_intronic && !tx_loc.is_exonic {
                tracing::debug!(
                    "Skipping purely intronic variant {:?} for sequence assembly on {}",
                    var,
                    tx.id
                );
                continue;
            }

            if !tx_loc.is_exonic {
                return Ok(None);
            }

            let proj = self.project_hgvs(&var_g, &var_g, tx, transcript_biotype)?;
            if proj.n.is_none() || proj.c.is_none() || !proj.is_within_cds_bounds() {
                return Ok(None);
            }

            projections.push(proj);
        }

        Ok(Some((projections, group_consequences)))
    }

    /// Assembles the compound `delins` sequence and constructs the final annotation field.
    fn predict_multiple_for_transcript(
        &self,
        sorted_vars: &[VcfVariant],
        sorted_originals: &[VcfVariant],
        tx_record: &TxForRegionRecord,
        chrom_acc: &str,
        hgvs_g: Option<String>,
    ) -> Result<Option<TxAnnotation>, SeqvarsError> {
        let tx = match self.provider.get_tx(&tx_record.tx_ac) {
            Some(t) => t,
            None => return Ok(None),
        };

        let transcript_biotype = TranscriptBiotype::try_from(tx.biotype).unwrap();
        if transcript_biotype != TranscriptBiotype::Coding || tx.start_codon.is_none() {
            return Ok(None);
        }

        let (projections, base_group_consequences) =
            match self.get_group_projections(sorted_vars, tx, chrom_acc, transcript_biotype)? {
                Some((p, csqs)) => (p, csqs),
                None => return Ok(None),
            };

        if projections.is_empty() {
            // All variants in this group were purely intronic and thus skipped.
            return Ok(None);
        }

        let ref_data = match hgvs::mapper::altseq::ref_transcript_data_cached(
            self.provider.clone(),
            &tx.id,
            None,
        ) {
            Ok(r) => r,
            Err(_) => return Ok(None),
        };

        let vars_n = projections
            .iter()
            .filter_map(|proj| proj.n.as_ref())
            .collect_vec();
        let Some(alt_seq) = apply_n_edits(&ref_data.transcript_sequence, &vars_n) else {
            return Ok(None);
        };

        let n_locs = vars_n.iter().filter_map(|var_n| match var_n {
            HgvsVariant::TxVariant { loc_edit, .. } => Some(loc_edit.loc.inner()),
            _ => None,
        });
        let (Some(n_min), Some(n_max)) = (
            n_locs.clone().map(|loc| loc.start.base).min(),
            n_locs.map(|loc| loc.end.base).max(),
        ) else {
            return Ok(None);
        };

        let tx_len = tx.tx_len();
        let total_delta = alt_seq.len() as i32 - ref_data.transcript_sequence.len() as i32;

        let new_length = (n_max - n_min + 1) + total_delta;
        if new_length < 0 {
            tracing::warn!("Calculated new_length {} is negative.", new_length);
            return Ok(None);
        }
        let start_idx = (n_min - 1) as usize;
        let end_idx = start_idx + new_length as usize;

        if start_idx > alt_seq.len() || end_idx > alt_seq.len() {
            tracing::warn!(
                "Slice index out of bounds: start={} end={} > {} (alt_seq length). Cannot assemble variant.",
                start_idx,
                end_idx,
                alt_seq.len()
            );
            return Ok(None);
        }

        let new_substring = &alt_seq[start_idx..end_idx];
        let ref_substring = &ref_data.transcript_sequence[(n_min - 1) as usize..n_max as usize];

        let compound_var_n = HgvsVariant::TxVariant {
            accession: projections[0].n.as_ref().unwrap().accession().clone(),
            gene_symbol: projections[0].n.as_ref().unwrap().gene_symbol().clone(),
            loc_edit: hgvs::parser::TxLocEdit {
                loc: Mu::Certain(hgvs::parser::TxInterval {
                    start: hgvs::parser::TxPos {
                        base: n_min,
                        offset: None,
                    },
                    end: hgvs::parser::TxPos {
                        base: n_max,
                        offset: None,
                    },
                }),
                edit: Mu::Certain(NaEdit::RefAlt {
                    reference: ref_substring.to_string(),
                    alternative: new_substring.to_string(),
                }),
            },
        };

        let compound_var_c = match self.mapper.n_to_c(&compound_var_n) {
            Ok(c) => c,
            Err(e) => {
                tracing::debug!("Failed to map compound n_loc to c_loc: {}", e);
                return Ok(None);
            }
        };

        let compound_var_n = HgvsVariant::TxVariant {
            accession: projections[0].n.as_ref().unwrap().accession().clone(),
            gene_symbol: projections[0].n.as_ref().unwrap().gene_symbol().clone(),
            loc_edit: hgvs::parser::TxLocEdit {
                loc: Mu::Certain(hgvs::parser::TxInterval {
                    start: hgvs::parser::TxPos {
                        base: n_min,
                        offset: None,
                    },
                    end: hgvs::parser::TxPos {
                        base: n_max,
                        offset: None,
                    },
                }),
                edit: Mu::Certain(NaEdit::RefAlt {
                    reference: ref_substring.to_string(),
                    alternative: new_substring.to_string(),
                }),
            },
        };

        let compound_var_p = self.safe_project_c_to_p(&compound_var_c)?;

        let compound_proj = HgvsProjectionContext {
            n: Some(compound_var_n.clone()),
            c: Some(compound_var_c.clone()),
            p: compound_var_p,
            region: self.variant_location(&compound_var_c),
            unshifted_c: None,
        };

        let tlc = TranscriptLocationContext {
            rank: Rank { ord: 1, total: 1 }, // dummy rank since we potentially span multiple exons
            distance: Some(0),
            is_exonic: true,
            is_intronic: false,
            is_upstream: false,
            is_downstream: false,
        };

        let c_ctx = self.analyze_transcript_consequences(
            &compound_proj,
            tx,
            tx_record,
            &tlc,
            tx_len,
            transcript_biotype,
        )?;

        let mut consequences =
            c_ctx.cds_consequences | c_ctx.protein_consequences | base_group_consequences;
        self.consequences_fix_special_cases(
            &mut consequences,
            c_ctx.cds_consequences,
            c_ctx.protein_consequences,
            &compound_proj,
            tx.is_incomplete_5p(),
        );

        if self.config.vep_consequence_terms {
            self.adjust_vep_terms(&mut consequences, Some(&compound_proj));
        }
        if consequences.is_empty() {
            consequences |= Consequence::GeneVariant;
        }

        let consequences_vec = consequences.iter().collect_vec();
        let putative_impact = (*consequences_vec.first().unwrap()).into();

        let hgvs_n = compound_proj.n.as_ref().map(|n| {
            format!("{}", NoRef(n))
                .split(':')
                .nth(1)
                .unwrap()
                .to_owned()
        });
        let hgvs_c = Some(
            format!("{}", NoRef(compound_proj.c.as_ref().unwrap()))
                .split(':')
                .nth(1)
                .unwrap()
                .to_owned(),
        );
        let hgvs_p = compound_proj
            .p
            .as_ref()
            .map(|p| format!("{}", p).split(':').nth(1).unwrap().to_owned());

        let ref_alts = sorted_originals
            .iter()
            .map(|v| v.reference.clone())
            .collect();
        let alt_alts = sorted_originals
            .iter()
            .map(|v| v.alternative.clone())
            .collect();

        let strand = match Strand::try_from(tx.genome_alignments.first().unwrap().strand) {
            Ok(Strand::Plus) => 1,
            Ok(Strand::Minus) => -1,
            _ => 0,
        };

        let feature_tags = tx
            .tags
            .iter()
            .map(|tag| TranscriptTag::try_from(*tag).expect("invalid transcript tag"))
            .filter(|tag| !matches!(tag, TranscriptTag::EnsemblGraft))
            .filter_map(|t| {
                if !matches!(FeatureTag::from(t), FeatureTag::Other(_)) {
                    Some(FeatureTag::from(t))
                } else {
                    None
                }
            })
            .collect_vec();

        let ann = AnnField {
            allele: Allele::Grouped(GroupedAlleles {
                references: ref_alts,
                alternatives: alt_alts,
            }),
            consequences: consequences_vec,
            putative_impact,
            gene_symbol: tx.gene_symbol.clone(),
            gene_id: tx.gene_id.clone(),
            feature_type: FeatureType::SoTerm {
                term: SoFeature::Transcript,
            },
            feature_id: tx.id.clone(),
            feature_biotype: vec![FeatureBiotype::Coding],
            feature_tags,
            rank: None,
            hgvs_g,
            hgvs_n,
            hgvs_c,
            hgvs_p,
            cdna_pos: c_ctx.cdna_pos,
            cds_pos: c_ctx.cds_pos,
            protein_pos: c_ctx.protein_pos,
            strand,
            distance: Some(0),
            messages: None,
            custom_fields: BTreeMap::new(),
        };

        Ok(Some(TxAnnotation {
            context: self
                .tx_context(Some(tx), Some(compound_var_c), Some(compound_var_n))
                .with_ref_data(ref_data),
            ann,
        }))
    }

    /// The region of the transcript that the c. variant `var_c` changes, as hgvs-rs `c_to_p`
    /// decides it, e.g., `FivePrimeUtr` for `c.1dup`. With `shift_over_boundary`, `c_to_p` can
    /// use another placement of `var_c`, so the region can differ then.
    ///
    /// Returns `None` if hgvs-rs has no reference data for the transcript.
    fn variant_location(&self, var_c: &HgvsVariant) -> Option<VariantLocation> {
        let ref_data = hgvs::mapper::altseq::ref_transcript_data_cached(
            self.provider.clone(),
            &var_c.accession().value,
            None,
        )
        .ok()?;
        AltSeqBuilder {
            ins_at_boundary_is_intronic: self.mapper.config.ins_at_boundary_is_intronic,
            ..AltSeqBuilder::new(var_c.clone(), &ref_data)
        }
        .variant_location()
        .ok()
    }

    /// An `AltSeqBuilder` for `var_c` with the settings that `c_to_p` of `self.mapper` uses.
    fn alt_seq_builder<'a>(
        &self,
        var_c: &HgvsVariant,
        ref_data: &'a RefTranscriptData,
    ) -> AltSeqBuilder<'a> {
        AltSeqBuilder {
            ins_at_boundary_is_intronic: self.mapper.config.ins_at_boundary_is_intronic,
            ..AltSeqBuilder::new(var_c.clone(), ref_data)
        }
    }

    /// Safely projects a CDS variant to a Protein variant, gracefully catching
    /// and swallowing expected incomplete-transcript errors as `None`.
    fn safe_project_c_to_p(
        &self,
        var_c: &HgvsVariant,
    ) -> Result<Option<HgvsVariant>, SeqvarsError> {
        self.mapper
            .variant_mapper()
            .c_to_p(var_c, None)
            .map_or_else(
                |e| {
                    if matches!(e, Error::DuplicationReferenceMissing(_)) {
                        // The intron bases of the dup are unknown, so is its protein change.
                        return Ok(Some(HgvsVariant::ProtVariant {
                            accession: var_c.accession().clone(),
                            gene_symbol: var_c.gene_symbol().clone(),
                            loc_edit: ProtLocEdit::Unknown,
                        }));
                    }
                    if matches!(
                        e,
                        Error::TranscriptLengthInvalid(_, _)
                            | Error::CannotConvertIntervalEnd(_)
                            | Error::MultipleAAVariants
                    ) {
                        tracing::debug!("c_to_p failed gracefully (typed error): {}", e);
                        return Ok(None);
                    }

                    let err_str = e.to_string();
                    if err_str.contains("does not contain a stop codon")
                        || err_str.contains("multiple of 3")
                        || err_str.contains("multiple of three")
                        || err_str.contains("out of bound")
                        || err_str.contains("outside of sequence bounds")
                    {
                        tracing::debug!(
                            "c_to_p failed gracefully (nested error string): {}",
                            err_str
                        );
                        Ok(None)
                    } else {
                        Err(SeqvarsError::HgvsProjection(format!(
                            "c_to_p mapping failed: {}",
                            e
                        )))
                    }
                },
                |v| Ok(Some(v)),
            )
    }
}

/// Whether the positions from `start` to `end` lie within c.1 to c.3. An intronic position such as
/// `c.1+50` is not in the start codon, also if the codon is split by that intron.
fn within_first_codon(start: CdsPos, end: CdsPos) -> bool {
    [&start, &end]
        .iter()
        .all(|pos| pos.cds_from == CdsFrom::Start && pos.offset.unwrap_or_default() == 0)
        && start.base >= 1
        && end.base <= 3
}

/// Whether the insertion or dup `var_c` inserts a number of bases that is not a multiple of 3.
fn inserts_frameshift(var_c: &HgvsVariant) -> bool {
    let HgvsVariant::CdsVariant { loc_edit, .. } = var_c else {
        return false;
    };
    let loc = loc_edit.loc.inner();
    let len = match loc_edit.edit.inner() {
        NaEdit::Ins { alternative } => alternative.len(),
        NaEdit::Dup { reference } if !reference.is_empty() => reference.len(),
        NaEdit::Dup { .. }
            if loc.start.cds_from == loc.end.cds_from
                && loc.start.offset.unwrap_or_default() == 0
                && loc.end.offset.unwrap_or_default() == 0 =>
        {
            loc.start.base.abs_diff(loc.end.base) as usize + 1
        }
        _ => return false,
    };
    len % 3 != 0
}

fn is_conservative_cds_variant(var_c: &HgvsVariant) -> bool {
    match var_c {
        HgvsVariant::CdsVariant { loc_edit, .. } => {
            // Handle the cases where the variant touches the start or stop codon based on `var_c`
            // coordinates. The cases where the start/stop codon is touched by the variant
            // directly is handled elsewhere based on the `var_p` prediction.
            let loc = loc_edit.loc.inner();
            let start_base = loc.start.base;
            let start_cds_from = loc.start.cds_from;
            let end_base = loc.end.base;
            let end_cds_from = loc.end.cds_from;
            // The range is "conservative" (regarding deletions and insertions) if
            // it does not start or end within codons.
            start_cds_from == CdsFrom::Start
                && end_cds_from == CdsFrom::Start
                && start_base % 3 == 1
                && (end_base + 1) % 3 == 1
        }
        _ => panic!("Expected CdsVariant, got {:#?}", var_c),
    }
}

/// Whether deleting `tx_seq[del]` leaves the CDS that starts at `cds_start` intact.
///
/// `del` must start at or after `cds_start`.  The CDS stays intact if it then starts
/// `del.len()` bases earlier, i.e. if the same bases could be deleted from the 5' UTR.
/// VEP's `_ins_del_start_altered` checks the same.
///
/// Only VEP terms use this check. With default terms, the region from hgvs-rs decides.
fn deletion_keeps_cds(tx_seq: &str, cds_start: usize, del: std::ops::Range<usize>) -> bool {
    let Some(new_cds_start) = cds_start.checked_sub(del.len()) else {
        return false;
    };
    // Behind the deletion, the edited sequence equals the reference.
    match (
        tx_seq.get(new_cds_start..del.start),
        tx_seq.get(cds_start..del.end),
    ) {
        (Some(new), Some(old)) => new == old,
        _ => false,
    }
}

#[inline]
fn overlaps(start_a: i32, end_a: i32, start_b: i32, end_b: i32) -> bool {
    (start_a < end_b) && (end_a > start_b)
}

impl ConsequencePredictor {
    /// Return data version string (if set).
    pub fn data_version(&self) -> Option<String> {
        self.provider.as_ref().tx_seq_db.version.clone()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::annotate::cli::{
        PickCriterion, PickUnit, PredictorSettings, TagCriterion, TranscriptSettings,
    };
    use crate::annotate::seqvars::consequence::ConfigBuilder;
    use crate::annotate::seqvars::consequence::load_tx_db;
    use crate::annotate::seqvars::consequence::plugin::{SequencePlugin, TranscriptPlugin};
    use crate::annotate::seqvars::consequence::terms::{
        ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_TX_SEQ_ALT, ANN_TX_SEQ_REF,
    };
    use crate::annotate::seqvars::{
        Args, AsyncAnnotatedVariantWriter, OutputFormat, run_with_writer,
    };
    use crate::common::noodles::{NoodlesVariantReader, open_variant_reader, open_variant_writer};
    use crate::db::transcripts::create::models::Reason;
    use crate::pbs::txs::ExonAlignment;
    use csv::ReaderBuilder;
    use enumflags2::BitFlags;
    use futures::TryStreamExt;
    use insta::assert_yaml_snapshot;
    use noodles::vcf::variant::Record as NoodlesRecord;
    use noodles::vcf::variant::record_buf::info::field::Value;
    use noodles::vcf::variant::record_buf::info::field::value::Array;
    use pretty_assertions::assert_eq;
    use serde::Deserialize;
    use std::collections::BTreeMap;
    use std::path::{Path, PathBuf};
    use std::str::FromStr;
    use std::{fs::File, io::BufReader, io::Write};
    use tempfile::NamedTempFile;

    #[test]
    fn test_sync() {
        fn is_sync<T: Sync>() {}
        is_sync::<super::ConsequencePredictor>();
    }

    /// One transcript per gene, MANE Plus Clinical first, then MANE Select, then the longest.
    fn pick_mane_by_gene() -> TranscriptSettings {
        TranscriptSettings {
            pick: Some(PickUnit::Gene),
            pick_order: vec![
                PickCriterion::ManePlusClinicalBackport,
                PickCriterion::ManeSelectBackport,
                PickCriterion::Length,
            ],
            ..Default::default()
        }
    }

    /// One coding transcript whose exon alignments have genome-only (`I`) and
    /// transcript-only (`D`) bases:
    ///
    /// ```text
    /// exon 1: g.1001_1012 (12 bp), CIGAR 6=2I4=, n.1_10
    /// exon 2: g.2001_2020 (20 bp), CIGAR 10=3D10=, n.11_33
    /// ```
    ///
    /// The genomic exon lengths add up to 32, but the transcript has 33 bases. The CDS is
    /// n.3_29. Like `db create` does for a complete CDS, the stored sequence ends with three
    /// extra `A`.
    fn tx_db_with_indels_in_exon_cigars() -> crate::pbs::txs::TxSeqDatabase {
        use crate::pbs::txs::{
            ExonAlignment, GeneToTxId, SequenceDb, SourceVersion, TranscriptDb, TxSeqDatabase,
        };

        let exon = |ord, alt_start_i, alt_end_i, tx_start, tx_end, cigar: &str| ExonAlignment {
            alt_start_i,
            alt_end_i,
            ord,
            alt_cds_start_i: Some(tx_start),
            alt_cds_end_i: Some(tx_end),
            cigar: cigar.into(),
        };

        TxSeqDatabase {
            tx_db: Some(TranscriptDb {
                transcripts: vec![Transcript {
                    id: "NM_000000.1".into(),
                    gene_symbol: "TEST".into(),
                    gene_id: "HGNC:0".into(),
                    biotype: TranscriptBiotype::Coding.into(),
                    protein: Some("NP_000000.1".into()),
                    start_codon: Some(2),
                    stop_codon: Some(29),
                    genome_alignments: vec![GenomeAlignment {
                        genome_build: "grch37".into(),
                        contig: "NC_000001.10".into(),
                        cds_start: Some(1002),
                        cds_end: Some(2016),
                        strand: Strand::Plus.into(),
                        exons: vec![
                            exon(0, 1000, 1012, 1, 10, "6=2I4="),
                            exon(1, 2000, 2020, 11, 33, "10=3D10="),
                        ],
                        ..Default::default()
                    }],
                    ..Default::default()
                }],
                gene_to_tx: vec![GeneToTxId {
                    gene_id: "HGNC:0".into(),
                    tx_ids: vec!["NM_000000.1".into()],
                    ..Default::default()
                }],
            }),
            seq_db: Some(SequenceDb {
                aliases: vec!["NM_000000.1".into()],
                aliases_idx: vec![0],
                seqs: vec!["GG".to_string() + "ATGGCCAAAGGGCCCTTTGGGAAATAA" + "CCCC" + "AAA"],
            }),
            source_version: vec![SourceVersion {
                assembly: "grch37".into(),
                ..Default::default()
            }],
            ..Default::default()
        }
    }

    /// A change of the first codon gives no start terms if the CDS lacks its start codon. c.2T>A
    /// changes the ATG at n.3_5.
    #[rstest::rstest]
    #[case::complete(false, Consequence::StartLost)]
    #[case::incomplete_5p(true, Consequence::CodingSequenceVariant)]
    fn first_codon_of_cds_without_start_codon(
        #[case] incomplete_5p: bool,
        #[case] expected: Consequence,
    ) -> Result<(), anyhow::Error> {
        use crate::db::transcripts::create::models::Reason;

        let mut tx_db = tx_db_with_indels_in_exon_cigars();
        if incomplete_5p {
            let tx = &mut tx_db.tx_db.as_mut().unwrap().transcripts[0];
            tx.filter_reason =
                Some(enumflags2::BitFlags::from(Reason::FivePrimeEndTruncated).bits());
        }
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(provider, Default::default());

        let anns = predictor
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1004,
                reference: "T".into(),
                alternative: "A".into(),
            })?
            .unwrap();

        assert_eq!(anns[0].hgvs_c.as_deref(), Some("c.2T>A"));
        assert!(anns[0].consequences.contains(&expected), "{:?}", anns[0]);
        assert_eq!(
            anns[0].consequences.contains(&Consequence::StartLost),
            !incomplete_5p
        );
        Ok(())
    }

    /// The cDNA total is the transcript length, both for single and for phased variants.
    #[test]
    fn annotate_totals_with_indels_in_exon_cigars() -> Result<(), anyhow::Error> {
        let provider = Arc::new(MehariProvider::new(
            tx_db_with_indels_in_exon_cigars(),
            None::<PathBuf>,
            true,
        ));
        let predictor = ConsequencePredictor::new(provider, Default::default());

        let var = |position, reference: &str, alternative: &str| VcfVariant {
            chromosome: "1".into(),
            position,
            reference: reference.into(),
            alternative: alternative.into(),
        };
        // c.14C>A, p.Pro5His
        let single = predictor.predict(&var(2006, "C", "A"))?.unwrap();
        // c.[14C>A;16T>G], p.Pro5_Phe6delinsHisVal
        let phased = predictor
            .predict_multiple(&[var(2006, "C", "A"), var(2008, "T", "G")])?
            .unwrap();

        for ann in [&single[0], &phased[0]] {
            let pos = |ord, total| Some(Pos { ord, total });
            assert_eq!(ann.feature_id, "NM_000000.1");
            assert_eq!(ann.cdna_pos, pos(16, Some(33)));
            assert_eq!(ann.cds_pos, pos(14, Some(27)));
            assert_eq!(ann.protein_pos, pos(5, Some(8)));
        }

        Ok(())
    }

    /// Like in VEP, the protein total does not count the stop codon. A variant in the stop
    /// codon therefore has a protein position one past the total. A transcript without a stop
    /// codon counts all of its codons.
    #[test]
    fn annotate_protein_total_excludes_stop_codon() -> Result<(), anyhow::Error> {
        let predictor = |db| {
            let provider = Arc::new(MehariProvider::new(db, None::<PathBuf>, true));
            ConsequencePredictor::new(provider, Default::default())
        };
        let var = |position, reference: &str, alternative: &str| VcfVariant {
            chromosome: "1".into(),
            position,
            reference: reference.into(),
            alternative: alternative.into(),
        };
        let pos = |ord, total| Some(Pos { ord, total });

        // c.26A>C in the stop codon TAA
        let ann = predictor(tx_db_with_indels_in_exon_cigars())
            .predict(&var(2015, "A", "C"))?
            .unwrap();
        assert_eq!(ann[0].protein_pos, pos(9, Some(8)));

        // c.22_24del deletes the last codon before the stop codon: no stop_gained.
        let ann = predictor(tx_db_with_indels_in_exon_cigars())
            .predict(&var(2010, "TAAA", "T"))?
            .unwrap();
        assert_eq!(ann[0].hgvs_p.as_deref(), Some("p.Lys8Ter"));
        assert_eq!(
            ann[0].consequences,
            vec![Consequence::ConservativeInframeDeletion]
        );

        let mut db = tx_db_with_indels_in_exon_cigars();
        if let Some(tx_db) = db.tx_db.as_mut() {
            tx_db.transcripts[0].filter_reason =
                Some(BitFlags::from(Reason::MissingStopCodon).bits());
        }
        // c.14C>A, p.Pro5His
        let ann = predictor(db).predict(&var(2006, "C", "A"))?.unwrap();
        assert_eq!(ann[0].protein_pos, pos(5, Some(9)));

        Ok(())
    }

    #[rstest::rstest]
    #[case("17:41197701:G:C", 0)] // exonic
    #[case("17:41196309:G:C", -3)] // 3bp 3' upstream
    #[case("17:41196310:G:C", -2)] // 2bp 3' upstream
    #[case("17:41196311:G:C", -1)] // 1bp 3' upstream
    #[case("17:41196312:G:C", 0)] // ex. 3' UTR
    #[case("17:41196313:G:C", 0)] // ex. 3' UTR
    #[case("17:41197818:G:C", 0)] // exonic
    #[case("17:41197819:G:C", 0)] // exonic
    #[case("17:41197820:G:C", 1)] // 1bp intronic
    #[case("17:41197821:G:C", 2)] // 2bp intronic
    #[case("17:41197822:G:C", 3)] // 3bp intronic
    #[case("17:41197823:G:C", 4)] // 4bp intronic
    #[case("17:41277379:A:C", 0)] // exonic
    #[case("17:41277380:G:C", 0)] // exonic
    #[case("17:41277381:G:T", 0)] // exonic
    #[case("17:41277382:G:C", 1)] // 1bp upstream
    #[case("17:41277383:A:C", 2)] // 2bp upstream
    #[case("17:41277384:G:C", 3)] // 3bp upstream
    fn annotate_snv_brca1_one_variant(
        #[case] spdi: &str,
        #[case] expected_dist: i32,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", spdi.replace(':', "-"));

        let variant: VcfVariant = spdi.parse()?;

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));

        let predictor = ConsequencePredictor::new(provider, Default::default());

        let res = predictor.predict(&variant)?.unwrap();

        assert_eq!(res.len(), 6);
        insta::assert_yaml_snapshot!(res);
        assert_eq!(res[0].distance, Some(expected_dist), "spdi = {}", spdi);

        Ok(())
    }

    /// `CDS.pos` must only be set for variants within the CDS, not for UTR variants.
    /// GRCh37, BRCA1, NM_007294.4 (MANE, reverse).
    #[rstest::rstest]
    #[case("17:41276119:T:C", "c.-", false)] // exonic, 5' UTR
    #[case("17:41196312:G:C", "c.*", false)] // exonic, 3' UTR
    #[case("17:41276049:A:T", "c.65", true)] // exonic, CDS
    fn annotate_brca1_cds_pos_only_in_cds(
        #[case] spdi: &str,
        #[case] expected_hgvs_c_prefix: &str,
        #[case] expect_cds_pos: bool,
    ) -> Result<(), anyhow::Error> {
        let spdi = spdi.split(':').map(|s| s.to_string()).collect::<Vec<_>>();

        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(provider, Default::default());

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].clone(),
                position: spdi[1].parse()?,
                reference: spdi[2].clone(),
                alternative: spdi[3].clone(),
            })?
            .unwrap();

        let ann = res
            .iter()
            .find(|ann| ann.feature_id == "NM_007294.4")
            .ok_or_else(|| anyhow::anyhow!("NM_007294.4 not annotated"))?;
        let hgvs_c = ann
            .hgvs_c
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("no hgvs_c"))?
            .to_string();
        assert!(
            hgvs_c.starts_with(expected_hgvs_c_prefix),
            "unexpected hgvs_c: {hgvs_c}"
        );
        assert_eq!(ann.cds_pos.is_some(), expect_cds_pos, "hgvs_c = {hgvs_c}");

        Ok(())
    }

    /// Test some intron specific variants, via the annotated consequences.
    /// GRCh37, BRCA1, NM_007294.4 (MANE, reverse).
    /// The order of the consequences is important: ordered by severity, descending.
    /// cf Consequences enum ordering.
    #[rstest::rstest]
    #[case("17:41197820:G:T", 1, vec![Consequence::SpliceAcceptorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 1bp intronic
    #[case("17:41197821:A:C", 2, vec![Consequence::SpliceAcceptorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 2bp intronic
    #[case("17:41197822:C:A", 3, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 3bp intronic
    #[case("17:41197823:C:A", 4, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 4bp intronic
    #[case("17:41197824:T:G", 5, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 5bp intronic
    #[case("17:41197825:C:A", 6, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 6bp intronic
    #[case("17:41197835:T:G", 16, vec![Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 16bp intronic
    #[case("17:41197836:G:A", 17, vec![Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 17bp intronic
    #[case("17:41197837:G:A", 18, vec![Consequence::CodingTranscriptIntronVariant]
    )] // 18bp intronic
    #[case("17:41199660:G:T", 0, vec![Consequence::MissenseVariant, Consequence::ExonicSpliceRegionVariant]
    )] // exonic
    #[case("17:41199659:G:T", -1, vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -1bp intronic
    #[case("17:41199658:T:G", -2, vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -2bp intronic
    #[case("17:41199657:G:T", -3, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -3bp intronic
    #[case("17:41199656:A:C", -4, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -4bp intronic
    #[case("17:41199655:G:T", -5, vec![Consequence::SpliceDonorFifthBaseVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -5bp intronic
    #[case("17:41199654:G:T", -6, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -6bp intronic
    #[case("17:41199653:T:G", -7, vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -7bp intronic
    #[case("17:41199652:G:T", -8, vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -8bp intronic
    #[case("17:41199651:C:A", -9, vec![Consequence::CodingTranscriptIntronVariant]
    )] // -9bp intronic
    fn annotate_snv_brca1_csq(
        #[case] spdi: &str,
        #[case] expected_dist: i32,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", spdi.replace(':', "-"));

        let variant: VcfVariant = spdi.parse()?;

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));

        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .transcript_settings(pick_mane_by_gene())
                .build()?,
        );

        let res = predictor.predict(&variant)?.unwrap();

        assert_eq!(res.len(), 1);
        assert_eq!(res[0].feature_id, "NM_007294.4");
        assert_eq!(res[0].distance, Some(expected_dist), "spdi = {}", spdi);
        assert_eq!(res[0].consequences, expected_csqs, "spdi = {}", spdi);
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Test some intron specific variants, via the annotated consequences.
    /// GRCh37, OPA1, NM_130837.3 (MANE, forward).
    /// The order of the consequences is important: ordered by severity, descending.
    /// cf Consequences enum ordering.
    #[rstest::rstest]
    #[case("3:193332512:T:G", 0, vec![Consequence::MissenseVariant, Consequence::ExonicSpliceRegionVariant]
    )] // exonic
    #[case("3:193332511:G:T", -1, vec![Consequence::SpliceAcceptorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -1bp intronic
    #[case("3:193332510:A:G", -2, vec![Consequence::SpliceAcceptorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -2bp intronic
    #[case("3:193332509:C:T", -3, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant,  Consequence::CodingTranscriptIntronVariant]
    )] // -3bp intronic
    #[case("3:193332508:T:C", -4, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant,  Consequence::CodingTranscriptIntronVariant]
    )] // -4bp intronic
    #[case("3:193332507:T:C", -5, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant,  Consequence::CodingTranscriptIntronVariant]
    )] // -5bp intronic
    #[case("3:193332506:T:C", -6, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -6bp intronic
    #[case("3:193332505:C:G", -7, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -7bp intronic
    #[case("3:193332504:T:C", -8, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -8bp intronic
    #[case("3:193332503:T:A", -9, vec![Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -9bp intronic
    #[case("3:193332831:G:T", 1, vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 1bp intronic
    #[case("3:193332832:T:C", 2, vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 2bp intronic
    #[case("3:193332833:G:A", 3, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 3bp intronic
    #[case("3:193332834:A:C", 4, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 4bp intronic
    #[case("3:193332835:A:T", 5, vec![Consequence::SpliceDonorFifthBaseVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 5bp intronic
    #[case("3:193332836:C:A", 6, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 6bp intronic
    #[case("3:193332837:T:G", 7, vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 7bp intronic
    #[case("3:193332838:T:G", 8, vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 8bp intronic
    #[case("3:193332839:G:A", 9, vec![Consequence::CodingTranscriptIntronVariant])] // 9bp intronic
    #[case("3:193332846:A:G", 16, vec![Consequence::CodingTranscriptIntronVariant])] // 16bp intronic
    #[case("3:193332847:G:A", 17, vec![Consequence::CodingTranscriptIntronVariant])] // 17bp intronic
    #[case("3:193332848:T:A", 18, vec![Consequence::CodingTranscriptIntronVariant])] // 18bp intronic
    fn annotate_snv_opa1_csq(
        #[case] spdi: &str,
        #[case] expected_dist: i32,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", spdi.replace(':', "-"));

        let variant: VcfVariant = spdi.parse()?;

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));

        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .transcript_settings(pick_mane_by_gene())
                .build()?,
        );

        let res = predictor.predict(&variant)?.unwrap();

        assert_eq!(res.len(), 1);
        assert_eq!(res[0].feature_id, "NM_130837.3");
        assert_eq!(res[0].distance, Some(expected_dist), "spdi = {}", spdi);
        assert_eq!(res[0].consequences, expected_csqs, "spdi = {}", spdi);
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Introns are numbered in transcript order on both strands.
    #[rstest::rstest]
    #[case("3:193332511:G:T", "NM_130837.3", 1, 30)] // c.33-1, plus strand
    #[case("3:193332831:G:T", "NM_130837.3", 2, 30)] // c.351+1
    #[case("17:41197820:G:T", "NM_007294.4", 22, 22)] // c.5468-1, minus strand
    #[case("17:41199659:G:T", "NM_007294.4", 22, 22)] // c.5467+1
    fn annotate_intron_rank(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] ord: i32,
        #[case] total: i32,
    ) -> Result<(), anyhow::Error> {
        let spdi = spdi.split(':').map(|s| s.to_string()).collect::<Vec<_>>();

        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(provider, Default::default());

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].clone(),
                position: spdi[1].parse()?,
                reference: spdi[2].clone(),
                alternative: spdi[3].clone(),
            })?
            .unwrap();
        let ann = res.iter().find(|ann| ann.feature_id == tx_id).unwrap();

        assert_eq!(ann.rank, Some(Rank { ord, total }));

        Ok(())
    }

    #[rstest::rstest]
    #[case("3:193311167:ATGT:T", vec![Consequence::StartLost, Consequence::ConservativeInframeDeletion]
    )]
    #[case("3:193311170:TGGC:C", vec![Consequence::ConservativeInframeDeletion])]
    #[case("3:193311170:TGGCG:G", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation]
    )]
    #[case("3:193311180:GTCG:G", vec![Consequence::DisruptiveInframeDeletion])]
    #[case("3:193409910:GAAA:G", vec![Consequence::ConservativeInframeDeletion])]
    #[case("3:193409913:ATAA:A", vec![Consequence::StopLost, Consequence::FeatureElongation, Consequence::ConservativeInframeDeletion]
    )]
    fn annotate_del_opa1_csqs(
        #[case] spdi: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", spdi.replace(':', "-"));

        let variant: VcfVariant = spdi.parse()?;

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));

        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .transcript_settings(pick_mane_by_gene())
                .build()?,
        );

        let res = predictor.predict(&variant)?.unwrap();

        assert_eq!(res.len(), 1);
        assert_eq!(res[0].feature_id, "NM_130837.3");
        assert_eq!(res[0].consequences, expected_csqs, "spdi = {}", spdi);
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Indels on Ensembl 108 chr22 transcripts.
    ///
    /// With VEP terms, the expected terms are those of VEP 116.3 with `--shift_3prime 1`. mehari
    /// adds `feature_elongation` and `protein_altering_variant`. At the stop codon, mehari keeps
    /// `frameshift_variant` next to `stop_lost` or `stop_retained_variant`, which VEP drops.
    #[rstest::rstest]
    // `p.Val170Ter`: VEP's changed codons hold no complete codon of the new frame
    #[case("22:19524002:AC:A", "ENST00000403084", true, vec![Consequence::FrameshiftVariant])]
    // `p.Tyr1910Ter`: insertion inside codon 1910
    #[case("22:17791223:T:TC", "ENST00000441493", true, vec![Consequence::StopGained, Consequence::FrameshiftVariant])]
    // Default terms: `frameshift_variant`, plus `stop_gained` if the first changed codon is a stop
    #[case("22:19524002:AC:A", "ENST00000403084", false, vec![Consequence::StopGained, Consequence::FrameshiftVariant])]
    #[case("22:17791223:T:TC", "ENST00000441493", false, vec![Consequence::StopGained, Consequence::FrameshiftVariant])]
    #[case("22:38112210:TCA:T", "ENST00000332509", false, vec![Consequence::StopGained, Consequence::FrameshiftVariant])]
    #[case("22:17191782:T:TTATG", "ENST00000262607", false, vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation])]
    // `p.Asp261AlafsTer2`: the new stop lies in the changed codons
    #[case("22:17191782:T:TTATG", "ENST00000262607", true, vec![Consequence::StopGained, Consequence::FrameshiftVariant])]
    // `p.Tyr790Ter`: deletion across codons 790 and 791
    #[case("22:38112210:TCA:T", "ENST00000332509", true, vec![Consequence::StopGained, Consequence::FrameshiftVariant])]
    // `p.Glu587Ter`: insertion between codons 586 and 587
    #[case("22:20112682:C:CA", "ENST00000252136", true, vec![Consequence::FrameshiftVariant])]
    // `p.Ter85ArgextTer9`
    #[case("22:22895417:CCT:C", "ENST00000531372", true, vec![Consequence::FrameshiftVariant, Consequence::StopLost, Consequence::FeatureElongation])]
    // `p.Ter204TrpextTer73`
    #[case("22:42571197:T:TG", "ENST00000340239", true, vec![Consequence::FrameshiftVariant, Consequence::StopLost, Consequence::FeatureElongation])]
    // `p.Ter704=`: VEP's changed codons end in a partial codon, and the bases at the stop codon
    // still read a stop codon
    #[case("22:45600444:TG:T", "ENST00000327858", true, vec![Consequence::FrameshiftVariant, Consequence::StopRetainedVariant])]
    // `p.Ter204=`: VEP's changed codons end in a partial codon, and VEP does not count the
    // insertion behind the first base of the stop codon as overlapping it
    #[case("22:42571196:T:TTAAC", "ENST00000340239", true, vec![Consequence::FrameshiftVariant])]
    // `p.Ter512=`: the new frame starts with a stop codon
    #[case("22:17181484:C:CT", "ENST00000262607", true, vec![Consequence::FrameshiftVariant, Consequence::StopRetainedVariant])]
    // `p.Ter704AspextTer1`: in-frame insertion in front of the stop codon
    #[case("22:45600443:C:CGAT", "ENST00000327858", true, vec![Consequence::FeatureElongation, Consequence::InframeInsertion])]
    // `p.Cys203Ter` (VEP: `p.Cys203del`): in-frame deletion in front of the stop codon
    #[case("22:42571193:CTGT:C", "ENST00000340239", true, vec![Consequence::InframeDeletion])]
    // in-frame deletion across the stop codon
    #[case("22:22895418:CTCT:C", "ENST00000531372", true, vec![Consequence::StopLost, Consequence::InframeDeletion, Consequence::ProteinAlteringVariant])]
    // `c.*2_*4del`: VEP checks no stop codon for a deletion behind the CDS, also if the input
    // position overlaps the stop codon
    #[case("22:22895421:TGAG:T", "ENST00000531372", true, vec![Consequence::ThreePrimeUtrVariant])]
    // `c.1449_*18del`: VEP builds no altered sequence for a deletion beyond the transcript end,
    // so the stop codon counts as kept
    #[case("22:50525751:GTTTCGCGGCAAAGGAGCTT:G", "ENST00000487577", true, vec![Consequence::FeatureElongation, Consequence::StopRetainedVariant, Consequence::ThreePrimeUtrVariant])]
    // `p.Ter704=`: in-frame insertion inside the stop codon that keeps it
    #[case("22:45600443:C:CTAA", "ENST00000327858", false, vec![Consequence::StopRetainedVariant])]
    #[case("22:45600443:C:CTAA", "ENST00000327858", true, vec![Consequence::InframeInsertion, Consequence::StopRetainedVariant])]
    // `p.Phe703=`: in-frame insertion that keeps codon 703 and forms a new stop codon in front
    // of the old one. VEP's changed codons translate to `F` and `F*H`, so VEP reports no
    // `synonymous_variant`.
    #[case("22:45600442:T:TTTAGCA", "ENST00000327858", false, vec![Consequence::SynonymousVariant])]
    #[case("22:45600442:T:TTTAGCA", "ENST00000327858", true, vec![Consequence::InframeInsertion, Consequence::StopRetainedVariant])]
    // `p.?`: in-frame insertion that creates a stop codon in a CDS without one (`cds_end_NF`).
    // VEP's changed codons keep the last amino acid and continue with the stop codon.
    #[case("22:38140065:C:CTAA", "ENST00000430886", false, vec![Consequence::DisruptiveInframeInsertion])]
    #[case("22:38140065:C:CTAA", "ENST00000430886", true, vec![Consequence::InframeInsertion, Consequence::StopRetainedVariant])]
    // `p.?`: frameshift in a CDS without a stop codon (`cds_end_NF`)
    #[case("22:38140065:C:CAG", "ENST00000430886", false, vec![Consequence::FrameshiftVariant])]
    // `c.498dup` (`p.?`): a dup of the last CDS base, but without a stop codon the CDS goes on
    #[case("22:38140065:C:CC", "ENST00000430886", false, vec![Consequence::FrameshiftVariant])]
    // `c.932_933dup` (`p.?`): the copy lands behind the stop codon
    #[case("22:21469819:T:TAA", "ENST00000432134", false, vec![Consequence::ThreePrimeUtrExonVariant])]
    #[case("22:21469819:T:TAA", "ENST00000432134", true, vec![Consequence::ThreePrimeUtrVariant])]
    // `c.933dup` (`p.?`): the copy lands behind the stop codon. VEP checks the overlap with the
    // stop codon at the input position, in front of the last two stop codon bases.
    #[case("22:21469819:T:TA", "ENST00000432134", true, vec![Consequence::StopRetainedVariant, Consequence::ThreePrimeUtrVariant])]
    // `p.Gln147AlafsTer?`: no stop codon in the new frame (`cds_end_NF` transcript)
    #[case("22:38140124:G:GC", "ENST00000430886", false, vec![Consequence::FrameshiftVariant])]
    // `p.Asp443ProfsTer?`: no stop codon in the new frame (complete transcript)
    #[case(
        "22:50525872:GCCGCTGAGCGCGGGGCCGTC:G",
        "ENST00000487577",
        false,
        vec![Consequence::FrameshiftVariant]
    )]
    // `p.Gln482ProfsTer?`: the new frame reads past the stop codon to the transcript end
    #[case("22:50525774:T:TG", "ENST00000487577", false, vec![Consequence::FrameshiftVariant, Consequence::FrameshiftElongation])]
    fn annotate_indel_csqs(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] vep_consequence_terms: bool,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let spdi = spdi.split(':').map(|s| s.to_string()).collect::<Vec<_>>();

        let tx_path = "tests/data/annotate/db/grch38/GRCh38-ensembl.frameshift-subset.txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .vep_consequence_terms(vep_consequence_terms)
                .build()?,
        );

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].clone(),
                position: spdi[1].parse()?,
                reference: spdi[2].clone(),
                alternative: spdi[3].clone(),
            })?
            .unwrap();

        let ann = res
            .iter()
            .find(|ann| ann.feature_id.starts_with(tx_id))
            .unwrap();
        assert_eq!(
            ann.consequences,
            expected_csqs,
            "spdi = {}, hgvs_p = {:?}",
            spdi.join(":"),
            ann.hgvs_p
        );

        Ok(())
    }

    /// Write chromosome 22 as an indexed FASTA file into `dir`: the windows from `path` and `N`
    /// elsewhere. Each record ID names the 1-based region of its window, e.g., `22:100-200`.
    fn write_chr22_from_windows(path: &str, dir: &Path) -> Result<PathBuf, anyhow::Error> {
        let mut seq = Vec::new();
        for record in bio::io::fasta::Reader::from_file(path)?.records() {
            let record = record?;
            let start = record
                .id()
                .split([':', '-'])
                .nth(1)
                .ok_or_else(|| anyhow::anyhow!("no window start in {}", record.id()))?
                .parse::<usize>()?
                - 1;
            let end = start + record.seq().len();
            if seq.len() < end {
                seq.resize(end, b'N');
            }
            seq[start..end].copy_from_slice(record.seq());
        }
        let fasta = dir.join("chr22.fa");
        let mut file = File::create(&fasta)?;
        file.write_all(b">22\n")?;
        file.write_all(&seq)?;
        file.write_all(b"\n")?;
        std::fs::write(
            dir.join("chr22.fa.fai"),
            format!("22\t{}\t4\t{}\t{}\n", seq.len(), seq.len(), seq.len() + 1),
        )?;
        Ok(fasta)
    }

    /// Annotation of `spdi` on the transcript `tx_id`, with Ensembl 108 chr22 transcripts and
    /// the chr22 reference around the variant, so that mehari shifts indels 3'.
    fn annotate_chr22_window(
        spdi: &str,
        tx_id: &str,
        vep_consequence_terms: bool,
    ) -> Result<AnnField, anyhow::Error> {
        let spdi = spdi.split(':').collect::<Vec<_>>();

        let dir = tempfile::tempdir()?;
        let reference = write_chr22_from_windows(
            "tests/data/annotate/seqvars/placement.chr22-windows.fa",
            dir.path(),
        )?;
        let tx_path = "tests/data/annotate/db/grch38/GRCh38-ensembl.placement-subset.txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            Some(reference),
            false,
            Default::default(),
        ));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .vep_consequence_terms(vep_consequence_terms)
                .build()?,
        );

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].to_string(),
                position: spdi[1].parse()?,
                reference: spdi[2].to_string(),
                alternative: spdi[3].to_string(),
            })?
            .unwrap();

        res.into_iter()
            .find(|ann| ann.feature_id.starts_with(tx_id))
            .ok_or_else(|| anyhow::anyhow!("no annotation on {}", tx_id))
    }

    /// Indels at exon edges on Ensembl 108 chr22 transcripts, with the reference so that mehari
    /// shifts them 3'. The terms come from the equivalent placement that keeps the essential
    /// splice site (and in the CDS the stop codon); the HGVS stays 3'-shifted.
    #[rstest::rstest]
    // deleting the acceptor G equals deleting the first exon G
    #[case("22:37713208:AG:A", "ENST00000644935", false, "c.257del", "p.Gly86AlafsTer33", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // deleting the donor G (`GCCG|GTGAGT`) equals deleting the last exon G
    #[case("22:23772977:CG:C", "ENST00000215743", false, "c.108+1del", "p.?", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    #[case("22:23772977:CG:C", "ENST00000215743", true, "c.108+1del", "p.?", vec![Consequence::FrameshiftVariant, Consequence::SpliceRegionVariant, Consequence::CodingSequenceVariant])]
    // minus strand: deleting donor +1 to +3 equals deleting the last three exon bases
    #[case("22:26481718:TCAC:T", "ENST00000336873", false, "c.41+1_41+3del", "p.?", vec![Consequence::DisruptiveInframeDeletion, Consequence::ExonicSpliceRegionVariant])]
    // minus strand: deleting acceptor -3 and -2 equals an intronic deletion that keeps the AG
    #[case("22:50246093:CTG:C", "ENST00000216271", false, "c.1651-3_1651-2del", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant])]
    // deletion of donor +3 to +9 that keeps the GT
    #[case("22:25763388:AGGTAAGT:A", "ENST00000335473", false, "c.198+3_198+9del", "p.?", vec![Consequence::SpliceDonorFifthBaseVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    // deleting one C of `CCCAG|GT`: a placement in front of the last three exon bases keeps them
    #[case("22:39101499:TC:T", "ENST00000442487", false, "c.416del", "p.Pro139GlnfsTer25", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation])]
    // every placement of the deletion removes donor +1 and +2
    #[case("22:28773775:ATGAG:A", "ENST00000249064", false, "c.238_239+2del", "p.?", vec![Consequence::SpliceDonorVariant, Consequence::ExonicSpliceRegionVariant])]
    // insertion between the last exon base and donor +1: the exon gains a base
    #[case("22:20996112:T:TA", "ENST00000646124", false, "c.2219_2219+1insA", "p.Tyr741LeufsTer41", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    #[case("22:20996112:T:TA", "ENST00000646124", true, "c.2219_2219+1insA", "p.Tyr741LeufsTer41", vec![Consequence::FrameshiftVariant, Consequence::SpliceRegionVariant, Consequence::CodingSequenceVariant])]
    // minus strand: insertion between donor +1 and the last exon base
    #[case("22:32150991:C:CTT", "ENST00000382097", false, "c.493_493+1insAA", "p.His165GlnfsTer19", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // insertion one base past the acceptor, between the first two exon bases
    #[case("22:50522343:G:GCTTTCTC", "ENST00000299821", false, "c.1234_1235insCTTTCTC", "p.Val412AlafsTer5", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // insertion between acceptor -1 and the first CDS base: the 5' UTR gains a base
    #[case("22:28987075:G:GC", "ENST00000402174", false, "c.1-1_1insC", "p.?", vec![Consequence::ExonicSpliceRegionVariant, Consequence::FivePrimeUtrExonVariant])]
    // insertion between the last stop codon base and donor +1: the 3' UTR gains a base
    #[case("22:31879752:G:GA", "ENST00000646998", false, "c.3858_3858+1insA", "p.?", vec![Consequence::ExonicSpliceRegionVariant, Consequence::ThreePrimeUtrExonVariant])]
    // minus strand: dup of the last exon base, in front of donor +1 or inside the exon
    #[case("22:28710005:C:CA", "ENST00000348295", false, "c.846dup", "p.Pro283SerfsTer7", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // dup inside the stop codon that equals one behind it (`TAA` to `TAAAA`)
    #[case("22:21469819:T:TAA", "ENST00000432134", false, "c.932_933dup", "p.?", vec![Consequence::ThreePrimeUtrExonVariant])]
    // frameshift inside the stop codon: every placement changes the stop codon
    #[case("22:42571197:T:TG", "ENST00000340239", false, "c.611dup", "p.Ter204TrpextTer73", vec![Consequence::StopLost, Consequence::FeatureElongation])]
    // dup of donor +1 to +7: the copy lands behind +7 and leaves the donor intact
    #[case("22:31604404:G:GGTAGGAA", "ENST00000400288", false, "c.1977+1_1977+7dup", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    #[case("22:31604404:G:GGTAGGAA", "ENST00000400288", true, "c.1977+1_1977+7dup", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::IntronVariant])]
    // dup of donor +5 and +6: the copy lands behind +6 and leaves donor +3 to +6 intact
    #[case("22:31616877:G:GGT", "ENST00000400288", false, "c.3433+5_3433+6dup", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    // dup of c.16 to c.20: a copy in front of c.16 leaves the last three exon bases intact
    #[case("22:42554033:G:GCCGCA", "ENST00000340239", false, "c.16_20dup", "p.Gly8HisfsTer4", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation])]
    fn annotate_indel_placement_csqs(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] vep_consequence_terms: bool,
        #[case] expected_hgvs_c: &str,
        #[case] expected_hgvs_p: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let ann = annotate_chr22_window(spdi, tx_id, vep_consequence_terms)?;
        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (Some(expected_hgvs_c), Some(expected_hgvs_p), &expected_csqs),
            "spdi = {spdi}"
        );

        Ok(())
    }

    /// An insertion at an exon edge adds its bases to the exon. So it gets the exon number and a
    /// CDS position.
    #[test]
    fn annotate_ins_at_exon_edge_location() -> Result<(), anyhow::Error> {
        let ann = annotate_chr22_window("22:32150991:C:CTT", "ENST00000382097", false)?;
        assert_eq!(ann.hgvs_c.as_deref(), Some("c.493_493+1insAA"));
        assert_eq!(ann.rank, Some(Rank { ord: 6, total: 9 }));
        assert_eq!(ann.distance, Some(0));
        assert_eq!(
            ann.cds_pos,
            Some(Pos {
                ord: 493,
                total: Some(756)
            })
        );

        Ok(())
    }

    /// Dups and insertions at the edges of a coding exon. A dup of the first or last exon base,
    /// or of the intron base next to the exon, equals an insertion of that base at the exon
    /// edge. There, it keeps the splice site and joins the exon, so the frame shifts. An
    /// insertion between the two donor site bases breaks the site instead.
    ///
    /// Exon 2 of ENST00000249064 (plus strand) spans c.186 to c.239 (`CAG|TGT`, `GATGA|GTAAGT`).
    /// On the minus strand, exon 7 of ENST00000348295 ends at c.846 (`CAT|GTAAGT`), and exon 18
    /// of ENST00000216271 starts at c.1651 (`CAG|ATG`).
    #[rstest::rstest]
    // codon 62 CG|T stays, codon 63 GTT (Val) becomes TGT (Cys)
    #[case("22:28773725:T:TT", "ENST00000249064", "c.186dup", "p.Val63CysfsTer6", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 62 becomes CGG (Arg, as CGT), codon 63 TGT (Cys) as for c.186dup
    #[case("22:28773724:G:GG", "ENST00000249064", "c.186-1dup", "p.Val63CysfsTer6", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 80 GA|T (Asp, its T in exon 3) becomes GAA (Glu)
    #[case("22:28773778:A:AA", "ENST00000249064", "c.239dup", "p.Asp80GlufsTer11", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 80 becomes GAG (Glu), the frame as for c.239dup
    #[case("22:28773779:G:GG", "ENST00000249064", "c.239+1dup", "p.Asp80GlufsTer11", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 80 becomes GAC (Asp, as GAT), codon 81 TGT (Cys) becomes TTG (Leu)
    #[case("22:28773778:A:AC", "ENST00000249064", "c.239_239+1insC", "p.Cys81LeufsTer10", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    #[case("22:28773779:G:GC", "ENST00000249064", "c.239+1_239+2insC", "p.?", vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant])]
    // minus strand: codon 282 CAT stays, codon 283 CCT (Pro) becomes TCC (Ser)
    #[case("22:28710005:C:CA", "ENST00000348295", "c.846dup", "p.Pro283SerfsTer7", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 283 becomes GCC (Ala)
    #[case("22:28710005:C:CC", "ENST00000348295", "c.846+1dup", "p.Pro283AlafsTer7", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 283 becomes CCC (Pro, as CCT), codon 284 TGC (Cys) becomes TTG (Leu)
    #[case("22:28710005:C:CG", "ENST00000348295", "c.846_846+1insC", "p.Cys284LeufsTer6", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    #[case("22:28710004:A:AG", "ENST00000348295", "c.846+1_846+2insC", "p.?", vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant])]
    // minus strand: codon 551 ATG (Met) becomes AAT (Asn)
    #[case("22:50246092:T:TT", "ENST00000216271", "c.1651dup", "p.Met551AsnfsTer23", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 551 becomes GAT (Asp)
    #[case("22:50246093:C:CC", "ENST00000216271", "c.1651-1dup", "p.Met551AspfsTer23", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    fn annotate_exon_boundary_dup_and_ins(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_hgvs_c: &str,
        #[case] expected_hgvs_p: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let ann = annotate_chr22_window(spdi, tx_id, false)?;
        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (Some(expected_hgvs_c), Some(expected_hgvs_p), &expected_csqs),
            "spdi = {spdi}"
        );

        Ok(())
    }

    /// A dup of the last CDS bases in the 3' UTR region gets 3' UTR terms. This also holds if
    /// the transcript is flagged as incomplete at its 3' end: hgvs-rs gives this region only
    /// behind a stop codon. An intronic dup stays intronic.
    #[rstest::rstest]
    #[case(
        "NM_000000.1:c.3857_3858dup",
        VariantLocation::ThreePrimeUtr,
        true,
        false,
        true
    )]
    #[case(
        "NM_000000.1:c.3857_3858dup",
        VariantLocation::ThreePrimeUtr,
        true,
        true,
        true
    )]
    #[case(
        "NM_000000.1:c.3858+1dup",
        VariantLocation::Intron,
        false,
        false,
        false
    )]
    fn analyze_cds_variant_dup_behind_stop(
        #[case] var_c: &str,
        #[case] region: VariantLocation,
        #[case] is_exonic: bool,
        #[case] incomplete_3p: bool,
        #[case] expected_utr: bool,
    ) -> Result<(), anyhow::Error> {
        let var_c = HgvsVariant::from_str(var_c)?;
        let csqs = ConsequencePredictor::analyze_cds_variant(
            &var_c,
            Some(region),
            is_exonic,
            false,
            incomplete_3p,
            Some(3858),
        );
        assert_eq!(
            csqs.intersects(
                Consequence::ThreePrimeUtrExonVariant | Consequence::ThreePrimeUtrIntronVariant
            ),
            expected_utr,
            "{:?}",
            csqs
        );

        Ok(())
    }

    /// An insertion or dup that can go before c.1 or after the last CDS base leaves the CDS
    /// intact, so it lies in a UTR. `c.1dup` inserts an A in front of the ATG, and `c.1_3dup` an
    /// ATG. The stop codon of NM_007294.4 ends at c.5592.
    #[rstest::rstest]
    #[case(
        "3:193311166:G:GC",
        "NM_130837.3",
        "c.-1_1insC",
        Consequence::FivePrimeUtrExonVariant
    )] // OPA1, forward
    #[case(
        "17:41258543:T:TA",
        "NM_007297.4",
        "c.-1_1insT",
        Consequence::FivePrimeUtrExonVariant
    )] // BRCA1, reverse
    #[case(
        "3:193311167:A:AA",
        "NM_130837.3",
        "c.1dup",
        Consequence::FivePrimeUtrExonVariant
    )]
    #[case(
        "17:41276112:A:AT",
        "NM_007294.4",
        "c.1dup",
        Consequence::FivePrimeUtrExonVariant
    )]
    #[case(
        "3:193311169:G:GATG",
        "NM_130837.3",
        "c.1_3dup",
        Consequence::FivePrimeUtrExonVariant
    )]
    #[case(
        "17:41276110:C:CCAT",
        "NM_007294.4",
        "c.1_3dup",
        Consequence::FivePrimeUtrExonVariant
    )]
    #[case(
        "17:41197695:T:TT",
        "NM_007294.4",
        "c.5592dup",
        Consequence::ThreePrimeUtrExonVariant
    )]
    fn annotate_ins_keeping_cds(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_hgvs_c: &str,
        #[case] expected_csq: Consequence,
    ) -> Result<(), anyhow::Error> {
        let spdi = spdi.split(':').collect::<Vec<_>>();

        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(provider, Default::default());

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].to_string(),
                position: spdi[1].parse()?,
                reference: spdi[2].to_string(),
                alternative: spdi[3].to_string(),
            })?
            .unwrap();

        let ann = res.iter().find(|ann| ann.feature_id == tx_id).unwrap();
        assert_eq!(ann.hgvs_c.as_deref(), Some(expected_hgvs_c));
        assert_eq!(ann.hgvs_p.as_deref(), Some("p.?"));
        assert_eq!(ann.consequences, vec![expected_csq]);

        Ok(())
    }

    /// With VEP terms, the expected terms are those of VEP 116.3 with `--shift_3prime 1`. VEP puts
    /// `c.1dup` between c.1 and c.2 and `c.1_3dup` between c.3 and c.4, so both change the CDS.
    /// `c.1_2insCCA` equals `c.-1_1insACC`.
    #[rstest::rstest]
    #[case("17:41276113:T:TG", "c.-1_1insC", vec![Consequence::FivePrimeUtrVariant])]
    #[case("17:41276112:A:AT", "c.1dup", vec![Consequence::FrameshiftVariant, Consequence::StartLost, Consequence::StartRetainedVariant])]
    #[case("17:41276110:C:CCAT", "c.1_3dup", vec![Consequence::InframeInsertion])]
    #[case("17:41276112:A:ATGG", "c.1_2insCCA", vec![Consequence::StartLost, Consequence::StartRetainedVariant])]
    fn annotate_ins_keeping_cds_vep(
        #[case] spdi: &str,
        #[case] expected_hgvs_c: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let spdi = spdi.split(':').collect::<Vec<_>>();

        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .vep_consequence_terms(true)
                .build()?,
        );

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].to_string(),
                position: spdi[1].parse()?,
                reference: spdi[2].to_string(),
                alternative: spdi[3].to_string(),
            })?
            .unwrap();

        let ann = res
            .iter()
            .find(|ann| ann.feature_id == "NM_007294.4")
            .unwrap();
        assert_eq!(ann.hgvs_c.as_deref(), Some(expected_hgvs_c));
        assert_eq!(ann.hgvs_p.as_deref(), Some("p.?"));
        assert_eq!(ann.consequences, expected_csqs);

        Ok(())
    }

    /// A deletion that can go before c.1 or after the last CDS base leaves the CDS intact. hgvs-rs
    /// then gives `p.?`, and the default terms are UTR terms. In NM_007294.4, c.-2, c.-1 and c.1
    /// are A bases. In NM_130837.3, c.3048 and c.*1 are A bases.
    ///
    /// VEP 116.3 with `--shift_3prime 1` keeps such a deletion at the start codon in the CDS and
    /// reports `start_retained_variant`. It shifts the one behind the stop codon into the 3' UTR.
    #[rstest::rstest]
    #[case("17:41276112:AT:A", "NM_007294.4", "c.1del", false, vec![Consequence::FivePrimeUtrExonVariant])]
    #[case("17:41276112:AT:A", "NM_007294.4", "c.1del", true, vec![Consequence::FrameshiftVariant, Consequence::StartRetainedVariant])]
    #[case("17:41276112:ATT:A", "NM_007294.4", "c.-1_1del", false, vec![Consequence::FivePrimeUtrExonVariant])]
    #[case("17:41276112:ATT:A", "NM_007294.4", "c.-1_1del", true, vec![Consequence::StartLost, Consequence::FivePrimeUtrVariant])]
    #[case("3:193409915:AA:A", "NM_130837.3", "c.3048del", false, vec![Consequence::ThreePrimeUtrExonVariant])]
    #[case("3:193409915:AA:A", "NM_130837.3", "c.3048del", true, vec![Consequence::ThreePrimeUtrVariant])]
    fn annotate_del_keeping_cds(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_hgvs_c: &str,
        #[case] vep_consequence_terms: bool,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let spdi = spdi.split(':').collect::<Vec<_>>();

        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .vep_consequence_terms(vep_consequence_terms)
                .build()?,
        );

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].to_string(),
                position: spdi[1].parse()?,
                reference: spdi[2].to_string(),
                alternative: spdi[3].to_string(),
            })?
            .unwrap();

        let ann = res.iter().find(|ann| ann.feature_id == tx_id).unwrap();
        assert_eq!(ann.hgvs_c.as_deref(), Some(expected_hgvs_c));
        assert_eq!(ann.hgvs_p.as_deref(), Some("p.?"));
        assert_eq!(ann.consequences, expected_csqs);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case("17:41197701:G:C", false, true)] // don't pick transcripts, report worst
    #[case("17:41197701:G:C", false, false)] // don't pick transcripts, report all
    #[case("17:41197701:G:C", true, true)] // pick transcripts, report worst
    #[case("17:41197701:G:C", true, false)] // pick transcripts, report all
    fn annotate_snv_brca1_transcript_picking_reporting(
        #[case] spdi: &str,
        #[case] pick_transcripts: bool,
        #[case] report_most_severe_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!(
            "{}-{}-{}",
            spdi.replace(':', "-"),
            pick_transcripts,
            !report_most_severe_consequence_only
        );

        let variant: VcfVariant = spdi.parse()?;

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        // After picking one transcript per gene, the most severe consequence has no choice left.
        let transcript_settings = match (pick_transcripts, report_most_severe_consequence_only) {
            (true, _) => pick_mane_by_gene(),
            (false, true) => TranscriptSettings::most_severe(PickUnit::Gene),
            (false, false) => TranscriptSettings::default(),
        };
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));

        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .transcript_settings(transcript_settings)
                .build()
                .unwrap(),
        );

        let res = predictor.predict(&variant)?.unwrap();

        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    // Test predictions on TTN where we have a ManeSelect and a ManePlusClinical
    // transcript.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case("2:179631246:G:A", false, true)] // don't pick transcripts, report worst
    #[case("2:179631246:G:A", false, false)] // don't pick transcripts, report all
    #[case("2:179631246:G:A", true, true)] // pick transcripts, report worst
    #[case("2:179631246:G:A", true, false)] // pick transcripts, report all
    fn annotate_snv_ttn_transcript_picking_reporting(
        #[case] spdi: &str,
        #[case] pick_transcripts: bool,
        #[case] report_most_severe_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!(
            "{}-{}-{}",
            spdi.replace(':', "-"),
            pick_transcripts,
            !report_most_severe_consequence_only
        );

        let variant: VcfVariant = spdi.parse()?;

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;

        // After picking one transcript per gene, the most severe consequence has no choice left.
        let transcript_settings = match (pick_transcripts, report_most_severe_consequence_only) {
            (true, _) => pick_mane_by_gene(),
            (false, true) => TranscriptSettings::most_severe(PickUnit::Gene),
            (false, false) => TranscriptSettings::default(),
        };
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));

        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .transcript_settings(transcript_settings)
                .build()
                .unwrap(),
        );

        let res = predictor.predict(&variant)?.unwrap();

        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// This is a set of variants where VEP and mehari to disagree,
    /// i.e. interesting/edge cases that are not as clear-cut as others.
    ///
    /// This test ensures we do not regress on these cases.
    #[tokio::test]
    async fn annotate_vep_disagreement_cases() -> Result<(), anyhow::Error> {
        let tx_path =
            "tests/data/annotate/db/grch38/GRCh38-ensembl.disagreement-subset.txs.bin.zst";

        let path_input_vcf = "tests/data/annotate/seqvars/vep.disagreement-cases.vcf";
        let output = NamedTempFile::new()?;
        let writer = open_variant_writer(output.as_ref()).await?;
        let mut seqvars_writer = crate::annotate::seqvars::SeqvarsVcfWriter::new(writer);
        run_with_writer(
            &mut seqvars_writer,
            &Args {
                threads: 1,
                reference: None,
                in_memory_reference: true,
                assembly: Some("grch38".into()),
                input: path_input_vcf.into(),
                output: output.as_ref().to_str().unwrap().into(),
                output_format: OutputFormat::Vcf,
                predictor_settings: PredictorSettings {
                    transcript_settings: TranscriptSettings {
                        restrict_transcripts: vec![TagCriterion::ManeSelect],
                        restrict_fallback: true,
                        ..TranscriptSettings::most_severe(PickUnit::Allele)
                    },
                    ..Default::default()
                },
                max_var_count: None,
                sources: crate::annotate::seqvars::Sources {
                    transcripts: Some(vec![tx_path.into()]),
                    frequencies: None,
                    clinvar: None,
                    ..Default::default()
                },
            },
        )
        .await?;
        seqvars_writer.shutdown().await?;

        let records_written = read_vcf(output).await?;

        let mut snapshot_data = BTreeMap::new();
        let header = noodles::vcf::io::reader::Builder::default()
            .build_from_path(path_input_vcf)?
            .read_header()?;

        for record in records_written {
            let key = format!(
                "{}:{}:{}:{}:{}",
                record.reference_sequence_name(),
                record
                    .variant_start()
                    .map_or_else(|| "0".into(), |s| s.to_string()),
                record
                    .variant_end(&header)
                    .map_or_else(|_| "0".into(), |s| s.to_string()),
                record.reference_bases(),
                record.alternate_bases().as_ref().join(",")
            );

            let ann_field = record.info().get("ANN").flatten().map(|v| match v {
                Value::Array(Array::String(inner)) => inner
                    .iter()
                    .map(|s| s.clone().unwrap_or_default())
                    .join("|"),
                _ => "".into(),
            });

            snapshot_data.insert(key, ann_field);
        }

        assert_yaml_snapshot!("vep_disagreement_cases_output", snapshot_data);

        Ok(())
    }

    async fn read_vcf(
        path: impl AsRef<Path>,
    ) -> Result<Vec<noodles::vcf::variant::RecordBuf>, anyhow::Error> {
        let mut output_reader = open_variant_reader(path.as_ref()).await?;
        let header = output_reader.read_header().await?;
        let mut record_iter = output_reader.records(&header).await;
        let mut records = Vec::new();
        while let Some(record) = record_iter.try_next().await? {
            records.push(record);
        }
        Ok(records)
    }

    #[derive(Debug, Deserialize)]
    struct Record {
        pub var: String,
        pub tx: String,
        pub csq: String,
    }

    // Compare to SnpEff annotated variants for OPA1, touching special cases.
    #[test]
    fn annotate_opa1_hand_picked_vars() -> Result<(), anyhow::Error> {
        annotate_opa1_vars("tests/data/annotate/seqvars/opa1.hand_picked.tsv", false)
    }

    // Compare to SnpEff annotated ClinVar variants for OPA1 (slow).
    #[test]
    fn annotate_opa1_clinvar_vars_snpeff() -> Result<(), anyhow::Error> {
        annotate_opa1_vars(
            "tests/data/annotate/seqvars/clinvar.excerpt.snpeff.opa1.tsv",
            false,
        )
    }

    // Compare to SnpEff annotated ClinVar variants for OPA1 (slow).
    #[test]
    fn annotate_opa1_clinvar_vars_vep() -> Result<(), anyhow::Error> {
        annotate_opa1_vars(
            "tests/data/annotate/seqvars/clinvar.excerpt.vep.opa1.tsv",
            false,
        )
    }

    fn annotate_opa1_vars(
        path_tsv: &str,
        report_most_severe_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        let txs = vec![
            String::from("NM_001354663.2"),
            String::from("NM_001354664.2"),
            String::from("NM_015560.3"),
            String::from("NM_130831.3"),
            String::from("NM_130832.3"),
            String::from("NM_130837.3"),
        ];

        annotate_vars(path_tsv, &txs, report_most_severe_consequence_only, false)
    }

    // Compare to SnpEff annotated variants for BRCA1, touching special cases.
    #[test]
    fn annotate_brca1_hand_picked_vars() -> Result<(), anyhow::Error> {
        annotate_brca1_vars("tests/data/annotate/seqvars/brca1.hand_picked.tsv", false)
    }

    // Compare to SnpEff annotated ClinVar variants for BRCA1 (slow).
    #[test]
    fn annotate_brca1_clinvar_vars_snpeff() -> Result<(), anyhow::Error> {
        annotate_brca1_vars(
            "tests/data/annotate/seqvars/clinvar.excerpt.snpeff.brca1.tsv",
            false,
        )
    }

    // Compare to SnpEff annotated ClinVar variants for BRCA1 (slow).
    #[test]
    fn annotate_brca1_clinvar_vars_vep() -> Result<(), anyhow::Error> {
        annotate_brca1_vars(
            "tests/data/annotate/seqvars/clinvar.excerpt.vep.brca1.tsv",
            false,
        )
    }

    fn annotate_brca1_vars(
        path_tsv: &str,
        report_most_severe_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        let txs = vec![
            String::from("NM_007294.4"),
            String::from("NM_007297.4"),
            String::from("NM_007298.3"),
            String::from("NM_007299.4"),
            String::from("NM_007300.4"),
        ];

        annotate_vars(path_tsv, &txs, report_most_severe_consequence_only, false)
    }

    fn annotate_vars(
        path_tsv: &str,
        txs: &[String],
        report_most_severe_consequence_only: bool,
        with_reference: bool,
    ) -> Result<(), anyhow::Error> {
        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let reference_path = "resources/GCF_000001405.25_GRCh37.p13_genomic.fna";
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            with_reference.then_some(reference_path),
            true,
        ));

        let transcript_settings = if report_most_severe_consequence_only {
            TranscriptSettings::most_severe(PickUnit::Gene)
        } else {
            TranscriptSettings::default()
        };

        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .transcript_settings(transcript_settings)
                .build()
                .unwrap(),
        );

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_reader(File::open(path_tsv).map(BufReader::new)?);

        // Read record with variant, transcript, and predicted consequences.
        //
        // We only have a limited set of transcripts for BRCA1 and this may not be the
        // same as the one in the file.  We thus limit ourselves to the transcripts
        // in `txs`.
        //
        // Also, the predicted ontology terms may not be the same as the ones in the
        // file.  We thus only check that we have a match in at least one term
        // of highest impact.
        let mut lineno = 0;
        for record in reader.deserialize() {
            lineno += 1;

            let record: Record = record?;

            // Because for this variant our highest impact is "HIGH" and vep's is not (and because we filter for the highest impact), skip it.
            // We predict FrameshiftVariant, FrameshiftTruncation, SpliceRegionVariant and ThreePrimeUtrExonVariant,
            // while vep calls "splice_region_variant", "coding_sequence_variant", "3_prime_UTR_variant".
            if record.var
                == "17-41196310-GGTGGAAGTGTTTGCTACCAAGTTTATTTGCAGTGTTAACAGCACAACATTTACAAAACGTATTTTGTACAATCAAGTCTTCACTGCCCTTGCACACTGGGGGGGCTAGGGAAGACCTAGTCCTTCCAACAGCTATAAACAGTCCTGGATAATGGGTTTATGAAAAACACTTTTTCTTCCTTCAGCAAGCAAAATTATTTATGAAGCTGTATGGTTTCAGCAACAGGGAGCAAAGGAAAAAAATCACCTCAAAGAAAGCAACAGCTTCCTTCCTGGTGGGATCTGTCATTTTATAGATATGAAATATTCATGCCAGAGGTCTTATATTTTAAGAGGAATGGATTATATACCAGAGCTACAACAATAAACATTTTACTTATTACTAATGAGGAATTAGAAGACTGTCTTTGGAAACCGGTTCTTGAAAATCTTCTGCTGTTTTAGAACACATTCTTTAGAAATCTAGCAAATATATCTCAGACTTTTAGAAATCTCTTCTAGTTTCATTTTCCTTTTTTTTTTTTTTTTTTTGAGCCACAGTCTCACTGTCACCCAGGCTGGAGTGCCGTGGTATGATCTTGGCTCACTGCAACCTCCACCTCCCGGGCTGAAGTGATTCTCCTGCCTTAGCCACCTGAGTAGCTGGGATTACAGGTGTCCACCACCATGACCGGCTAATTTCTGTATTTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTTTCGAACTCCTGACCTCCAGTGATCTGCCCACCTTGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCATGCCCAGGTTTCAAGTTTCCTTTTCATTTCTAATACCTGCCTCAGAATTTCCTCCCCAATGTTCCACTCCAACATTTGAGAACTGCCCAAGGACTATTCTGACTTTAAGTCACATAATCGATCCCAAGCACTCTCCTTCCATTGAAGGGTCTGACTCTCTGCCTTTGTGAACACAGGGTTTTAGAGAAGTAAACTTAGGGAAACCAGCTATTCTCTTGAGGCCAAGCCACTCTGTGCTTCCAGCCCTAAGCCAACAACAGCCTGAATAGAAAGAATAGGGCTGATAAATAATGAATCAGCATCTTGCTCAATTGGTGGCGTTTAAATGGTTTTAAAATCTTCTCAGGTGAAAAATTACCATAATTTTGTGCTCATGGCAGATTTCCAAGGGAGACTTCAAGCAGAAAATCTTTAAGGGACCCTTGCATAGCCAGAAGTCCTTTTCAGGCTGATGTACATAAAATATTTAGTAGCCAGGACAGTAGAAGGACTGAAGAGTGAGAGGAGCTCCCAGGGCCTGGAAAGGCCACTTTGTAAGCTCATTCTTGGGGTCCTGTGGCTCTGTACCTGTGGCTGGCTGCAGTCAGTAGTGGCTGTGGGGGATCTGGGGTATCAGGTAGGTGTCCAGCTCCTGGCACTGGTAGAGTGCTACACTGTCCAACACCCACTCTCGGGTCACCACAGGTGCCTCACACATCTGCCCAATT-G"
            {
                continue;
            }

            if txs.contains(&record.tx) {
                // "Parse" out the variant.
                let arr = record.var.split('-').collect::<Vec<_>>();
                // Predict consequences for the variant using Mehari.
                let anns = predictor
                    .predict(&VcfVariant {
                        chromosome: arr[0].to_string(),
                        position: arr[1].parse::<i32>()?,
                        reference: arr[2].to_string(),
                        alternative: arr[3].to_string(),
                    })?
                    .unwrap();
                // Now, for the overlapping transcripts, check that we have a match with the
                // consequences in the highest impact category.
                for ann in anns.iter().filter(|ann| ann.feature_id == record.tx) {
                    // We perform a comparison based on strings because we may not be able to parse out
                    // all consequences from the other tool.
                    let record_csqs = record.csq.split('&').collect::<Vec<_>>();

                    let highest_impact = ann.consequences.first().unwrap().impact();
                    let expected_one_of = ann
                        .consequences
                        .iter()
                        .filter(|csq| csq.impact() == highest_impact)
                        .map(|csq| csq.to_string())
                        .collect::<Vec<_>>();
                    let mut expected_one_of = expected_one_of
                        .iter()
                        .map(|s| s.as_str())
                        .collect::<Vec<_>>();

                    // Map effects a bit for VEP.
                    if path_tsv.contains(".vep")
                        && (expected_one_of.contains(&"disruptive_inframe_deletion")
                            || expected_one_of.contains(&"conservative_inframe_deletion"))
                    {
                        expected_one_of.push("inframe_deletion");
                    }

                    let found_one = [
                        // Try to find a direct match.
                        record_csqs.iter().any(|csq| expected_one_of.contains(csq)),
                        // vep sometimes only reports a coding_sequence_variant, so we accept anything
                        path_tsv.contains(".vep")
                            && (record_csqs == ["coding_sequence_variant"]
                                && !expected_one_of.is_empty()),
                        // It is common that the other tool predicts a frameshift variant while the actual prediction
                        // is stop_gained or stop_lost.  We thus also check for this case and allow it.
                        (record_csqs.contains(&"frameshift_variant")
                            || record_csqs.contains(&"frameshift_truncation")
                            || record_csqs.contains(&"frameshift_elongation"))
                            && (expected_one_of.contains(&"stop_gained"))
                            || expected_one_of.contains(&"stop_lost"),
                        // … or vice-versa
                        (expected_one_of.contains(&"frameshift_variant")
                            || expected_one_of.contains(&"frameshift_truncation")
                            || expected_one_of.contains(&"frameshift_elongation"))
                            && (record_csqs.contains(&"stop_gained"))
                            || record_csqs.contains(&"stop_lost"),
                        // VEP does not differentiate between disruptive and conservative inframe deletions and insertions.
                        (record_csqs.contains(&"inframe_deletion")
                            && (expected_one_of.contains(&"disruptive_inframe_deletion")
                                || expected_one_of.contains(&"conservative_inframe_deletion")))
                            || (record_csqs.contains(&"inframe_insertion")
                                && (expected_one_of.contains(&"disruptive_inframe_insertion")
                                    || expected_one_of
                                        .contains(&"conservative_inframe_insertion"))),
                        // delins on protein level are sometimes erroneously reported as inframe_insertion/deletion instead
                        (expected_one_of.contains(&"protein_altering_variant")
                            && ([
                                "disruptive_inframe_deletion",
                                "disruptive_inframe_insertion",
                                "conservative_inframe_deletion",
                                "conservative_inframe_insertion",
                            ]
                            .iter()
                            .any(|c| record_csqs.contains(c)))),
                        // NB: We cannot predict 5_prime_UTR_premature_start_codon_gain_variant yet. For now, we
                        // also accept 5_prime_UTR_variant.
                        ((expected_one_of.contains(&"5_prime_UTR_exon_variant")
                            || expected_one_of.contains(&"5_prime_UTR_intron_variant"))
                            && (record_csqs
                                .contains(&"5_prime_UTR_premature_start_codon_gain_variant"))),
                        // We accept 5_prime_UTR_exon_variant and 5_prime_UTR_intron_variant if the
                        // other tool predicts upstream_gene_variant.
                        ((expected_one_of.contains(&"5_prime_UTR_exon_variant")
                            || expected_one_of.contains(&"5_prime_UTR_intron_variant"))
                            && (record_csqs.contains(&"upstream_gene_variant"))),
                        // A coding_transcript_intron_variant is a more specific intron_variant
                        (expected_one_of.contains(&"coding_transcript_intron_variant")
                            && (record_csqs.contains(&"intron_variant"))),
                        // VEP predicts `splice_donor_5th_base_variant` rather than `splice_region_variant`.
                        // Same for `splice_donor_region_variant`.
                        (expected_one_of.contains(&"splice_region_variant")
                            && (record_csqs.contains(&"splice_donor_5th_base_variant")
                                || record_csqs.contains(&"splice_donor_region_variant"))),
                        // In the case of insertions at the end of an exon, VEP predicts `splice_region_variant`
                        // while we predict `splice_donor_variant`, same for start.
                        (expected_one_of.contains(&"splice_donor_variant")
                            || expected_one_of.contains(&"splice_acceptor_variant"))
                            && (record_csqs.contains(&"splice_region_variant")),
                        // VEP sometimes mispredicts disruptive inframe deletion as missense...
                        // cf. https://github.com/Ensembl/ensembl-vep/issues/1388
                        expected_one_of.contains(&"disruptive_inframe_deletion")
                            && (record_csqs.contains(&"missense_variant")),
                        // VEP does not provide `exon_loss_variant`, so we also accept `inframe_deletion` and
                        // `splice_region_variant` (BRA1 test case).
                        expected_one_of.contains(&"exon_loss_variant")
                            && (record_csqs.contains(&"inframe_deletion")
                                || record_csqs.contains(&"splice_region_variant")),
                        // On BRCA1, there is a case where VEP predicts `protein_altering_variant` rather than
                        // `disruptive_inframe_deletion`.  We accept this as well.
                        (expected_one_of.contains(&"disruptive_inframe_deletion")
                            || expected_one_of.contains(&"inframe_indel"))
                            && (record_csqs.contains(&"protein_altering_variant")),
                        // We have specialized {5,3}_prime_UTR_{exon,intron}_variant handling, while
                        // vep and snpEff do not
                        record_csqs.contains(&"5_prime_UTR_variant")
                            && (expected_one_of.contains(&"5_prime_UTR_exon_variant")
                                || expected_one_of.contains(&"5_prime_UTR_intron_variant")),
                        record_csqs.contains(&"3_prime_UTR_variant")
                            && (expected_one_of.contains(&"3_prime_UTR_exon_variant")
                                || expected_one_of.contains(&"3_prime_UTR_intron_variant")),
                        // an inframe_indel can be a missense_variant if it is an MNV (which we do not explicitly check here)
                        expected_one_of.contains(&"inframe_indel")
                            && (record_csqs.contains(&"missense_variant")),
                        // inframe_indel also is a superclass of *_inframe_{deletion, insertion}
                        expected_one_of.contains(&"inframe_indel")
                            && [
                                "disruptive_inframe_deletion",
                                "conservative_inframe_deletion",
                                "inframe_deletion",
                                "disruptive_inframe_insertion",
                                "conservative_inframe_insertion",
                                "inframe_insertion",
                            ]
                            .iter()
                            .any(|c| record_csqs.contains(c)),
                        // SnpEff has a different interpretation of disruptive/conservative inframe deletions.
                        // We thus allow both.
                        expected_one_of.contains(&"disruptive_inframe_deletion")
                            && (record_csqs.contains(&"conservative_inframe_deletion"))
                            || expected_one_of.contains(&"disruptive_inframe_insertion")
                                && (record_csqs.contains(&"conservative_inframe_insertion")),
                        // SnpEff may not predict `splice_region_variant` for 5' UTR correctly, so we
                        // allow this.
                        (expected_one_of.contains(&"splice_region_variant")
                            || expected_one_of.contains(&"exonic_splice_region_variant"))
                            && (record_csqs.contains(&"5_prime_UTR_variant")),
                        // SnpEff does not predict `splice_polypyrimidine_tract_variant`
                        expected_one_of.contains(&"splice_polypyrimidine_tract_variant")
                            && (record_csqs.contains(&"splice_region_variant")
                                || record_csqs.contains(&"intron_variant")),
                        // For `GRCh37:3:193366573:A:ATATTGCCTAGAATGAACT`, SnpEff predicts
                        // `stop_gained` while this rather is a intron variant.  We skip this variant.
                        record_csqs.contains(&"stop_gained")
                            && record.var == "3-193366573-A-ATATTGCCTAGAATGAACT",
                        // For `GRCh37:3:193409913:ATAAAT:A`, there appears to be a model error
                        // in SnpEff as it predicts `exon_loss`.  We skip this variant.
                        record_csqs.contains(&"exon_loss_variant")
                            && record.var == "3-193409913-ATAAAT-A",
                        // SnpEff may predict `pMet1.?` as `initiator_codon_variant` rather than `start_lost`.
                        expected_one_of.contains(&"start_lost")
                            && (record_csqs.contains(&"initiator_codon_variant")),
                        // SnpEff predicts `c.-1_1ins` as `start_retained` while VEP and we predict a
                        // 5' UTR variant.
                        expected_one_of.contains(&"5_prime_UTR_exon_variant")
                            && (record_csqs.contains(&"start_retained_variant")),
                        // SnpEff calls this insertion at c.5193+2_5193+3insT a splice donor variant
                        // even though the third intronic base is affected, not the first or second
                        record_csqs.contains(&"splice_donor_variant")
                            && expected_one_of.contains(&"splice_region_variant")
                            && [
                                "17-41215347-T-TA",
                                "17-41215888-T-TA",
                                "17-41242958-T-TA",
                                "17-41256882-T-TA",
                                "17-41276031-T-TA",
                                "17-41277285-T-TA",
                            ]
                            .contains(&record.var.as_str()),
                        // SnpEff calls these insertions at an exon edge a splice donor or acceptor
                        // variant. The splice site stays intact, and the inserted bases join the
                        // exon: the CDS of OPA1 and the 5' UTR of BRCA1.
                        record_csqs.contains(&"splice_donor_variant")
                            && expected_one_of.contains(&"frameshift_variant")
                            && record.var == "3-193363589-A-AG",
                        record_csqs.contains(&"splice_acceptor_variant")
                            && expected_one_of.contains(&"exonic_splice_region_variant")
                            && record.var == "17-41276132-A-ACT",
                        // we call exonic_splice_region_variant, while the others only call splice_region_variant
                        record_csqs.contains(&"splice_region_variant")
                            && expected_one_of.contains(&"exonic_splice_region_variant"),
                    ]
                    .iter()
                    .any(|b| *b);

                    assert!(
                        found_one,
                        "line no. {}, variant: {}, tx: {}, hgvs_c: {:?}, hgvs_p: {:?}, \
                        their_csqs: {:?}, expected_one_of: {:?}, our_csqs: {:?}",
                        lineno,
                        record.var,
                        record.tx,
                        ann.hgvs_c.as_ref(),
                        ann.hgvs_p.as_ref(),
                        record_csqs,
                        expected_one_of,
                        ann.consequences,
                    );
                }
            }
        }

        Ok(())
    }

    /// A predictor for `tx`, a transcript on chr1 (GRCh38) with the sequence `seq`. It reports
    /// the transcript sequences.
    fn predictor_for(tx: Transcript, seq: String) -> ConsequencePredictor {
        use crate::pbs::txs::{GeneToTxId, SequenceDb, SourceVersion, TranscriptDb, TxSeqDatabase};

        // hgvs-rs caches the reference protein by the database version, so each database
        // needs its own version.
        let version = format!("{tx:?} {seq}");
        let tx_seq_db = TxSeqDatabase {
            tx_db: Some(TranscriptDb {
                gene_to_tx: vec![GeneToTxId {
                    gene_id: tx.gene_id.clone(),
                    tx_ids: vec![tx.id.clone()],
                    filtered: Some(false),
                    filter_reason: None,
                }],
                transcripts: vec![tx.clone()],
            }),
            seq_db: Some(SequenceDb {
                aliases: vec![tx.id],
                aliases_idx: vec![0],
                seqs: vec![seq],
            }),
            version: Some(version),
            source_version: vec![SourceVersion {
                assembly: "GRCh38".into(),
                ..Default::default()
            }],
        };
        let provider = Arc::new(MehariProvider::new(tx_seq_db, None::<PathBuf>, true));
        let config = ConfigBuilder::default()
            .plugins(vec![
                Arc::new(SequencePlugin::TxRef),
                Arc::new(SequencePlugin::TxAlt),
            ])
            .build()
            .unwrap();
        ConsequencePredictor::new(provider, config)
    }

    /// A coding transcript `NM_000001.1` with one exon at chr1:1001 on the plus strand. Its
    /// sequence is `utr5`, `cds` and `utr3`.
    fn one_exon_tx(utr5: &str, cds: &str, utr3: &str) -> (Transcript, String) {
        use crate::pbs::txs::ExonAlignment;

        let seq = format!("{utr5}{cds}{utr3}");
        let len = i32::try_from(seq.len()).unwrap();
        let start = i32::try_from(utr5.len()).unwrap();
        let stop = start + i32::try_from(cds.len()).unwrap();
        let tx = Transcript {
            id: "NM_000001.1".into(),
            gene_symbol: "GENE1".into(),
            gene_id: "1".into(),
            biotype: TranscriptBiotype::Coding.into(),
            protein: Some("NP_000001.1".into()),
            start_codon: Some(start),
            stop_codon: Some(stop),
            genome_alignments: vec![GenomeAlignment {
                genome_build: "grch38".into(),
                contig: "NC_000001.11".into(),
                cds_start: Some(1000 + start),
                cds_end: Some(1000 + stop),
                strand: Strand::Plus.into(),
                exons: vec![ExonAlignment {
                    alt_start_i: 1000,
                    alt_end_i: 1000 + len,
                    ord: 0,
                    alt_cds_start_i: Some(1),
                    alt_cds_end_i: Some(len),
                    cigar: format!("{len}M"),
                }],
                ..Default::default()
            }],
            filtered: Some(false),
            ..Default::default()
        };
        (tx, seq)
    }

    /// A predictor for one selenoprotein transcript on the plus strand of chr1 (GRCh38).
    ///
    /// Its one exon spans chr1:1001-1073: 20 bases of 5' UTR, the CDS
    /// `ATG GCC AAG CTG TGG GAA CCA TGA CGC GTT TAA` (`MAKLWEPURV*`) at chr1:1021-1053, and 20
    /// bases of 3' UTR. The eighth codon is the Sec codon.
    fn selenoprotein_predictor(tagged: bool, positions: Vec<u32>) -> ConsequencePredictor {
        use crate::pbs::txs::TranslationException;

        let utr = "CAGTCAGTCAGTCAGTCAGT";
        let (tx, seq) = one_exon_tx(utr, "ATGGCCAAGCTGTGGGAACCATGACGCGTTTAA", utr);
        let tx = Transcript {
            gene_symbol: "SELENOX".into(),
            tags: tagged
                .then(|| TranscriptTag::Selenoprotein.into())
                .into_iter()
                .collect(),
            translation_exceptions: positions
                .into_iter()
                .map(|position| TranslationException {
                    position,
                    amino_acid: "U".into(),
                })
                .collect(),
            ..tx
        };
        predictor_for(tx, seq)
    }

    /// With Sec positions, UGA reads as selenocysteine only at them. Without, a transcript
    /// tagged as selenoprotein reads every UGA as selenocysteine.
    #[rstest::rstest]
    #[case::new_uga_is_a_stop(false, vec![8], "1035:G:A", "p.Trp5Ter", Consequence::StopGained)]
    #[case::new_uga_is_a_stop_despite_tag(
        true,
        vec![8],
        "1035:G:A",
        "p.Trp5Ter",
        Consequence::StopGained
    )]
    #[case::sec_codon_reads_sec(false, vec![8], "1044:A:G", "p.Sec8Trp", Consequence::MissenseVariant)]
    #[case::codon_after_sec_codon(
        false,
        vec![8],
        "1046:G:A",
        "p.Arg9His",
        Consequence::MissenseVariant
    )]
    #[case::sec_codon_moves_with_deletion(
        false,
        vec![8],
        "1035:GGAA:G",
        "p.Glu6del",
        Consequence::ConservativeInframeDeletion
    )]
    #[case::tag_without_positions(true, vec![], "1035:G:A", "p.Trp5Sec", Consequence::MissenseVariant)]
    fn selenocysteine_only_at_its_positions(
        #[case] tagged: bool,
        #[case] positions: Vec<u32>,
        #[case] var: &str,
        #[case] hgvs_p: &str,
        #[case] consequence: Consequence,
    ) -> Result<(), anyhow::Error> {
        let var = var.split(':').collect::<Vec<_>>();
        let anns = selenoprotein_predictor(tagged, positions)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: var[0].parse()?,
                reference: var[1].into(),
                alternative: var[2].into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(ann.hgvs_p.as_deref(), Some(hgvs_p));
        assert!(
            ann.consequences.contains(&consequence),
            "{:?}",
            ann.consequences
        );

        Ok(())
    }

    /// A selenoprotein without Sec positions reads UGA as selenocysteine inside its CDS, but a
    /// UGA at the CDS end is its stop codon, as in SEPHS2.
    ///
    /// The transcript is the one of `selenoprotein_predictor`, but its CDS
    /// `ATG GCC AAG CTG TGG GAA CCA TGA CGC GTT TGA` (`MAKLWEPURV*`) ends with UGA.
    #[rstest::rstest]
    #[case::stop_codon("1051:T:C", "p.Ter11ArgextTer?", Consequence::StopLost)]
    #[case::sec_codon("1044:A:G", "p.Sec8Trp", Consequence::MissenseVariant)]
    fn last_uga_is_the_stop_codon(
        #[case] var: &str,
        #[case] hgvs_p: &str,
        #[case] consequence: Consequence,
    ) -> Result<(), anyhow::Error> {
        let utr = "CAGTCAGTCAGTCAGTCAGT";
        let (tx, seq) = one_exon_tx(utr, "ATGGCCAAGCTGTGGGAACCATGACGCGTTTGA", utr);
        let tx = Transcript {
            tags: vec![TranscriptTag::Selenoprotein.into()],
            ..tx
        };
        let var = var.split(':').collect::<Vec<_>>();
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: var[0].parse()?,
                reference: var[1].into(),
                alternative: var[2].into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(ann.hgvs_p.as_deref(), Some(hgvs_p));
        assert!(
            ann.consequences.contains(&consequence),
            "{:?}",
            ann.consequences
        );

        Ok(())
    }

    /// A mitochondrial transcript reads UGA as Trp, as the vertebrate mitochondrial code does.
    /// `c.15G>A` turns TGG into TGA.
    #[test]
    fn mitochondrial_code_reads_uga_as_trp() -> Result<(), anyhow::Error> {
        let (mut tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CGCATCA");
        tx.genome_alignments[0].contig = "NC_012920.1".into();
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "MT".into(),
                position: 1035,
                reference: "G".into(),
                alternative: "A".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (Some("c.15G>A"), &vec![Consequence::SynonymousVariant])
        );
        Ok(())
    }

    /// Without a stop codon, a change of the last amino acid is an ordinary change of it, not an
    /// extension. A frameshift that leaves no complete codon there gives `p.?`.
    ///
    /// The transcript has 20 bases of 5' UTR and the CDS `ATG GCC AAG CTG TGG GAA` (`MAKLWE`)
    /// at chr1:1021-1038, which ends at the transcript end. `db create` flags it with
    /// `MissingStopCodon`.
    #[rstest::rstest]
    #[case::synonymous("1038:A:G", "p.Glu6=", &[Consequence::SynonymousVariant])]
    #[case::missense("1038:A:C", "p.Glu6Asp", &[Consequence::MissenseVariant])]
    // c.16_18del
    #[case::deletion(
        "1035:GGAA:G",
        "p.Glu6del",
        &[Consequence::ConservativeInframeDeletion]
    )]
    // c.17_18insC
    #[case::frameshift("1037:A:AC", "p.Glu6AspfsTer?", &[Consequence::FrameshiftVariant])]
    // c.18del
    #[case::frameshift_in_last_codon("1036:GA:G", "p.?", &[Consequence::FrameshiftVariant])]
    // c.16del
    #[case::frameshift_at_last_codon("1034:GG:G", "p.?", &[Consequence::FrameshiftVariant])]
    fn last_amino_acid_without_stop_codon(
        #[case] var: &str,
        #[case] hgvs_p: &str,
        #[case] expected: &[Consequence],
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGGAA", "");
        let tx = Transcript {
            filter_reason: Some(BitFlags::from(Reason::MissingStopCodon).bits()),
            ..tx
        };
        let var = var.split(':').collect::<Vec<_>>();
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: var[0].parse()?,
                reference: var[1].into(),
                alternative: var[2].into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(ann.hgvs_p.as_deref(), Some(hgvs_p), "{ann:?}");
        for consequence in [
            Consequence::SynonymousVariant,
            Consequence::MissenseVariant,
            Consequence::ConservativeInframeDeletion,
            Consequence::FrameshiftVariant,
            Consequence::FeatureElongation,
        ] {
            assert_eq!(
                ann.consequences.contains(&consequence),
                expected.contains(&consequence),
                "{consequence:?}: {ann:?}"
            );
        }
        Ok(())
    }

    /// A CDS whose end the annotation marks as incomplete keeps its partial last codon, and
    /// `db create` flags it with `MissingStopCodon`. The protein ends with the last complete
    /// codon. The amino acids after it are unknown, so a variant after it gives `p.?`.
    ///
    /// The transcript has 20 bases of 5' UTR and the CDS `ATG GCC AAG CTG TGG GAA CC`
    /// (`MAKLWE` and 2 bases) at chr1:1021-1040, which ends at the transcript end.
    #[rstest::rstest]
    #[case::stop_gained("1035:G:A", "p.Trp5Ter", &[Consequence::StopGained])]
    #[case::last_complete_codon("1038:A:G", "p.Glu6=", &[Consequence::SynonymousVariant])]
    #[case::partial_codon(
        "1039:C:T",
        "p.?",
        &[Consequence::IncompleteTerminalCodonVariant]
    )]
    #[case::frameshift_into_partial_codon(
        "1039:C:CGTGGGAAC",
        "p.?",
        &[Consequence::FrameshiftVariant, Consequence::IncompleteTerminalCodonVariant]
    )]
    // c.13_14insCCTA: the new stop takes the place of the last complete codon
    #[case::frameshift_truncation(
        "1033:T:TCCTA",
        "p.Trp5SerfsTer2",
        &[Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation]
    )]
    // c.16_18dup
    #[case::dup_of_last_complete_codon(
        "1035:G:GGAA",
        "p.?",
        &[Consequence::ConservativeInframeInsertion]
    )]
    fn incomplete_cds_end_has_no_stop_codon(
        #[case] var: &str,
        #[case] hgvs_p: &str,
        #[case] expected: &[Consequence],
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGGAACC", "");
        let tx = Transcript {
            filter_reason: Some(BitFlags::from(Reason::MissingStopCodon).bits()),
            ..tx
        };
        let var = var.split(':').collect::<Vec<_>>();
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: var[0].parse()?,
                reference: var[1].into(),
                alternative: var[2].into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(ann.hgvs_p.as_deref(), Some(hgvs_p), "{ann:?}");
        for consequence in expected {
            assert!(ann.consequences.contains(consequence), "{ann:?}");
        }
        for consequence in [
            Consequence::StopLost,
            Consequence::StopRetainedVariant,
            Consequence::FeatureElongation,
            Consequence::IncompleteTerminalCodonVariant,
        ] {
            if !expected.contains(&consequence) {
                assert!(!ann.consequences.contains(&consequence), "{ann:?}");
            }
        }
        Ok(())
    }

    /// `db create` completes a stop codon at the transcript end with `A` bases. They follow the
    /// last exon, and no total or sequence counts them.
    ///
    /// The transcript has 20 bases of 5' UTR and the CDS `ATG GCC AAG CTG TGG GAA T`
    /// (`MAKLWE` and the stop codon `T`) at chr1:1021-1039, which ends at the transcript end.
    /// `db create` completes the stop codon to `TAA`.
    #[rstest::rstest]
    #[case::missense("1033:T:C", "p.Trp5Arg", Consequence::MissenseVariant)]
    #[case::stop_lost("1039:T:C", "p.Ter7GlnextTer?", Consequence::StopLost)]
    fn totals_without_poly_a_padding(
        #[case] var: &str,
        #[case] hgvs_p: &str,
        #[case] consequence: Consequence,
    ) -> Result<(), anyhow::Error> {
        let (utr5, cds) = ("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGGAAT");
        let (tx, seq) = one_exon_tx(utr5, cds, "");
        let tx = Transcript {
            stop_codon: tx.stop_codon.map(|stop| stop + 2),
            ..tx
        };
        let var = var.split(':').collect::<Vec<_>>();
        let position: usize = var[0].parse()?;
        let anns = predictor_for(tx, format!("{seq}AA"))
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: position.try_into()?,
                reference: var[1].into(),
                alternative: var[2].into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(ann.hgvs_p.as_deref(), Some(hgvs_p), "{ann:?}");
        assert!(ann.consequences.contains(&consequence), "{ann:?}");
        let totals = [&ann.cdna_pos, &ann.cds_pos].map(|pos| pos.as_ref().and_then(|p| p.total));
        assert_eq!(totals, [Some(39), Some(19)], "{ann:?}");
        assert_eq!(custom_field(ann, ANN_TX_SEQ_REF), Some(seq.as_str()));
        let mut seq_alt = seq.clone();
        seq_alt.replace_range(position - 1001..position - 1000, var[2]);
        assert_eq!(custom_field(ann, ANN_TX_SEQ_ALT), Some(seq_alt.as_str()));
        Ok(())
    }

    /// Deleting the last base of a stop codon at the transcript end loses the stop codon. No
    /// appended `A` bases complete a new one.
    ///
    /// The transcripts come from `db create` with the GENCODE 42 GFF3 records of FAM246B. Its
    /// transcript `ENST00000652395.1` has one exon, which ends with the stop codon `TGA` at
    /// chr22:18634680-18634682.
    #[rstest::rstest]
    #[case(false, vec![Consequence::StopLost, Consequence::FeatureElongation])]
    // VEP 116.3 reports `stop_lost` only
    #[case(true, vec![Consequence::FrameshiftVariant, Consequence::StopLost, Consequence::FeatureElongation])]
    fn stop_codon_at_transcript_end_without_padding(
        #[case] vep_consequence_terms: bool,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let tx_db = load_tx_db(
            "tests/data/annotate/db/grch38/GRCh38-gencode42.stop-at-tx-end.txs.bin.zst",
        )?;
        // The test needs no reference: a `G` follows the deleted `A`, so the deletion cannot
        // shift.
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .vep_consequence_terms(vep_consequence_terms)
                .build()?,
        );

        let anns = predictor
            .predict(&VcfVariant {
                chromosome: "22".into(),
                position: 18634681,
                reference: "GA".into(),
                alternative: "G".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "ENST00000652395.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for ENST00000652395.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (Some("c.699del"), Some("p.Ter233XaaextTer?"), &expected_csqs)
        );
        Ok(())
    }

    #[test]
    fn test_predict_multiple_bounds_checking() -> Result<(), anyhow::Error> {
        let alt_seq = String::from("ATGCGTACGTAGCTAGCT");
        let n_min = 5;
        let n_max = 7;

        let total_delta_negative = -10;
        let new_length_negative = (n_max - n_min + 1) + total_delta_negative;
        assert!(
            new_length_negative < 0,
            "Length should be negative and caught by the guard"
        );

        let total_delta_past_end = 100;
        let new_length = (n_max - n_min + 1) + total_delta_past_end;
        let start_idx = (n_min - 1) as usize;
        let end_idx = start_idx + new_length as usize;

        assert!(
            end_idx > alt_seq.len(),
            "End index should exceed alt_seq length and be caught by the guard"
        );

        let total_delta_valid = 2; // e.g., del 1 base, ins 3 bases
        let new_length_valid = (n_max - n_min + 1) + total_delta_valid;
        let end_idx_valid = start_idx + new_length_valid as usize;

        assert!(new_length_valid >= 0);
        assert!(end_idx_valid <= alt_seq.len());
        let _new_substring = &alt_seq[start_idx..end_idx_valid];

        Ok(())
    }

    fn brca1_sequence_predictor() -> Result<ConsequencePredictor, anyhow::Error> {
        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        Ok(ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .plugins(vec![
                    Arc::new(SequencePlugin::TxRef),
                    Arc::new(SequencePlugin::TxAlt),
                    Arc::new(SequencePlugin::AaRef),
                    Arc::new(SequencePlugin::AaAlt),
                ])
                .build()?,
        ))
    }

    fn vcf_variant(spdi: &str) -> VcfVariant {
        spdi.parse().unwrap()
    }

    /// Length and FNV-1a hash per custom field, to keep full sequences out of snapshots.
    fn summarize_custom_fields(
        anns: &[AnnField],
    ) -> Vec<(String, Vec<Consequence>, BTreeMap<String, String>)> {
        anns.iter()
            .map(|ann| {
                let fields = ann
                    .custom_fields
                    .iter()
                    .map(|(key, value)| {
                        let value = value.as_deref().unwrap_or_default();
                        let hash = value.bytes().fold(0xcbf2_9ce4_8422_2325_u64, |h, b| {
                            (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
                        });
                        (
                            key.clone(),
                            format!("len={} fnv1a={hash:016x}", value.len()),
                        )
                    })
                    .collect();
                (ann.feature_id.clone(), ann.consequences.clone(), fields)
            })
            .collect()
    }

    #[rstest::rstest]
    #[case("17:41215920:G:T")] // missense
    #[case("17:41215920:G:GA")] // frameshift insertion
    #[case("17:41197818:GG:G")] // frameshift deletion
    fn report_sequences(#[case] spdi: &str) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", spdi.replace(':', "-"));

        let predictor = brca1_sequence_predictor()?;
        let res = predictor.predict(&vcf_variant(spdi))?.unwrap();

        insta::assert_yaml_snapshot!(summarize_custom_fields(&res));

        Ok(())
    }

    #[test]
    fn report_sequences_compound() -> Result<(), anyhow::Error> {
        let predictor = brca1_sequence_predictor()?;
        let res = predictor
            .predict_multiple(&[
                vcf_variant("17:41197701:G:C"),
                vcf_variant("17:41215920:G:T"),
            ])?
            .unwrap();

        insta::assert_yaml_snapshot!(summarize_custom_fields(&res));

        Ok(())
    }

    /// Records how many transcripts each call receives.
    #[derive(Debug, Default)]
    struct CountingPlugin {
        calls: std::sync::Mutex<Vec<usize>>,
    }

    impl TranscriptPlugin for CountingPlugin {
        fn columns(&self) -> Vec<String> {
            vec!["tx_count".into()]
        }

        fn annotate(
            &self,
            _variants: &[VcfVariant],
            txs: &mut Vec<TxAnnotation>,
        ) -> Result<(), SeqvarsError> {
            self.calls.lock().unwrap().push(txs.len());
            let count = txs.len().to_string();
            for tx in txs {
                tx.ann
                    .custom_fields
                    .insert("tx_count".into(), Some(count.clone()));
            }
            Ok(())
        }
    }

    #[test]
    fn plugin_runs_once_per_variant_on_reported_transcripts() -> Result<(), anyhow::Error> {
        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let plugin = Arc::new(CountingPlugin::default());
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .transcript_settings(TranscriptSettings::most_severe(PickUnit::Gene))
                .plugins(vec![Arc::new(SequencePlugin::AaAlt), plugin.clone()])
                .build()?,
        );

        assert_eq!(predictor.custom_columns(), [ANN_AA_SEQ_ALT, "tx_count"]);

        let res = predictor.predict(&vcf_variant("17:41215920:G:T"))?.unwrap();

        assert_eq!(res.len(), 1);
        assert_eq!(*plugin.calls.lock().unwrap(), [1]);
        assert_eq!(res[0].custom_fields["tx_count"].as_deref(), Some("1"));

        Ok(())
    }

    #[test]
    fn plugins_must_fill_distinct_columns() {
        let config = ConfigBuilder::default()
            .plugins(vec![
                Arc::new(CountingPlugin::default()),
                Arc::new(CountingPlugin::default()),
            ])
            .build();

        assert!(config.is_err());
    }

    #[test]
    fn flag_pick_column_follows_the_plugin_columns() -> Result<(), anyhow::Error> {
        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let plugin = Arc::new(CountingPlugin::default());
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .transcript_settings(TranscriptSettings {
                    flag_pick: true,
                    ..TranscriptSettings::most_severe(PickUnit::Gene)
                })
                .plugins(vec![plugin.clone()])
                .build()?,
        );

        assert_eq!(predictor.custom_columns(), ["tx_count", ANN_PICK]);

        let res = predictor.predict(&vcf_variant("17:41215920:G:T"))?.unwrap();
        let picked = res
            .iter()
            .filter(|ann| ann.custom_fields.get(ANN_PICK) == Some(&Some("1".into())))
            .count();

        assert_eq!(picked, 1);
        assert_eq!(*plugin.calls.lock().unwrap(), [res.len()]);

        Ok(())
    }

    /// Fills the column of `--flag-pick`.
    #[derive(Debug)]
    struct PickPlugin;

    impl TranscriptPlugin for PickPlugin {
        fn columns(&self) -> Vec<String> {
            vec![ANN_PICK.into()]
        }

        fn annotate(
            &self,
            _variants: &[VcfVariant],
            _txs: &mut Vec<TxAnnotation>,
        ) -> Result<(), SeqvarsError> {
            Ok(())
        }
    }

    #[test]
    fn flag_pick_and_plugins_must_fill_distinct_columns() {
        let settings = TranscriptSettings {
            flag_pick: true,
            ..TranscriptSettings::most_severe(PickUnit::Gene)
        };
        let with_flag = ConfigBuilder::default()
            .transcript_settings(settings)
            .plugins(vec![Arc::new(PickPlugin)])
            .build();
        let without_flag = ConfigBuilder::default()
            .plugins(vec![Arc::new(PickPlugin)])
            .build();

        assert!(with_flag.is_err());
        assert!(without_flag.is_ok());
    }

    /// A phased insertion in front of an SNV on a plus strand transcript. Both edits start at
    /// the same transcript index, and the inserted bases stay in front of the changed base.
    #[test]
    fn annotate_multiple_insertion_before_snv() -> Result<(), anyhow::Error> {
        // The genome around the exons, with the two genome-only bases of exon 1 and without the
        // three transcript-only bases of exon 2. The transcript ends at n.33, before the `A`
        // bases that the stored sequence appends.
        let tx_db = tx_db_with_indels_in_exon_cigars();
        let tx_seq = &tx_db.seq_db.as_ref().unwrap().seqs[0][..33];
        let genome = [
            "T".repeat(1000),
            tx_seq[..6].to_string(),
            "TT".into(),
            tx_seq[6..10].to_string(),
            "T".repeat(988),
            tx_seq[10..20].to_string(),
            tx_seq[23..].to_string(),
            "T".repeat(100),
        ]
        .concat();
        let mut reference = NamedTempFile::new()?;
        std::io::Write::write_all(&mut reference, format!(">1\n{genome}\n").as_bytes())?;
        let provider = Arc::new(MehariProvider::new(tx_db, Some(reference.path()), true));
        let predictor = ConsequencePredictor::new(provider, Default::default());

        let var = |position, reference: &str, alternative: &str| VcfVariant {
            chromosome: "1".into(),
            position,
            reference: reference.into(),
            alternative: alternative.into(),
        };
        // c.13_14insGGG and c.14C>A turn CCC (Pro5) into CGG GAC (Arg Asp).
        let res = predictor
            .predict_multiple(&[var(2005, "C", "CGGG"), var(2006, "C", "A")])?
            .unwrap();

        assert_eq!(res[0].hgvs_g.as_deref(), Some("g.2005_2006delinsCGGGA"));
        assert_eq!(res[0].hgvs_p.as_deref(), Some("p.Pro5delinsArgAsp"));

        Ok(())
    }

    /// Phased annotation reports the transcript strand like `predict` does.
    /// GRCh37, BRCA1, NM_007294.4 (MANE, reverse).
    #[test]
    fn annotate_multiple_brca1_minus_strand() -> Result<(), anyhow::Error> {
        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(provider, Default::default());

        let vars = [(41197707, "G", "T"), (41197711, "G", "A")].map(
            |(position, reference, alternative)| VcfVariant {
                chromosome: "17".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            },
        );
        let res = predictor.predict_multiple(&vars)?.unwrap();

        assert!(!res.is_empty());
        for ann in &res {
            assert_eq!(ann.strand, -1, "feature_id = {}", ann.feature_id);
        }
        // Combines c.5576C>T and c.5580C>A in transcript orientation.
        let mane = res
            .iter()
            .find(|ann| ann.feature_id == "NM_007294.4")
            .unwrap();
        assert_eq!(mane.hgvs_c.as_deref(), Some("c.5576_5580delinsTCCAA"));
        assert_eq!(
            mane.hgvs_p.as_deref(),
            Some("p.Pro1859_His1860delinsLeuGln")
        );

        Ok(())
    }

    #[rstest::rstest]
    #[case("n.2C>T", Some("ATGTAC"))]
    #[case("n.2_3del", Some("ATAC"))]
    #[case("n.2_3insGG", Some("ACGGGTAC"))]
    #[case("n.2_3dupCG", Some("ACGCGTAC"))]
    #[case("n.2+1C>T", None)] // intronic offset
    #[case("n.7C>T", None)] // outside the sequence
    fn apply_n_edits_to_sequence(
        #[case] edit: &str,
        #[case] expected: Option<&str>,
    ) -> Result<(), anyhow::Error> {
        let var_n = HgvsVariant::from_str(&format!("NM_000000.1:{edit}"))?;
        assert_eq!(apply_n_edits("ACGTAC", &[&var_n]).as_deref(), expected);

        Ok(())
    }

    /// An insertion in front of a changed base starts at the same index as the change. The
    /// inserted bases stay in front of the changed bases, whatever the input order.
    #[rstest::rstest]
    #[case(&["n.2_3insGG", "n.3G>T"], "ACGGTTAC")]
    #[case(&["n.3G>T", "n.2_3insGG"], "ACGGTTAC")]
    #[case(&["n.2_3insGG", "n.3_4del"], "ACGGAC")]
    fn apply_n_edits_insertion_before_change(
        #[case] edits: &[&str],
        #[case] expected: &str,
    ) -> Result<(), anyhow::Error> {
        let vars_n = edits
            .iter()
            .map(|edit| HgvsVariant::from_str(&format!("NM_000000.1:{edit}")))
            .collect::<Result<Vec<_>, _>>()?;
        let vars_n = vars_n.iter().collect::<Vec<_>>();
        assert_eq!(apply_n_edits("ACGTAC", &vars_n).as_deref(), Some(expected));

        Ok(())
    }

    /// Exons at n.1_10, n.11_20 and n.21_30.
    #[rstest::rstest]
    #[case("n.5_10del", false)] // ends with exon 1
    #[case("n.10_11del", true)] // last base of exon 1, first base of exon 2
    #[case("n.8_25del", true)]
    #[case("n.10_11insA", false)] // no reference bases
    #[case("n.10+1G>A", true)] // intronic offset
    #[case("n.28_32del", false)] // from the last exon past the transcript end
    fn alt_depends_on_splicing_in_exons(
        #[case] edit: &str,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let alignment = GenomeAlignment {
            exons: [(1, 10), (11, 20), (21, 30)]
                .into_iter()
                .enumerate()
                .map(|(ord, (start, end))| ExonAlignment {
                    ord: ord as i32,
                    alt_cds_start_i: Some(start),
                    alt_cds_end_i: Some(end),
                    ..Default::default()
                })
                .collect(),
            ..Default::default()
        };
        let var_n = HgvsVariant::from_str(&format!("NM_000000.1:{edit}"))?;
        assert_eq!(var_n.alt_depends_on_splicing(&alignment), expected);

        Ok(())
    }

    /// Annotates `spdis` (several: as phased variants) with reference and alternative sequences
    /// reported, and returns the annotation for `tx_id`. NM_007294.4 is BRCA1 (minus strand, CDS
    /// at n.114_5705).
    fn annotate_with_sequences(spdis: &[&str], tx_id: &str) -> Result<AnnField, anyhow::Error> {
        let vars = spdis
            .iter()
            .map(|spdi| {
                let spdi = spdi.split(':').collect::<Vec<_>>();
                Ok(VcfVariant {
                    chromosome: spdi[0].to_string(),
                    position: spdi[1].parse()?,
                    reference: spdi[2].to_string(),
                    alternative: spdi[3].to_string(),
                })
            })
            .collect::<Result<Vec<_>, anyhow::Error>>()?;
        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .plugins(vec![
                    Arc::new(SequencePlugin::TxRef),
                    Arc::new(SequencePlugin::TxAlt),
                    Arc::new(SequencePlugin::AaRef),
                    Arc::new(SequencePlugin::AaAlt),
                ])
                .build()?,
        );

        let anns = match vars.as_slice() {
            [var] => predictor.predict(var)?,
            _ => predictor.predict_multiple(&vars)?,
        };
        anns.unwrap_or_default()
            .into_iter()
            .find(|ann| ann.feature_id == tx_id)
            .ok_or_else(|| anyhow::anyhow!("no annotation for {tx_id}"))
    }

    fn custom_field<'a>(ann: &'a AnnField, key: &str) -> Option<&'a str> {
        ann.custom_fields.get(key).and_then(|v| v.as_deref())
    }

    #[test]
    fn alt_sequences_cds_snv() -> Result<(), anyhow::Error> {
        let ann = annotate_with_sequences(&["17:41197701:G:C"], "NM_007294.4")?;
        assert_eq!(ann.hgvs_n.as_deref(), Some("n.5699C>G"));
        assert_eq!(ann.hgvs_p.as_deref(), Some("p.His1862Gln"));

        let mut tx_alt = custom_field(&ann, ANN_TX_SEQ_REF).unwrap().to_string();
        tx_alt.replace_range(5698..5699, "G");
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), Some(tx_alt.as_str()));
        let mut aa_alt = custom_field(&ann, ANN_AA_SEQ_REF).unwrap().to_string();
        aa_alt.replace_range(1861..1862, "Q");
        assert_eq!(custom_field(&ann, ANN_AA_SEQ_ALT), Some(aa_alt.as_str()));

        Ok(())
    }

    #[test]
    fn alt_sequences_utr5_snv() -> Result<(), anyhow::Error> {
        let ann = annotate_with_sequences(&["17:41277340:T:C"], "NM_007294.4")?;
        assert_eq!(ann.hgvs_n.as_deref(), Some("n.42A>G"));
        assert_eq!(ann.hgvs_c.as_deref(), Some("c.-72A>G"));

        let mut tx_alt = custom_field(&ann, ANN_TX_SEQ_REF).unwrap().to_string();
        tx_alt.replace_range(41..42, "G");
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), Some(tx_alt.as_str()));
        assert_eq!(
            custom_field(&ann, ANN_AA_SEQ_ALT),
            custom_field(&ann, ANN_AA_SEQ_REF)
        );

        Ok(())
    }

    /// An insertion between c.-1 and c.1 lies in the 5' UTR: the start codon stays intact.
    #[test]
    fn alt_sequences_utr5_insertion_before_cds() -> Result<(), anyhow::Error> {
        let ann = annotate_with_sequences(&["17:41276113:T:TG"], "NM_007294.4")?;
        assert_eq!(ann.hgvs_n.as_deref(), Some("n.113_114insC"));
        assert_eq!(ann.hgvs_c.as_deref(), Some("c.-1_1insC"));

        let mut tx_alt = custom_field(&ann, ANN_TX_SEQ_REF).unwrap().to_string();
        tx_alt.insert(113, 'C');
        assert_eq!(&tx_alt[113..117], "CATG");
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), Some(tx_alt.as_str()));
        assert_eq!(
            custom_field(&ann, ANN_AA_SEQ_ALT),
            custom_field(&ann, ANN_AA_SEQ_REF)
        );

        Ok(())
    }

    /// A dup of the A at c.1 or at the last CDS base c.5592 leaves the CDS intact, see
    /// `annotate_ins_keeping_cds`. The alternative transcript has one more A at `index`.
    #[rstest::rstest]
    #[case("17:41276112:A:AT", "c.1dup", 113)]
    #[case("17:41197695:T:TT", "c.5592dup", 5704)]
    fn alt_sequences_dup_keeping_cds(
        #[case] spdi: &str,
        #[case] hgvs_c: &str,
        #[case] index: usize,
    ) -> Result<(), anyhow::Error> {
        let ann = annotate_with_sequences(&[spdi], "NM_007294.4")?;
        assert_eq!(ann.hgvs_c.as_deref(), Some(hgvs_c));

        let mut tx_alt = custom_field(&ann, ANN_TX_SEQ_REF).unwrap().to_string();
        tx_alt.insert(index, 'A');
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), Some(tx_alt.as_str()));
        assert_eq!(
            custom_field(&ann, ANN_AA_SEQ_ALT),
            custom_field(&ann, ANN_AA_SEQ_REF)
        );

        Ok(())
    }

    #[test]
    fn alt_sequences_utr3_deletion() -> Result<(), anyhow::Error> {
        let ann = annotate_with_sequences(&["17:41196499:CAAA:C"], "NM_007294.4")?;
        assert_eq!(ann.hgvs_n.as_deref(), Some("n.6898_6900del"));
        assert_eq!(ann.hgvs_c.as_deref(), Some("c.*1193_*1195del"));

        let mut tx_alt = custom_field(&ann, ANN_TX_SEQ_REF).unwrap().to_string();
        tx_alt.replace_range(6897..6900, "");
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), Some(tx_alt.as_str()));
        assert_eq!(
            custom_field(&ann, ANN_AA_SEQ_ALT),
            custom_field(&ann, ANN_AA_SEQ_REF)
        );

        Ok(())
    }

    /// With an intronic offset, the alternative transcript is unknown.
    #[rstest::rstest]
    #[case("17:41197837:G:A", "c.5468-18C>T")] // CDS intron
    #[case("17:41277000:G:C", "c.-20+288C>G")] // 5' UTR intron
    fn alt_sequences_intronic(
        #[case] spdi: &str,
        #[case] hgvs_c: &str,
    ) -> Result<(), anyhow::Error> {
        let ann = annotate_with_sequences(&[spdi], "NM_007294.4")?;
        assert_eq!(ann.hgvs_c.as_deref(), Some(hgvs_c));

        assert!(custom_field(&ann, ANN_TX_SEQ_REF).is_some());
        assert!(custom_field(&ann, ANN_AA_SEQ_REF).is_some());
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), None);
        assert_eq!(custom_field(&ann, ANN_AA_SEQ_ALT), None);

        Ok(())
    }

    /// A deletion from exon 13 into exon 14 (0-based `ord`) of NM_130837.3 (OPA1, plus strand)
    /// covers the whole intron between them. Its c. positions have no intronic offset, but the
    /// alternative transcript is unknown. The test database has no genome sequence, so the
    /// intron bases of the VCF REF allele are made up.
    #[test]
    fn alt_sequences_deletion_covers_intron() -> Result<(), anyhow::Error> {
        let reference = format!("TAAT{}ACT", "A".repeat(83));
        let ann = annotate_with_sequences(&[&format!("3:193361230:{reference}:T")], "NM_130837.3")?;
        assert_eq!(ann.hgvs_n.as_deref(), Some("n.1545_1550del"));

        assert!(custom_field(&ann, ANN_TX_SEQ_REF).is_some());
        assert!(custom_field(&ann, ANN_AA_SEQ_REF).is_some());
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), None);
        assert_eq!(custom_field(&ann, ANN_AA_SEQ_ALT), None);

        Ok(())
    }

    /// The phased path applies the n. edits of all variants to the transcript sequence.
    #[test]
    fn alt_sequences_phased_cds() -> Result<(), anyhow::Error> {
        let ann = annotate_with_sequences(&["17:41197701:G:C", "17:41197709:GG:G"], "NM_007294.4")?;
        assert_eq!(ann.hgvs_n.as_deref(), Some("n.5691_5699delinsACAGCCAG"));
        assert_eq!(ann.hgvs_p.as_deref(), Some("p.His1860ThrfsTer62"));

        let mut tx_alt = custom_field(&ann, ANN_TX_SEQ_REF).unwrap().to_string();
        tx_alt.replace_range(5698..5699, "G");
        tx_alt.replace_range(5690..5691, "");
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), Some(tx_alt.as_str()));

        Ok(())
    }

    /// Deletions within the start codon, after 5' UTRs of 0 to 3 bases.
    #[rstest::rstest]
    #[case::utr0_c1_3del("ATGTAG", 0, 0..3, false)]
    #[case::utr1_c1del_after_a("AATGGCCTAA", 1, 1..2, true)]
    #[case::utr1_c1del_after_c("CATGGCCTAA", 1, 1..2, false)]
    #[case::utr2_c1del_after_a("CAATGGCCTAA", 2, 2..3, true)]
    #[case::utr2_c3del("CCATGCCTAA", 2, 4..5, false)]
    #[case::utr3_c1del_after_a("GCAATGGCCTAA", 3, 3..4, true)]
    #[case::utr3_c1del_after_c("GCCATGGCCTAA", 3, 3..4, false)]
    #[case::utr3_c1_2del_repeat("GATATGGCCTAA", 3, 3..5, true)]
    #[case::utr3_c1_3del("GCAATGTGCTAA", 3, 3..6, false)]
    fn deletion_in_start_codon_keeps_cds(
        #[case] tx_seq: &str,
        #[case] cds_start: usize,
        #[case] del: std::ops::Range<usize>,
        #[case] expected: bool,
    ) {
        assert_eq!(deletion_keeps_cds(tx_seq, cds_start, del), expected);
    }

    /// A coding transcript `NM_000001.1` on the plus strand of chr1. Exon 1 spans
    /// chr1:1001-1060, the intron chr1:1061-1160 and exon 2 chr1:1161-1220. The CDS is `ATG`,
    /// 18 `GCC` codons and `TAA`, and starts `cds_start` bases into the transcript.
    fn two_exon_tx(cds_start: usize) -> (Transcript, String) {
        let utr = "CAGT".repeat(15);
        let cds = format!("ATG{}TAA", "GCC".repeat(18));
        let seq = format!("{}{cds}{}", &utr[..cds_start], &utr[cds_start..]);
        let start = i32::try_from(cds_start).unwrap();
        let genomic = |tx_pos: i32| {
            if tx_pos < 60 {
                1000 + tx_pos
            } else {
                1100 + tx_pos
            }
        };
        let exon = |ord: i32| ExonAlignment {
            alt_start_i: 1000 + 160 * ord,
            alt_end_i: 1060 + 160 * ord,
            ord,
            alt_cds_start_i: Some(1 + 60 * ord),
            alt_cds_end_i: Some(60 + 60 * ord),
            cigar: "60M".into(),
        };
        let tx = Transcript {
            id: "NM_000001.1".into(),
            gene_symbol: "GENE1".into(),
            gene_id: "1".into(),
            biotype: TranscriptBiotype::Coding.into(),
            protein: Some("NP_000001.1".into()),
            start_codon: Some(start),
            stop_codon: Some(start + 60),
            genome_alignments: vec![GenomeAlignment {
                genome_build: "grch38".into(),
                contig: "NC_000001.11".into(),
                cds_start: Some(genomic(start)),
                cds_end: Some(genomic(start + 59) + 1),
                strand: Strand::Plus.into(),
                exons: vec![exon(0), exon(1)],
                ..Default::default()
            }],
            filtered: Some(false),
            ..Default::default()
        };
        (tx, seq)
    }

    /// The range of an insertion is the base in front of it. The insertion moves the bases
    /// behind it, so a splice window that counts from its left end, e.g. donor +3 to +8 on the
    /// plus strand, contains an insertion if it contains the base behind it. A window that
    /// counts from its right end contains an insertion if it contains the base in front of it.
    ///
    /// `two_exon_tx(20)` has the donor +1 at chr1:1061 and the acceptor -1 at chr1:1160.
    #[rstest::rstest]
    #[case::exon_and_donor_1(1060, &[Consequence::ExonicSpliceRegionVariant])]
    #[case::donor_1_2(1061, &[Consequence::SpliceDonorVariant])]
    #[case::donor_2_3(
        1062,
        &[Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant]
    )]
    #[case::donor_4_5(
        1064,
        &[
            Consequence::SpliceRegionVariant,
            Consequence::SpliceDonorRegionVariant,
            Consequence::SpliceDonorFifthBaseVariant,
        ]
    )]
    #[case::donor_5_6(
        1065,
        &[Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant]
    )]
    #[case::donor_6_7(1066, &[Consequence::SpliceRegionVariant])]
    #[case::donor_8_9(1068, &[])]
    #[case::acceptor_9_8(1152, &[Consequence::SplicePolypyrimidineTractVariant])]
    #[case::acceptor_3_2(
        1158,
        &[Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant]
    )]
    #[case::acceptor_1_and_exon(1160, &[Consequence::ExonicSpliceRegionVariant])]
    #[case::exon_2_3(1162, &[Consequence::ExonicSpliceRegionVariant])]
    #[case::exon_3_4(1163, &[])]
    fn splice_windows_for_insertions(
        #[case] position: i32,
        #[case] expected: &[Consequence],
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = two_exon_tx(20);
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: "C".into(),
                alternative: "CTTT".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        let splice_terms = [
            Consequence::SpliceDonorVariant,
            Consequence::SpliceAcceptorVariant,
            Consequence::SpliceRegionVariant,
            Consequence::ExonicSpliceRegionVariant,
            Consequence::SpliceDonorRegionVariant,
            Consequence::SpliceDonorFifthBaseVariant,
            Consequence::SplicePolypyrimidineTractVariant,
        ];
        let found = splice_terms
            .iter()
            .filter(|term| ann.consequences.contains(term))
            .collect::<Vec<_>>();
        let expected = splice_terms
            .iter()
            .filter(|term| expected.contains(term))
            .collect::<Vec<_>>();
        assert_eq!(found, expected, "{ann:?}");
        Ok(())
    }

    /// hgvs-rs takes the intron bases of a dup next to an exon from the variant. Without them,
    /// the protein change is unknown.
    #[test]
    fn safe_project_c_to_p_dup_without_intron_bases() -> Result<(), anyhow::Error> {
        let (tx, seq) = two_exon_tx(20);
        let var_c = HgvsVariant::from_str("NM_000001.1:c.40+1dup")?;
        let var_p = predictor_for(tx, seq).safe_project_c_to_p(&var_c)?;
        assert_eq!(
            var_p.map(|var_p| FormattedLoc(&var_p).to_string()),
            Some("p.?".to_string())
        );
        Ok(())
    }

    /// An insertion or dup shifts the frame if it inserts a number of bases that is not a
    /// multiple of 3.
    #[rstest::rstest]
    #[case("c.10_11insTAA", false)]
    #[case("c.10_11insTA", true)]
    #[case("c.10_12dup", false)]
    #[case("c.10_23dup", true)]
    fn insertion_shifts_frame(
        #[case] var_c: &str,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let var_c = HgvsVariant::from_str(&format!("NM_000001.1:{var_c}"))?;
        assert_eq!(inserts_frameshift(&var_c), expected);
        Ok(())
    }

    /// Positions in the start codon are c.1 to c.3 without intron offsets.
    #[rstest::rstest]
    #[case("c.1_3del", true)]
    #[case("c.2del", true)]
    #[case("c.-1_1del", false)]
    #[case("c.3_4del", false)]
    #[case("c.1+50del", false)]
    #[case("c.1+4_1+7del", false)]
    fn deletion_within_first_codon(
        #[case] var_c: &str,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let HgvsVariant::CdsVariant { loc_edit, .. } =
            HgvsVariant::from_str(&format!("NM_000001.1:{var_c}"))?
        else {
            anyhow::bail!("not a c. variant: {var_c}");
        };
        let loc = loc_edit.loc.inner();
        assert_eq!(
            within_first_codon(loc.start.clone(), loc.end.clone()),
            expected
        );
        Ok(())
    }
}
