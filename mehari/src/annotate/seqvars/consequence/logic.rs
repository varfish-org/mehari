//! Compute molecular consequence of variants.
use super::pick::{self, Candidate};
use super::plugin::{TxAnnotation, TxContext};
use super::terms::{ANN_PICK, FeatureTag, GroupedAlleles};
use super::terms::{
    Allele, AnnField, Consequence, FeatureBiotype, FeatureType, Message, Pos, Rank, SoFeature,
};
use crate::annotate::cli::HgvsParentheses;
use crate::annotate::seqvars::consequence::{Config, Consequences, FormattedLoc, VcfVariant};
use crate::annotate::seqvars::provider::PbsTranscriptExt;
use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::errors::{GroupValidationError, SeqvarsError};
use crate::pbs::txs::{GenomeAlignment, Strand, Transcript, TranscriptBiotype, TranscriptTag};
use hgvs::mapper::altseq::{RefTranscriptData, VariantLocation, ref_transcript_data_shared};
use hgvs::parser::{NoRef, ProteinEdit};
use hgvs::sequences::{revcomp, translate_cds_with_exceptions};
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
    /// The region of the transcript that `c` changes, as hgvs-rs uses it to build `p`.
    region: Option<VariantLocation>,
    /// In VEP compatibility mode, the g. variant at the input position, before the 3' shift. VEP
    /// checks with it whether an indel overlaps the stop codon.
    unshifted_g: Option<HgvsVariant>,
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

    /// Whether VEP builds codons for the c. variant, see `vep_indel_protein_terms`. VEP does so
    /// if the variant has no intron offset, or if its bases go to an exon edge, which VEP adds to
    /// the exon. These are an insertion at the edge, and a dup that ends next to the edge, as VEP
    /// inserts the copy behind the dup.
    fn vep_has_codons(&self) -> bool {
        let Some(HgvsVariant::CdsVariant { loc_edit, .. }) = &self.c else {
            return false;
        };
        let loc = loc_edit.loc.inner();
        let offsets = (loc.start.offset.unwrap_or(0), loc.end.offset.unwrap_or(0));
        match loc_edit.edit.inner() {
            NaEdit::Ins { .. } => matches!(offsets, (0, 0) | (0, 1) | (-1, 0)),
            NaEdit::Dup { .. } => matches!(offsets.1, 0 | -1),
            _ => offsets == (0, 0),
        }
    }

    /// Whether the bases that the c. variant puts in hold another base than A, C, G and T, e.g.
    /// the `N` of an unknown stretch. VEP translates no such bases.
    fn has_ambiguous_alt(&self) -> bool {
        let Some(HgvsVariant::CdsVariant { loc_edit, .. }) = &self.c else {
            return false;
        };
        match loc_edit.edit.inner() {
            NaEdit::RefAlt { alternative, .. }
            | NaEdit::NumAlt { alternative, .. }
            | NaEdit::Ins { alternative } => !is_unambiguous(alternative),
            _ => false,
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

    /// Check if the variant is an insertion between an exon's edge base and the intron base next
    /// to it, e.g. `c.80_80+1insA` or `c.81-1_81insA`.
    fn is_exon_edge_insertion(&self) -> bool;

    /// The exon base in front of the exon edge where an n. variant can put its bases. These are
    /// an insertion at the edge, e.g. `n.80_80+1insA` or `n.81-1_81insA`, and a dup of intron
    /// bases next to it, e.g. `n.80+1dup` or `n.81-2_81-1dup`, whose copy can go there. The
    /// result is `80` in these examples, and `None` for other variants.
    fn exon_edge_gap(&self) -> Option<i32>;

    /// VEP's `stop_gained`, `stop_lost`, `stop_retained_variant` and `synonymous_variant` for an
    /// indel in the CDS, as in VEP 116.3.
    ///
    /// VEP translates only the codons that the variant changes
    /// (`TranscriptVariationAllele::codon`): the reference codons, and the same stretch of the
    /// altered CDS and 3' UTR. A partial codon at the end of this stretch gives an `X`. A stop
    /// that a new frame reaches behind these codons thus gives no `stop_gained`. After an `X`,
    /// VEP reads the bases at the old stop codon instead (`_ins_del_stop_altered`). For
    /// `stop_retained_variant`, VEP then checks the overlap with the stop codon at `unshifted_g`,
    /// the g. variant at the input position, before the 3' shift, see `vep_overlaps_stop_codon`.
    ///
    /// VEP puts an insertion at an exon edge into the exon (`Mapper::map_insert`): `c.N_N+1ins`
    /// behind c.N, `c.N-1_Nins` in front of c.N. For another indel with an intronic offset, VEP
    /// builds no codons. Of the stop terms, VEP then reports only `stop_retained_variant`: if the
    /// indel covers a CDS base, and if `unshifted_g` overlaps the stop codon.
    ///
    /// VEP translates no bases other than A, C, G and T (`seq_is_unambiguous_dna`), e.g. the `N`
    /// of an unknown stretch, and finds no stop codon altered by them. Of the stop terms, such an
    /// indel can then only get `stop_retained_variant`, from the overlap at `unshifted_g`.
    ///
    /// VEP also takes `inframe_insertion` and `protein_altering_variant` of an in-frame insertion
    /// from these codons, see `VepCodonTerms::inframe_insertion`.
    ///
    /// Returns `None` for an edit other than an insertion, duplication, deletion or delins, and
    /// for an edit that keeps the length.
    fn vep_indel_protein_terms(
        &self,
        unshifted_g: &HgvsVariant,
        ref_data: &RefTranscriptData,
        tx: &Transcript,
    ) -> Option<VepCodonTerms>;
}

/// VEP's terms of an indel from the codons that it changes, see `vep_indel_protein_terms`.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct VepCodonTerms {
    /// `stop_gained`, `stop_lost`, `stop_retained_variant` or `synonymous_variant`, or none.
    stops: Consequences,
    /// For an in-frame insertion, VEP's `inframe_insertion` and `protein_altering_variant`. The
    /// first applies if the altered peptide, cut behind its first stop, starts or ends with the
    /// reference peptide. Else the second applies, unless the uncut altered peptide starts with a
    /// stop, or starts or ends with the reference peptide. Else neither applies, e.g. for `Y` to
    /// `*Y`.
    ///
    /// `None` for another indel, and where other rules decide: in or right behind the start
    /// codon, in the stop codon, and for a peptide with an `X`.
    inframe_insertion: Option<Consequences>,
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

    fn is_exon_edge_insertion(&self) -> bool {
        let (start, end, edit) = match self {
            HgvsVariant::CdsVariant { loc_edit, .. } => {
                let loc = loc_edit.loc.inner();
                (loc.start.offset, loc.end.offset, loc_edit.edit.inner())
            }
            HgvsVariant::TxVariant { loc_edit, .. } => {
                let loc = loc_edit.loc.inner();
                (loc.start.offset, loc.end.offset, loc_edit.edit.inner())
            }
            _ => return false,
        };
        edit.is_ins() && matches!((start.unwrap_or(0), end.unwrap_or(0)), (0, 1) | (-1, 0))
    }

    fn exon_edge_gap(&self) -> Option<i32> {
        let HgvsVariant::TxVariant { loc_edit, .. } = self else {
            return None;
        };
        let loc = loc_edit.loc.inner();
        let offsets = (loc.start.offset.unwrap_or(0), loc.end.offset.unwrap_or(0));
        match (loc_edit.edit.inner(), offsets) {
            (NaEdit::Ins { .. }, (0, 1)) | (NaEdit::Dup { .. }, (1, _)) => Some(loc.start.base),
            (NaEdit::Ins { .. }, (-1, 0)) | (NaEdit::Dup { .. }, (_, -1)) => Some(loc.end.base - 1),
            _ => None,
        }
    }

    fn vep_indel_protein_terms(
        &self,
        unshifted_g: &HgvsVariant,
        ref_data: &RefTranscriptData,
        tx: &Transcript,
    ) -> Option<VepCodonTerms> {
        let codon = |cds_pos: i32| (cds_pos + 2) / 3;
        let alignment = tx.genome_alignments.first()?;
        let incomplete_3p = tx.is_incomplete_3p();
        let overlaps_stop_cil = !incomplete_3p && vep_overlaps_stop_codon(unshifted_g, alignment);

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
            let edit = loc_edit.edit.inner();
            let behind_intron = match (loc.start.offset.unwrap_or(0), loc.end.offset.unwrap_or(0)) {
                (0, 0) => false,
                (0, 1) if edit.is_ins() => false,
                (-1, 0) if edit.is_ins() => true,
                _ => return None,
            };
            let (start, end) = (cds_pos(&loc.start), cds_pos(&loc.end));
            let (start, end, alt) = match edit {
                // `c.N-1_Nins` lies behind c.N-1, the end of the exon in front.
                NaEdit::Ins { alternative } => {
                    let behind = start - i32::from(behind_intron);
                    (behind + 1, behind, alternative.clone())
                }
                NaEdit::Dup { .. } => (
                    end + 1,
                    end,
                    cds_and_utr.get(index(start - 1)?..index(end)?)?.to_string(),
                ),
                NaEdit::DelRef { .. } | NaEdit::DelNum { .. } => (start, end, String::new()),
                NaEdit::RefAlt { alternative, .. } => (start, end, alternative.clone()),
                _ => return None,
            };
            Some((start, end, alt))
        };
        // The CDS positions of the exon ends.
        let exon_ends = (alignment.exons.iter())
            .filter_map(|exon| Some(exon.alt_cds_end_i? - ref_data.cds_start + 1))
            .collect::<Vec<_>>();
        let Some((mut start, mut end, mut alt)) = changed(self) else {
            // VEP builds no codons for another indel with an intronic offset.
            let HgvsVariant::CdsVariant { loc_edit, .. } = self else {
                return None;
            };
            let indel = match loc_edit.edit.inner() {
                NaEdit::RefAlt {
                    reference,
                    alternative,
                } => reference.len() != alternative.len(),
                NaEdit::NumAlt { count, alternative } => {
                    usize::try_from(*count).ok() != Some(alternative.len())
                }
                NaEdit::InvRef { .. } | NaEdit::InvNum { .. } => false,
                NaEdit::DelRef { .. }
                | NaEdit::DelNum { .. }
                | NaEdit::Ins { .. }
                | NaEdit::Dup { .. } => true,
            };
            // Whether a deletion or delins covers an exon base in `from..=to`. A position with
            // an intronic offset lies between two exon bases.
            let covers = |variant: &HgvsVariant, from: i32, to: i32| {
                let HgvsVariant::CdsVariant { loc_edit, .. } = variant else {
                    return false;
                };
                if loc_edit.edit.inner().is_ins() || loc_edit.edit.inner().is_dup() {
                    return false;
                }
                let loc = loc_edit.loc.inner();
                let start = cds_pos(&loc.start) + i32::from(loc.start.offset.unwrap_or(0) > 0);
                let end = cds_pos(&loc.end) - i32::from(loc.end.offset.unwrap_or(0) < 0);
                start <= to && end >= from
            };
            let stop_retained = covers(self, 1, cds_len) && overlaps_stop_cil;
            return indel.then(|| VepCodonTerms {
                stops: if stop_retained {
                    Consequence::StopRetainedVariant.into()
                } else {
                    Consequences::empty()
                },
                inframe_insertion: None,
            });
        };
        // VEP shifts an insertion or deletion 3' on the genome, also across the CDS end, where
        // hgvs-rs stops. On the genome, an intron follows the end of an exon.
        let base = |cds_pos: i32| cds_and_utr.as_bytes().get(index(cds_pos - 1)?).copied();
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
            return Some(VepCodonTerms::default());
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
        let unambiguous = is_unambiguous(alt);
        let peptides = if in_cds && unambiguous {
            let alt_len = index(3 * (last - first + 1) + len_change)?;
            // VEP takes the reference codons from the CDS alone, so a partial last codon stays
            // partial.
            let ref_codons = cds_and_utr.get(from..index((3 * last).min(cds_len))?)?;
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

        let overlaps_stop = !incomplete_3p && end >= cds_len - 2 && start <= cds_len;
        // VEP's `_ins_del_stop_altered`: the altered CDS and 3' UTR read no stop codon at the
        // position of the old one. Without an altered sequence, VEP counts the stop codon as
        // kept.
        let stop_altered = |overlaps: bool| {
            let stop_codon = index(cds_len - 3)?..index(cds_len)?;
            Some(
                overlaps
                    && unambiguous
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
        // See `VepCodonTerms::inframe_insertion`. An insertion lies between `end` and `start`.
        let inframe_insertion = peptides
            .as_ref()
            .filter(|(ref_pep, alt_pep)| {
                start > end
                    && len_change % 3 == 0
                    && end > 3
                    && !ref_pep.contains(['*', 'X'])
                    && !alt_pep.contains('X')
            })
            .map(|(ref_pep, alt_pep)| {
                let up_to_stop = alt_pep.split_inclusive('*').next().unwrap_or_default();
                let keeps_ref = |pep: &str| {
                    pep.starts_with(ref_pep.as_str()) || pep.ends_with(ref_pep.as_str())
                };
                if keeps_ref(up_to_stop) {
                    Consequence::InframeInsertion.into()
                } else if !alt_pep.starts_with('*') && !keeps_ref(alt_pep) {
                    Consequence::ProteinAlteringVariant.into()
                } else {
                    Consequences::empty()
                }
            });

        Some(VepCodonTerms {
            stops: if stop_retained {
                Consequence::StopRetainedVariant.into()
            } else if stop_lost {
                Consequence::StopLost.into()
            } else if stop_gained {
                Consequence::StopGained.into()
            } else if synonymous {
                Consequence::SynonymousVariant.into()
            } else {
                Consequences::empty()
            },
            inframe_insertion,
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

#[derive(Debug, Clone)]
struct TranscriptLocationContext {
    rank: Rank,
    distance: Option<i32>,
    is_exonic: bool,
    is_intronic: bool,
    is_upstream: bool,
    is_downstream: bool,
}

/// The terms of a variant on a coding transcript by their source. `consequences_fix_special_cases`
/// reconciles them.
#[derive(Debug)]
struct ConsequenceContext {
    /// Terms of the c. variant, e.g. `frameshift_variant` by its length.
    cds_consequences: Consequences,
    /// Terms read off the p. variant, e.g. `frameshift_variant` for `fs`.
    protein_consequences: Consequences,
    /// Terms of the protein the variant produces compared to the reference, see
    /// `sequence_consequences()`, for the variants whose terms depend on it.
    sequence_consequences: Consequences,
    /// In VEP compatibility mode, VEP's terms of an indel from its changed codons, see
    /// `vep_indel_protein_terms()`.
    vep_codon_consequences: Option<VepCodonTerms>,
    cdna_pos: Option<Pos>,
    cds_pos: Option<Pos>,
    protein_pos: Option<Pos>,
}

/// The consequences of a variant on a transcript, with the location and the HGVS projection they
/// come from.
struct TxConsequences {
    consequences: Consequences,
    location: TranscriptLocationContext,
    /// `None` up- or downstream of the transcript.
    projection: Option<(HgvsProjectionContext, ConsequenceContext)>,
}

/// The terms of the essential splice sites.
const SPLICE_SITES: Consequences =
    enumflags2::make_bitflags!(Consequence::{SpliceDonorVariant | SpliceAcceptorVariant});

/// The terms of the splice sites and the regions around them.
const SPLICE_TERMS: Consequences = enumflags2::make_bitflags!(Consequence::{
    SpliceDonorVariant
        | SpliceAcceptorVariant
        | SpliceRegionVariant
        | ExonicSpliceRegionVariant
        | SpliceDonorRegionVariant
        | SpliceDonorFifthBaseVariant
        | SplicePolypyrimidineTractVariant
});

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
            // The bases of an insertion or dup at an exon edge are intronic if they keep the
            // splice site right next to the exon, else they extend the exon, also where the
            // intron bases are unknown. VEP compatibility keeps them in the exon.
            ins_at_boundary_by_sequence: !config.vep_compatibility_mode,
            inferred_p_is_uncertain: config.hgvs_parentheses == HgvsParentheses::Predicted,
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

    /// Whether to report VEP terms. VEP compatibility mode implies them.
    fn vep_consequence_terms(&self) -> bool {
        self.config.vep_consequence_terms || self.config.vep_compatibility_mode
    }

    /// The HGNC ID of the gene of transcript `tx_id`, if it has one.
    pub fn hgnc_id(&self, tx_id: &str) -> Option<&str> {
        self.provider.get_tx(tx_id)?.hgnc_id.as_deref()
    }

    /// The context of a transcript annotation; empty if no plugin reads it.
    fn tx_context(
        &self,
        tx_id: Option<String>,
        var_c: Option<HgvsVariant>,
        var_n: Option<HgvsVariant>,
    ) -> TxContext {
        if self.config.plugins.is_empty() {
            TxContext::empty()
        } else {
            TxContext::new(self.provider.clone(), tx_id, var_c, var_n)
        }
    }

    /// The context of the annotation of `projection` on transcript `tx_id`, which takes the
    /// variants of `projection` over; empty if no plugin reads it.
    fn projection_context(
        &self,
        tx_id: String,
        projection: Option<HgvsProjectionContext>,
        alignment: &GenomeAlignment,
    ) -> TxContext {
        if self.config.plugins.is_empty() {
            return TxContext::empty();
        }
        let edge_in_intron = (projection.as_ref()).and_then(|p| self.edge_bases_in_intron(p));
        let ins_at_boundary_is_intronic =
            edge_in_intron.unwrap_or(self.mapper.config.ins_at_boundary_is_intronic);
        let alt_depends_on_splicing = (projection.as_ref())
            .and_then(|p| p.n.as_ref())
            .is_some_and(|var_n| {
                // Unless they count as intronic, the bases of an insertion at an exon edge join
                // the exon and leave the splice site intact. Where hgvs-rs decides this per
                // variant, so do those of the copy of a dup that can go there.
                let at_edge = match edge_in_intron {
                    Some(_) => var_n.exon_edge_gap().is_some(),
                    None => var_n.is_exon_edge_insertion(),
                };
                var_n.alt_depends_on_splicing(alignment)
                    && (ins_at_boundary_is_intronic || !at_edge)
            });
        let (var_c, var_n, region) =
            projection.map_or((None, None, None), |p| (p.c, p.n, p.region));
        TxContext::new(self.provider.clone(), Some(tx_id), var_c, var_n).with_alt_inputs(
            region,
            alt_depends_on_splicing,
            ins_at_boundary_is_intronic,
        )
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
        // Normalize variant by stripping common prefix and suffix. An `N` in ALT is one unknown
        // base, also an ALT of `N` alone.
        let norm_var = self.normalize_variant(var);
        // A variant without change has nothing to annotate. hgvs-rs would put the empty change
        // at both ends of the contig.
        if norm_var.reference.is_empty() && norm_var.alternative.is_empty() {
            return Ok(Some(
                self.run_plugins(std::slice::from_ref(var), Vec::new())?,
            ));
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
        let (var_g_fwd, var_g_rev) = self.normalize_g(&var_g, chrom_acc)?;
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
        txs.retain(|tx| {
            let var_g = if tx.alt_strand == -1 {
                &var_g_rev
            } else {
                &var_g_fwd
            };
            self.in_reach(&tx.tx_ac, var_g)
        });

        // Handle case of no overlapping transcripts -> intergenic.
        if txs.is_empty() {
            if !self.config.keep_intergenic {
                // Plugins may want predictions without annotations, too.
                return Ok(Some(
                    self.run_plugins(std::slice::from_ref(var), Vec::new())?,
                ));
            }
            let hgvs_g = self
                .writes_hgvs(&var.alternative)
                .then(|| self.format_hgvs_g(&var_g));

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

        // An indel in a repeat has two ends, one 3'-shifted on each strand.
        let shifts = var_g_fwd != var_g_rev;
        // The same for all transcripts on a strand.
        let hgvs_g_fwd = self.format_hgvs_g(&var_g_fwd);
        let hgvs_g_rev = if shifts {
            self.format_hgvs_g(&var_g_rev)
        } else {
            hgvs_g_fwd.clone()
        };
        let anns = self.select_annotations(txs, |tx| {
            if tx.alt_strand == -1 {
                let other_g = shifts.then_some(&var_g_fwd);
                self.build_ann_field(
                    var,
                    &var_g_rev,
                    &hgvs_g_rev,
                    &var_g,
                    other_g,
                    tx.clone(),
                    var_start_rev,
                    var_end_rev,
                )
            } else {
                let other_g = shifts.then_some(&var_g_rev);
                self.build_ann_field(
                    var,
                    &var_g_fwd,
                    &hgvs_g_fwd,
                    &var_g,
                    other_g,
                    tx.clone(),
                    var_start_fwd,
                    var_end_fwd,
                )
            }
        })?;
        Ok(Some(self.run_plugins(std::slice::from_ref(var), anns)?))
    }

    /// Whether the transcript `tx_ac` is in reach of `var_g`, shifted 3' on its strand. The
    /// transcripts come from `PADDING` around the variant at both ends of its repeat, but the 3'
    /// shift, or the copy that a dup inserts behind its bases, can put the bases farther away. A
    /// variant that lies outside the transcript is then out of reach, as if the transcript were
    /// not nearby.
    fn in_reach(&self, tx_ac: &str, var_g: &HgvsVariant) -> bool {
        let Some(alignment) =
            (self.provider.get_tx(tx_ac)).and_then(|tx| tx.genome_alignments.first())
        else {
            return true;
        };
        let Ok(strand) = Strand::try_from(alignment.strand) else {
            return true;
        };
        let location = |var_g: &HgvsVariant| {
            let (start, end) = Self::get_var_start_end(var_g);
            (self.determine_transcript_context(alignment, strand, var_g, start, end, false)).0
        };
        let outside =
            |location: &TranscriptLocationContext| location.is_upstream || location.is_downstream;
        let own = location(var_g);
        if !outside(&own) {
            return true;
        }
        let bases = Self::dup_as_ins(var_g, strand).map_or(own, |ins_g| location(&ins_g));
        !outside(&bases)
            || bases
                .distance
                .is_none_or(|distance| distance.abs() <= PADDING)
    }

    /// `var_g` shifted 3' for a transcript on the plus strand and on the minus strand, if
    /// `renormalize_g` is set.
    fn normalize_g(
        &self,
        var_g: &HgvsVariant,
        chrom_acc: &str,
    ) -> Result<(HgvsVariant, HgvsVariant), SeqvarsError> {
        if !self.mapper.config.renormalize_g {
            return Ok((var_g.clone(), var_g.clone()));
        }
        // On circular contigs, normalize one lap further so that shifting 5' can cross the
        // origin, then shift back.  Positions may then lie just beyond either contig end.
        let lap = self.circular_length(chrom_acc).unwrap_or(0);
        let var_g_lap = Self::shift_var_g(var_g, lap);
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
        Ok((
            Self::shift_var_g(&right, -lap),
            Self::shift_var_g(&left, -lap),
        ))
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
                Some(tx) => Some(Candidate::new(record, tx, self.provider.gene_key(tx))),
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

    /// The insertion of the copy that the dup `var_g` adds behind its bases on a transcript on
    /// `strand`.
    ///
    /// Returns `None` if `var_g` is no dup.
    fn dup_as_ins(var_g: &HgvsVariant, strand: Strand) -> Option<HgvsVariant> {
        let HgvsVariant::GenomeVariant {
            accession,
            gene_symbol,
            loc_edit,
        } = var_g
        else {
            return None;
        };
        let NaEdit::Dup { reference } = loc_edit.edit.inner() else {
            return None;
        };
        let loc = loc_edit.loc.inner();
        // Behind the bases on the transcript is in front of them on the genome on the minus
        // strand.
        let pos = if strand == Strand::Minus {
            loc.start? - 1
        } else {
            loc.end?
        };
        Some(HgvsVariant::GenomeVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit: GenomeLocEdit {
                loc: Mu::Certain(GenomeInterval {
                    start: Some(pos),
                    end: Some(pos + 1),
                }),
                edit: Mu::Certain(NaEdit::Ins {
                    alternative: reference.clone(),
                }),
            },
        })
    }

    /// The insertion at an exon edge that equals the dup `var_g` whose bases span that edge: its
    /// copy, rotated to start with the bases behind the edge, e.g. `c.186-1_186insTAG` for
    /// `c.186-2_186dup` of `AGT`. `None` for another variant, and for a dup across no edge or
    /// more than one.
    fn dup_across_edge_as_ins(
        var_g: &HgvsVariant,
        alignment: &GenomeAlignment,
    ) -> Option<HgvsVariant> {
        let HgvsVariant::GenomeVariant {
            accession,
            gene_symbol,
            loc_edit,
        } = var_g
        else {
            return None;
        };
        let NaEdit::Dup { reference } = loc_edit.edit.inner() else {
            return None;
        };
        let loc = loc_edit.loc.inner();
        let (start, end) = (loc.start?, loc.end?);
        // The 0-based positions between two bases of the dup that lie at an exon edge.
        let mut edges = (alignment.exons.iter().tuple_windows())
            .flat_map(|(prev, next)| [prev.alt_end_i, next.alt_start_i])
            .filter(|&pos| (start..end).contains(&pos));
        let (Some(edge), None) = (edges.next(), edges.next()) else {
            return None;
        };
        let split = usize::try_from(edge - start + 1).ok()?;
        let copy = [reference.get(split..)?, reference.get(..split)?].concat();
        Some(HgvsVariant::GenomeVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit: GenomeLocEdit {
                loc: Mu::Certain(GenomeInterval {
                    start: Some(edge),
                    end: Some(edge + 1),
                }),
                edit: Mu::Certain(NaEdit::Ins { alternative: copy }),
            },
        })
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
                // The range of an insertion is the base in front of it, as its end is the base
                // behind it. A dup covers its bases.
                let end = if edit.is_ins() { end - 1 } else { end };
                (start, end)
            }
            _ => unreachable!(),
        }
    }

    /// Whether the 0-based position `pos` between two bases lies between an exon and an intron
    /// of `alignment`.
    fn at_intron_edge(alignment: &GenomeAlignment, pos: i32) -> bool {
        (alignment.exons.iter().tuple_windows())
            .any(|(prev, next)| prev.alt_end_i == pos || next.alt_start_i == pos)
    }

    /// The location of `var_g` on the transcript of `alignment`, and the terms of that location.
    /// With `ins_in_intron`, an insertion between an exon and an intron puts its bases into the
    /// intron, else into the exon.
    fn determine_transcript_context(
        &self,
        alignment: &GenomeAlignment,
        strand: Strand,
        var_g: &HgvsVariant,
        var_start: i32,
        var_end: i32,
        ins_in_intron: bool,
    ) -> (TranscriptLocationContext, Consequences) {
        let mut consequences = Consequences::empty();
        let mut rank = Rank::default();
        let mut is_exonic = false;
        let mut is_intronic = false;
        let mut distance: Option<i32> = None;

        // The range of an insertion is the base in front of it. An insertion at the start of
        // an exon adds its bases to the exon, so its range is the first exon base instead. This
        // includes the start of the transcript on the genome, so that an insertion at either
        // transcript end joins the exon, as its HGVS description has it. With `ins_in_intron`,
        // an insertion between an exon and an intron adds them to the intron, so its range is
        // empty, at the intron edge. In VEP compatibility mode, an insertion in front of the
        // first or behind the last exon on the genome lies outside the transcript, as in VEP
        // 116.3, so its range is the base outside. Without a reference, mehari does not shift an
        // insertion on the genome, so an insertion that VEP shifts behind the last exon stays
        // inside.
        let is_ins = matches!(var_g, HgvsVariant::GenomeVariant { loc_edit, .. }
            if matches!(loc_edit.edit.inner(), NaEdit::Ins { .. }));
        let in_intron = is_ins && ins_in_intron && Self::at_intron_edge(alignment, var_end);
        let at_exon_start = (alignment.exons.iter())
            .skip(usize::from(self.config.vep_compatibility_mode))
            .any(|exon| exon.alt_start_i == var_end);
        let behind_last_exon = self.config.vep_compatibility_mode
            && alignment.exons.iter().map(|exon| exon.alt_end_i).max() == Some(var_end);
        let (var_start, var_end) = if in_intron {
            (var_end, var_end)
        } else if is_ins && (at_exon_start || behind_last_exon) {
            (var_start + 1, var_end + 1)
        } else {
            (var_start, var_end)
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
                    self.config.vep_compatibility_mode,
                );
            }

            min_start = Some(std::cmp::min(min_start.unwrap_or(exon_start), exon_start));
            max_end = Some(std::cmp::max(max_end.unwrap_or(exon_end), exon_end));
            prev_end = Some(exon_end);
        }

        // An MNV changes only the bases where REF and ALT differ, so only these get splice terms.
        if let Some(runs) = var_g
            .na_edit()
            .and_then(|edit| changed_runs(edit, var_start))
        {
            consequences &= !SPLICE_TERMS;
            for (start, end) in runs {
                consequences |= self.splice_terms(alignment, strand, start, end);
            }
        }

        // Bases in the intron at an exon edge keep the splice site, but they lie in the splice
        // region.
        if in_intron {
            consequences |= Consequence::SpliceRegionVariant;
        }

        // VEP 116.3 reports no polypyrimidine tract for a variant that overlaps an exon.
        if self.config.vep_compatibility_mode && is_exonic {
            consequences.remove(Consequence::SplicePolypyrimidineTractVariant);
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
    ///
    /// `ins_c` is the c. insertion of the copy of a dup, if `var_g` is one, see `dup_as_ins`.
    /// In VEP compatibility mode, `p.` comes from it.
    fn project_hgvs(
        &self,
        var_g: &HgvsVariant,
        unshifted_g: &HgvsVariant,
        ins_c: Option<&HgvsVariant>,
        tx: &Transcript,
        transcript_biotype: TranscriptBiotype,
    ) -> Result<HgvsProjectionContext, SeqvarsError> {
        let mut projection = HgvsProjectionContext {
            n: None,
            c: None,
            p: None,
            region: None,
            unshifted_g: None,
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
                // In VEP compatibility mode, `p.` of a dup comes from the insertion of its copy, as
                // in VEP 116.3. VEP puts the copy behind the dup, also next to an exon, while
                // hgvs-rs puts it in front of the dup.
                let vep = self.config.vep_compatibility_mode;
                let var_c = if vep { ins_c.unwrap_or(var_c) } else { var_c };
                if let Some((var_p, region)) = self.safe_project_c_to_p(var_c)? {
                    projection.p = Some(var_p);
                    projection.region = region;
                }
                self.mark_unknown_protein(&mut projection, &tx.id);
            }
        }

        if self.config.vep_compatibility_mode && transcript_biotype == TranscriptBiotype::Coding {
            projection.unshifted_g = Some(unshifted_g.clone());
        }

        Ok(projection)
    }

    /// The positions and consequences of `projection` at `transcript_location`. The terms of a
    /// dup come from the insertion of its copy: `copy_location` is where the copy goes, `ins_c`
    /// the insertion as a c. variant.
    #[allow(clippy::too_many_arguments)]
    fn analyze_transcript_consequences(
        &self,
        projection: &HgvsProjectionContext,
        tx: &Transcript,
        tx_record: &TxForRegionRecord,
        transcript_location: &TranscriptLocationContext,
        copy_location: Option<&TranscriptLocationContext>,
        ins_c: Option<&HgvsVariant>,
        tx_len: i32,
        transcript_biotype: TranscriptBiotype,
    ) -> Result<ConsequenceContext, SeqvarsError> {
        let mut context = ConsequenceContext {
            cds_consequences: Consequences::empty(),
            protein_consequences: Consequences::empty(),
            sequence_consequences: Consequences::empty(),
            vep_codon_consequences: None,
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
            // before c.1, e.g., `c.1del` or `c.1_3dup`. In VEP compatibility mode, the position of
            // such an indel counts, not its 5' UTR region.
            let region = projection.region.filter(|&region| {
                !(self.config.vep_compatibility_mode && region == VariantLocation::FivePrimeUtr)
            });
            context.cds_consequences = Self::analyze_cds_variant(
                ins_c.unwrap_or(var_c),
                region,
                copy_location.unwrap_or(transcript_location).is_exonic,
                conservative,
                incomplete_3p,
                available_cds_len,
            );
            // A change of a start codon base that hgvs-rs places in the 5' UTR keeps the start
            // codon in front of the CDS, as the same change fits in front of c.1, e.g. `c.-2_1del`
            // or `c.1dup`, whose copy goes between c.1 and c.2. VEP compatibility mode decides
            // this as VEP, see `consequences_fix_special_cases`.
            if !self.config.vep_compatibility_mode
                && projection.region == Some(VariantLocation::FivePrimeUtr)
                && changes_start_codon(var_c)
            {
                context.cds_consequences |= Consequence::StartRetainedVariant;
            }

            if let Some(var_p) = &projection.p {
                // The terms of a frameshift, of a delins with a stop, of a change of the stop
                // codon and of a change that leaves no CDS codon behind it depend on the protein,
                // see `consequences_fix_special_cases`. Building the protein costs time, so only
                // these get it. A frameshift counts by p. (`fs`) or by its length, as p. can
                // describe it otherwise, e.g. `p.(Trp5Ser)`. Without the start codon
                // (`p.Met1?`), the protein is unknown.
                let start_lost = matches!(
                    var_p,
                    HgvsVariant::ProtVariant {
                        loc_edit: ProtLocEdit::InitiationUncertain,
                        ..
                    }
                );
                let fs = matches!(protein_edit(projection), Some(ProteinEdit::Fs { .. }))
                    || (context
                        .cds_consequences
                        .contains(Consequence::FrameshiftVariant)
                        && !start_lost);
                let delins_with_stop = matches!(protein_edit(projection),
                    Some(ProteinEdit::DelIns { alternative }) if has_stop(alternative));
                let cds_len = tx
                    .start_codon
                    .zip(tx.stop_codon)
                    .map(|(start, stop)| stop - start);
                // A protein that reads unknown bases gives no terms, see `project_hgvs`.
                let unknown_protein = !self.config.vep_compatibility_mode
                    && matches!(
                        var_p,
                        HgvsVariant::ProtVariant {
                            loc_edit: ProtLocEdit::Unknown,
                            ..
                        }
                    )
                    && projection.has_ambiguous_alt();
                // The terms that need the amino acids give way to their SO parent, unless the
                // frameshift or in-frame term of the c. variant stands for it. The bases of
                // `c.-1_1ins` go into the 5' UTR, also where hgvs-rs counts them to the CDS, in
                // front of a transcript without 5' UTR.
                if unknown_protein
                    && projection.region == Some(VariantLocation::Cds)
                    && context.cds_pos.is_some()
                    && projection.is_within_coding_sequence()
                    && !context.cds_consequences.intersects(
                        Consequence::FrameshiftVariant
                            | Consequence::ConservativeInframeInsertion
                            | Consequence::DisruptiveInframeInsertion
                            | Consequence::ConservativeInframeDeletion
                            | Consequence::DisruptiveInframeDeletion,
                    )
                {
                    context.cds_consequences |= Consequence::CodingSequenceVariant;
                }
                if (fs
                    || delins_with_stop
                    || touches_stop_codon(var_c, cds_len)
                    || self.leaves_no_codon_behind(projection, tx))
                    && !unknown_protein
                    && let Ok(ref_data) = self.ref_data(&tx_record.tx_ac)
                {
                    context.sequence_consequences = self.sequence_consequences(var_c, &ref_data);
                }
                if matches!(
                    var_p,
                    HgvsVariant::ProtVariant {
                        loc_edit: ProtLocEdit::Unknown,
                        ..
                    }
                ) {
                    // protein_pos and protein_consequences remain empty.
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

                    let deletes_codons = conservative
                        && matches!(
                            var_c.na_edit(),
                            Some(NaEdit::DelRef { .. } | NaEdit::DelNum { .. })
                        );
                    context.protein_consequences = Self::analyze_protein_variant(
                        var_p,
                        conservative,
                        deletes_codons,
                        incomplete_3p,
                    );

                    let indel_terms = Consequence::FrameshiftVariant
                        | Consequence::ConservativeInframeInsertion
                        | Consequence::DisruptiveInframeInsertion
                        | Consequence::ConservativeInframeDeletion
                        | Consequence::DisruptiveInframeDeletion;
                    let indel_in_cds = context.cds_consequences.intersects(indel_terms);
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
                    // An indel with a silent change keeps only the term of that change. VEP
                    // reports the indel terms instead. hgvs-rs writes no `p.=` for an indel in
                    // the CDS: an alt sequence that keeps the CDS makes a UTR variant, `p.?`.
                    if indel_in_cds && silent && !self.config.vep_compatibility_mode {
                        context.cds_consequences.remove(indel_terms);
                    }
                }
            }

            // VEP takes the stop terms and `synonymous_variant` of an indel from the changed
            // codons, whatever `p.` says. For a dup, these are the codons that the insertion of
            // its copy changes.
            if self.config.vep_compatibility_mode
                && let Ok(ref_data) = self.ref_data(&tx_record.tx_ac)
            {
                let var_c = ins_c.unwrap_or(var_c);
                // A phased change has no input position, so it stands for itself, its copy for a
                // dup.
                let unshifted_g = projection.unshifted_g.clone().or_else(|| {
                    let alignment = tx.genome_alignments.first()?;
                    let strand = Strand::try_from(alignment.strand).ok()?;
                    let var_g = (self.mapper.variant_mapper())
                        .c_to_g(var_c, &alignment.contig, ALT_ALN_METHOD)
                        .ok()?;
                    Some(Self::dup_as_ins(&var_g, strand).unwrap_or(var_g))
                });
                context.vep_codon_consequences = unshifted_g.and_then(|unshifted_g| {
                    var_c.vep_indel_protein_terms(&unshifted_g, &ref_data, tx)
                });
            }
        }

        Ok(context)
    }

    /// Consequences of `var_g` on `tx`, with the location and the HGVS projection they come
    /// from. `unshifted_g` is the variant at the input position, before the 3' shift. Returns
    /// `None` if the projection onto the transcript fails.
    #[allow(clippy::too_many_arguments)]
    fn tx_consequences(
        &self,
        var_g: &HgvsVariant,
        unshifted_g: &HgvsVariant,
        tx: &Transcript,
        tx_record: &TxForRegionRecord,
        alignment: &GenomeAlignment,
        strand: Strand,
        transcript_biotype: TranscriptBiotype,
        var_start: i32,
        var_end: i32,
    ) -> Result<Option<TxConsequences>, SeqvarsError> {
        let (mut transcript_location, transcript_consequences) =
            self.determine_transcript_context(alignment, strand, var_g, var_start, var_end, false);
        let outside = transcript_location.is_upstream || transcript_location.is_downstream;
        let coding = transcript_biotype == TranscriptBiotype::Coding;

        let mut consequences = transcript_consequences;

        // A dup inserts its copy behind its bases, so its terms come from there, e.g., the exon
        // terms for a dup of the first exon base. Its positions stay those of `var_g`.
        let mut ins_g = Self::dup_as_ins(var_g, strand);
        let mut copy_location = match &ins_g {
            Some(ins_g) => {
                let (ins_start, ins_end) = Self::get_var_start_end(ins_g);
                let (location, copy_consequences) = self.determine_transcript_context(
                    alignment, strand, ins_g, ins_start, ins_end, false,
                );
                consequences = copy_consequences;
                Some(location)
            }
            None => None,
        };

        // The insertion of the copy of a dup on a coding transcript, as a c. variant. The CDS
        // terms of the dup come from it, too.
        let ins_to_c = |ins_g: &HgvsVariant| {
            self.mapper
                .variant_mapper()
                .g_to_c(ins_g, &tx.id, ALT_ALN_METHOD)
                .inspect_err(|e| {
                    tracing::debug!("{}: cannot map {}: {}", &tx.id, ins_g, e);
                })
                .ok()
        };
        let mut ins_c = ins_g
            .as_ref()
            .filter(|_| coding && !outside)
            .and_then(ins_to_c);
        let projection = if outside {
            None
        } else {
            let projection =
                self.project_hgvs(var_g, unshifted_g, ins_c.as_ref(), tx, transcript_biotype)?;
            if projection.n.is_none() {
                return Ok(None); // Early exit if g->n projection failed.
            }
            Some(projection)
        };

        // hgvs-rs decides whether the bases of an insertion at an exon edge, or the copy of a dup
        // that can go there, lie in the intron or join the exon. The variant takes the location
        // of these bases at the edge, and their terms. The copy of a dup across the edge can go
        // there, rotated, so that the terms agree with the region of the dup and with its p.
        if let Some(in_intron) = projection
            .as_ref()
            .and_then(|p| self.edge_bases_in_intron(p))
        {
            let at_edge = |ins_g: &HgvsVariant| {
                Self::at_intron_edge(alignment, Self::get_var_start_end(ins_g).1)
            };
            // The copy can also go in front of the dup, which is behind it as seen from the other
            // strand.
            let opposite = if strand == Strand::Minus {
                Strand::Plus
            } else {
                Strand::Minus
            };
            let edge_g = match &ins_g {
                Some(behind) if at_edge(behind) => Some(behind.clone()),
                Some(_) => Self::dup_as_ins(var_g, opposite)
                    .filter(|front| at_edge(front))
                    .or_else(|| Self::dup_across_edge_as_ins(var_g, alignment)),
                None if var_g.na_edit().is_some_and(NaEdit::is_ins) && at_edge(var_g) => {
                    Some(var_g.clone())
                }
                None => None,
            };
            if let Some(edge_g) = edge_g {
                let (edge_start, edge_end) = Self::get_var_start_end(&edge_g);
                let (edge_location, edge_consequences) = self.determine_transcript_context(
                    alignment, strand, &edge_g, edge_start, edge_end, in_intron,
                );
                // The copy behind a dup keeps its terms if it lies on that side away from the
                // edge, e.g., behind a dup of donor +1 to +7.
                let keeps_copy = ins_g.as_ref().is_some_and(|behind| !at_edge(behind))
                    && copy_location
                        .as_ref()
                        .is_some_and(|location| location.is_intronic == in_intron);
                if !keeps_copy {
                    consequences = edge_consequences;
                    if ins_g.is_some() {
                        ins_c = ins_to_c(&edge_g);
                        copy_location = Some(edge_location.clone());
                        ins_g = Some(edge_g);
                    }
                }
                transcript_location = edge_location;
            }
        }

        let location = copy_location.as_ref().unwrap_or(&transcript_location);
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

        let projection = if let Some(projection) = projection {
            // VEP 116.3 puts an insertion between an exon and an intron outside the exon, also
            // the copy of a dup.
            if self.config.vep_compatibility_mode
                && transcript_biotype == TranscriptBiotype::NonCoding
            {
                let ins_n = match &ins_g {
                    Some(ins_g) => self
                        .mapper
                        .variant_mapper()
                        .g_to_n(ins_g, &tx.id, ALT_ALN_METHOD)
                        .ok(),
                    None => projection.n.clone(),
                };
                if ins_n.is_some_and(|var_n| var_n.is_exon_edge_insertion()) {
                    consequences.remove(Consequence::NonCodingTranscriptExonVariant);
                    consequences |= Consequence::NonCodingTranscriptVariant;
                }
            }

            let consequence_ctx = self.analyze_transcript_consequences(
                &projection,
                tx,
                tx_record,
                &transcript_location,
                copy_location.as_ref(),
                ins_c.as_ref(),
                tx.tx_len(),
                transcript_biotype,
            )?;

            // In VEP compatibility mode, the copy of a dup of the last transcript base lies
            // outside the transcript, so the dup gets no CDS or protein terms.
            let copy_outside = copy_location
                .as_ref()
                .is_some_and(|location| location.is_upstream || location.is_downstream);
            let (cds, protein, sequence, vep_codons) =
                if self.config.vep_compatibility_mode && copy_outside {
                    Default::default()
                } else {
                    (
                        consequence_ctx.cds_consequences,
                        consequence_ctx.protein_consequences,
                        consequence_ctx.sequence_consequences,
                        consequence_ctx.vep_codon_consequences,
                    )
                };
            consequences |= cds;

            self.consequences_fix_special_cases(
                &mut consequences,
                // exon_alignment_consequences, // TODO include these as well
                cds,
                protein,
                sequence,
                vep_codons,
                &projection,
                tx,
            );

            Some((projection, consequence_ctx))
        } else {
            None
        };

        Ok(Some(TxConsequences {
            consequences,
            location: transcript_location,
            projection,
        }))
    }

    /// Whether the same change at the other end of its repeat, `other_g`, leaves the essential
    /// splice sites intact, as an equivalent position of `var_g`.
    fn splice_sites_intact_elsewhere(
        &self,
        var_g: &HgvsVariant,
        other_g: Option<&HgvsVariant>,
        alignment: &GenomeAlignment,
        strand: Strand,
    ) -> bool {
        // The copy of a dup can also go in front of the dup on the transcript, which is behind it
        // as seen from the other strand.
        let opposite = if strand == Strand::Minus {
            Strand::Plus
        } else {
            Strand::Minus
        };
        let front_ins_g = Self::dup_as_ins(other_g.unwrap_or(var_g), opposite);
        front_ins_g.as_ref().or(other_g).is_some_and(|other_g| {
            let (other_start, other_end) = Self::get_var_start_end(other_g);
            let (_, other_consequences) = self.determine_transcript_context(
                alignment,
                strand,
                other_g,
                other_start,
                other_end,
                false,
            );
            !other_consequences.intersects(SPLICE_SITES)
        })
    }

    /// Annotation of `var_g` on the transcript of `tx_record`. `unshifted_g` is `var_g` at the
    /// input position, before the 3' shift.
    ///
    /// The HGVS descriptions, the positions and the consequence terms all come from `var_g`, the
    /// 3'-shifted variant, so that they agree. `hgvs_g` is `var_g` formatted, see
    /// `format_hgvs_g`. `other_g` is the same indel at the other end of its repeat, if that
    /// differs, see `Message::InfoSpliceSiteIntactAtEquivalentPosition`.
    #[allow(clippy::too_many_arguments)]
    fn build_ann_field(
        &self,
        orig_var: &VcfVariant,
        var_g: &HgvsVariant,
        hgvs_g: &str,
        unshifted_g: &HgvsVariant,
        other_g: Option<&HgvsVariant>,
        tx_record: TxForRegionRecord,
        var_start: i32,
        var_end: i32,
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

        let Some(TxConsequences {
            mut consequences,
            location: transcript_location,
            projection,
        }) = self.tx_consequences(
            var_g,
            unshifted_g,
            tx,
            &tx_record,
            alignment,
            strand,
            transcript_biotype,
            var_start,
            var_end,
        )?
        else {
            return Ok(None); // Early exit if g->n projection failed.
        };
        let (rank, cdna_pos, cds_pos, protein_pos) = match &projection {
            Some((_, ctx)) => (
                Some(transcript_location.rank.clone()),
                ctx.cdna_pos.clone(),
                ctx.cds_pos.clone(),
                ctx.protein_pos.clone(),
            ),
            None => (None, None, None, None),
        };
        let owned_projection = projection.map(|(projection, _)| projection);
        let projection = owned_projection.as_ref();

        let hgvs = self.writes_hgvs(&orig_var.alternative);
        let hgvs_g = hgvs.then(|| hgvs_g.to_string());
        let hgvs_projection = projection.filter(|_| hgvs);
        let hgvs_n = hgvs_projection
            .and_then(|p| p.n.as_ref())
            .map(|var| FormattedLoc(var).to_string());
        let hgvs_c = hgvs_projection
            .and_then(|p| p.c.as_ref())
            .map(|var| FormattedLoc(var).to_string());
        let hgvs_p = hgvs_projection
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

        if self.vep_consequence_terms() {
            self.adjust_vep_terms(&mut consequences, projection);
        }

        if consequences.is_empty() {
            tracing::warn!(
                "No consequences for {}:{} on {}, adding `gene_variant`; most likely the \
                transcript has multiple stop codons and the variant lies behind the first.",
                orig_var.chromosome,
                orig_var.position,
                &tx_record.tx_ac,
            );
            consequences |= Consequence::GeneVariant;
        }

        // The 3' shift can move an indel onto an essential splice site, although the same change
        // at the other end of its repeat leaves the site intact. The terms follow the HGVS
        // description; a message points this out.
        let messages = (consequences.intersects(SPLICE_SITES)
            && self.splice_sites_intact_elsewhere(var_g, other_g, alignment, strand))
        .then(|| vec![Message::InfoSpliceSiteIntactAtEquivalentPosition]);

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
            messages,
            custom_fields: BTreeMap::new(),
        };

        Ok(Some(TxAnnotation {
            context: self.projection_context(tx_record.tx_ac, owned_projection, alignment),
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

        // VEP's `missense_variant` excludes a lost start codon and the stop terms.
        if self.config.vep_compatibility_mode
            && consequences.intersects(StartLost | StopGained | StopLost | StopRetainedVariant)
        {
            consequences.remove(MissenseVariant);
        }
        // VEP's `frameshift` excludes `stop_retained_variant`.
        if self.config.vep_compatibility_mode && consequences.contains(StopRetainedVariant) {
            consequences.remove(FrameshiftVariant | FrameshiftElongation | FrameshiftTruncation);
        }

        let exonic_splice_region = consequences.contains(ExonicSpliceRegionVariant);
        if exonic_splice_region {
            consequences.remove(ExonicSpliceRegionVariant);
            consequences.insert(SpliceRegionVariant);
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

        // VEP builds no codons for a change of the CDS with an end in an intron, so it reports no
        // in-frame term there, but `coding_sequence_variant`, unless one of the terms below
        // applies (`coding_unknown` in VEP's `VariationEffect.pm`).
        let is_inframe = consequences.intersects(InframeDeletion | InframeInsertion);
        if self.config.vep_compatibility_mode
            && (exonic_splice_region || is_inframe)
            && projection_context.is_some_and(|p| p.is_within_cds_bounds() && !p.vep_has_codons())
        {
            consequences.remove(InframeDeletion | InframeInsertion);
            if !consequences.intersects(
                FrameshiftVariant
                    | ProteinAlteringVariant
                    | StartRetainedVariant
                    | StartLost
                    | StopRetainedVariant
                    | StopLost,
            ) {
                consequences.insert(CodingSequenceVariant);
            }
        }

        // VEP translates no bases other than A, C, G and T, so it gives such bases in the CDS no
        // terms of the changed amino acids. An in-frame deletion needs the changed codons to
        // match, which these bases never do. VEP then reports `coding_sequence_variant`, unless
        // one of the terms below applies (`coding_unknown`). Bases that go into a UTR, e.g. with
        // `c.-1_1ins`, are not in the CDS.
        let in_utr = |p: &HgvsProjectionContext| {
            matches!(
                p.region,
                Some(VariantLocation::FivePrimeUtr | VariantLocation::ThreePrimeUtr)
            )
        };
        if self.config.vep_compatibility_mode
            && projection_context.is_some_and(|p| {
                p.is_within_cds_bounds()
                    && !in_utr(p)
                    && p.vep_has_codons()
                    && p.has_ambiguous_alt()
            })
        {
            consequences.remove(
                InframeDeletion
                    | InframeInsertion
                    | ProteinAlteringVariant
                    | MissenseVariant
                    | StopGained
                    | SynonymousVariant,
            );
            if !consequences.intersects(
                FrameshiftVariant
                    | StartRetainedVariant
                    | StartLost
                    | StopRetainedVariant
                    | StopLost,
            ) {
                consequences.insert(CodingSequenceVariant);
            }
        }

        if consequences.contains(ExonLossVariant)
            && *consequences != Into::<Consequences>::into(ExonLossVariant)
        {
            consequences.remove(ExonLossVariant);
        }

        // VEP reports `feature_elongation` only for structural variants.
        if self.config.vep_compatibility_mode
            && consequences.contains(FeatureElongation)
            && *consequences != Into::<Consequences>::into(FeatureElongation)
        {
            consequences.remove(FeatureElongation);
        }

        // VEP's splice and intron terms exclude each other where SO's regions overlap.
        if !self.config.vep_compatibility_mode {
            return;
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

        // VEP counts the intron bases from +3 and -3 on as intronic (`_intron_effects`), so a
        // change of donor +5 stays an `intron_variant`.
        let suppress_intron = SpliceDonorVariant | SpliceAcceptorVariant;

        if consequences.intersects(suppress_intron) {
            consequences.remove(IntronVariant);
        }
    }

    /// Add the protein terms to `consequences`, which hold the terms of the location and of the
    /// c. variant (`consequences_cds`), and fix special cases. The protein terms are those read
    /// off the p. variant (`consequences_protein`), reconciled with those of the protein that
    /// the variant produces (`consequences_sequence`). In VEP compatibility mode, VEP's stop
    /// terms and `synonymous_variant` of an indel (`vep_codons`) replace those of the protein,
    /// and its in-frame insertion terms those of the c. and the p. variant. See
    /// `ConsequenceContext`.
    #[allow(clippy::too_many_arguments, unused_variables)]
    fn consequences_fix_special_cases(
        &self,
        consequences: &mut Consequences,
        consequences_cds: Consequences,
        consequences_protein: Consequences,
        consequences_sequence: Consequences,
        vep_codons: Option<VepCodonTerms>,
        projection: &HgvsProjectionContext,
        tx: &Transcript,
    ) {
        let vep_compat = self.config.vep_compatibility_mode;
        let no_codon_behind = self.leaves_no_codon_behind(projection, tx);
        let mut consequences_protein = consequences_protein;
        if no_codon_behind {
            // A change that leaves no CDS codon behind it has the terms of the protein it
            // produces, without `missense_variant` in VEP compatibility mode. HGVS can still
            // describe it as a frameshift, e.g., a deletion through the stop codon, or as `p.?`
            // past the transcript end.
            consequences_protein = consequences_sequence;
            if vep_compat {
                consequences_protein.remove(Consequence::MissenseVariant);
            }
        } else if consequences_protein.contains(Consequence::FrameshiftVariant) {
            // An `fs` that ends the protein in front of its stop codon or behind it.
            if consequences_sequence.contains(Consequence::StopGained) {
                consequences_protein |= Consequence::FrameshiftTruncation;
            }
            if consequences_sequence.contains(Consequence::FeatureElongation) {
                consequences_protein |= Consequence::FrameshiftElongation;
            }
            // An `fs` whose protein keeps its length changes residues, and another stop codon
            // takes the old position. VEP keeps only the frameshift.
            if consequences_sequence.contains(Consequence::MissenseVariant) && !vep_compat {
                consequences_protein |= consequences_sequence
                    & (Consequence::MissenseVariant | Consequence::StopRetainedVariant);
            }
        }
        // Inserted residues that end with a stop codon in front of the old one shorten the
        // protein, also where p. keeps the number of residues, e.g. `p.Lys3_Leu4delinsAsnTer`,
        // so no `missense_variant`. A new stop codon at the old position keeps the length: then
        // the changed residues are `missense_variant`, not `protein_altering_variant`, and the
        // stop is retained, not gained. VEP keeps the terms of p.
        if !vep_compat
            && matches!(protein_edit(projection), Some(ProteinEdit::DelIns { alternative })
                if has_stop(alternative))
        {
            if consequences_sequence.contains(Consequence::StopRetainedVariant) {
                consequences_protein
                    .remove(Consequence::StopGained | Consequence::ProteinAlteringVariant);
                consequences_protein |= consequences_sequence
                    & (Consequence::MissenseVariant | Consequence::StopRetainedVariant);
            } else {
                consequences_protein.remove(Consequence::MissenseVariant);
            }
        }
        // A change of stop codon bases that keeps the stop codon is `stop_retained_variant`,
        // also where p. does not show it, e.g. `p.(Trp5Cys)` for `c.15_17delinsCTG`, where TGA
        // takes the place of TAA.
        if !vep_compat
            && consequences_sequence.contains(Consequence::StopRetainedVariant)
            && !consequences_protein.intersects(
                Consequence::StopGained | Consequence::StopLost | Consequence::StopRetainedVariant,
            )
        {
            consequences_protein |= Consequence::StopRetainedVariant;
        }
        // An in-frame insertion in front of the stop codon keeps it, also where p. describes the
        // longer protein as an extension, e.g. `p.Ter6AspextTer1` for `c.15_16insGAT`. So it is
        // no `stop_lost`. VEP keeps it.
        if !vep_compat
            && matches!(protein_edit(projection), Some(ProteinEdit::Ext { .. }))
            && let (Some(var_c), Some(start), Some(stop)) =
                (&projection.c, tx.start_codon, tx.stop_codon)
            && inserts_in_frame_before_stop_codon(var_c, stop - start)
        {
            consequences_protein.remove(Consequence::StopLost);
        }
        if let Some(vep_codons) = vep_codons {
            consequences_protein.remove(
                Consequence::StopGained
                    | Consequence::StopLost
                    | Consequence::StopRetainedVariant
                    | Consequence::SynonymousVariant,
            );
            consequences_protein |= vep_codons.stops;
            if let Some(inframe_insertion) = vep_codons.inframe_insertion {
                let terms = Consequence::ConservativeInframeInsertion
                    | Consequence::DisruptiveInframeInsertion
                    | Consequence::ProteinAlteringVariant;
                consequences.remove(terms);
                consequences_protein.remove(terms);
                consequences_protein |= inframe_insertion;
            }
        }
        *consequences |= consequences_protein;

        // If we have a transcript_ablation, we can remove all other consequences
        if consequences.contains(Consequence::TranscriptAblation) {
            *consequences = Consequence::TranscriptAblation.into();
            return;
        }

        // Do not report splice variants in UTRs.
        if self.config.discard_utr_splice_variants {
            let splice_variants = Consequence::ExonicSpliceRegionVariant
                | Consequence::SpliceDonorVariant
                | Consequence::SpliceAcceptorVariant
                | Consequence::SpliceRegionVariant
                | Consequence::SpliceDonorFifthBaseVariant
                | Consequence::SpliceDonorRegionVariant
                | Consequence::SplicePolypyrimidineTractVariant;
            let utr_intron_variants =
                Consequence::FivePrimeUtrIntronVariant | Consequence::ThreePrimeUtrIntronVariant;
            let utr_exon_variants =
                Consequence::FivePrimeUtrExonVariant | Consequence::ThreePrimeUtrExonVariant;
            // As for the UTR terms, the region decides, not the c. description. Without a
            // region, e.g., if p. failed, both ends of the c. description must lie in one UTR, so
            // that a deletion from the 5' UTR into the CDS keeps its splice terms. Without c.,
            // the UTR terms decide.
            let in_utr = match projection.region {
                Some(region) => matches!(
                    region,
                    VariantLocation::FivePrimeUtr | VariantLocation::ThreePrimeUtr
                ),
                None => projection.c.as_ref().is_none_or(lies_in_one_utr),
            };
            if in_utr
                && consequences.intersects(splice_variants)
                && consequences.intersects(utr_intron_variants | utr_exon_variants)
            {
                *consequences &= !splice_variants;
            }
        }

        // `frameshift_variant` of the c. variant needs a codon of the CDS behind the change that
        // reads in another frame. Without the start codon (`p.Met1?`), no reading frame of the
        // protein is known to shift. VEP keeps it.
        if !vep_compat && (no_codon_behind || consequences_protein.contains(Consequence::StartLost))
        {
            *consequences &= !(consequences_cds & Consequence::FrameshiftVariant);
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

        // VEP reports `start_retained_variant` for a deletion in the start codon, or from the 5'
        // UTR into it, that leaves the CDS intact. It also reports `start_lost` there, as it reads
        // the codon at the old place of the start codon; mehari does not. Outside VEP
        // compatibility mode, hgvs-rs puts such a deletion in the 5' UTR.
        if self.config.vep_compatibility_mode
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

            // If edit occurs within the first 3 bases of the CDS, or a deletion from the 5' UTR
            // into them,
            let (start, end) = (c_loc.start.base, c_loc.end.base);
            let from_utr = CdsPos {
                base: start.max(1),
                ..c_loc.start.clone()
            };
            if within_first_codon(from_utr, c_loc.end.clone()) {
                // … then we need to check whether this is a start lost or a start retained.
                // To that end, check the deletion against the transcript sequence (0-based), in
                // which c.1 follows c.-1.
                let start_retained = match c_edit {
                    NaEdit::DelRef { .. } => match (
                        self.provider.get_seq_part(&accession.value, None, None),
                        usize::try_from(n_loc.start.base - start - i32::from(start < 0)),
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
            }
            if within_first_codon(c_loc.start.clone(), c_loc.end.clone()) {
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
        if tx.is_incomplete_5p() {
            *consequences &= !(Consequence::StartLost | Consequence::StartRetainedVariant);
            if consequences.is_empty() {
                *consequences |= Consequence::CodingSequenceVariant;
            }
        }
    }

    /// The splice terms of the bases from `start` to `end` (0-based, end-exclusive) on the
    /// transcript of `alignment`, as `determine_transcript_context` finds them for a variant.
    fn splice_terms(
        &self,
        alignment: &GenomeAlignment,
        strand: Strand,
        start: i32,
        end: i32,
    ) -> Consequences {
        let mut consequences = Consequences::empty();
        let mut prev_end = None;
        for exon in &alignment.exons {
            if overlaps(start, end, exon.alt_start_i, exon.alt_end_i) {
                let rank = Rank {
                    ord: exon.ord + 1,
                    total: alignment.exons.len() as i32,
                };
                consequences |= Self::analyze_exonic_variant(
                    0,
                    strand,
                    start,
                    end,
                    exon.alt_start_i,
                    exon.alt_end_i,
                    &rank,
                    false,
                );
            }
            if let Some(intron_start) = prev_end {
                consequences |= Self::analyze_intronic_variant(
                    0,
                    alignment,
                    strand,
                    start,
                    end,
                    intron_start,
                    exon.alt_start_i,
                    false,
                    self.config.vep_compatibility_mode,
                );
            }
            prev_end = Some(exon.alt_end_i);
        }
        consequences & SPLICE_TERMS
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
                } else if strand == Strand::Minus && !rank.is_first() {
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
        vep_compatibility_mode: bool,
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
        // VEP 116.3 puts an insertion into the donor region or at the fifth base only if the
        // bases on both sides of it are inside. So no insertion is at the fifth base.
        let var_overlaps_donor = |start: i32, end: i32| -> bool {
            if vep_compatibility_mode {
                var_overlaps(start, end) && var_overlaps_behind(start, end)
            } else if strand == Strand::Plus {
                var_overlaps_behind(start, end)
            } else {
                var_overlaps(start, end)
            }
        };
        // VEP 116.3 puts an insertion into the polypyrimidine tract if the base on either side of
        // it is inside.
        let var_overlaps_tract = |start: i32, end: i32| -> bool {
            if vep_compatibility_mode {
                var_overlaps(start, end) || var_overlaps_behind(start, end)
            } else if strand == Strand::Plus {
                var_overlaps(start, end)
            } else {
                var_overlaps_behind(start, end)
            }
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
        if strand == Strand::Plus && var_overlaps_tract(intron_end - 17, intron_end - 2) {
            consequences |= Consequence::SplicePolypyrimidineTractVariant;
        }
        if strand == Strand::Minus && var_overlaps_tract(intron_start + 2, intron_start + 17) {
            consequences |= Consequence::SplicePolypyrimidineTractVariant;
        }

        // Check conditions for splice_donor_region_variant
        // (A sequence variant that falls in the region between the 3rd and 6th base after splice junction (5' end of intron))
        // Note that this is two bases short of the intronic part of splice_region_variant
        if strand == Strand::Plus && var_overlaps_donor(intron_start + 2, intron_start + 6) {
            consequences |= Consequence::SpliceDonorRegionVariant;
        }
        if strand == Strand::Minus && var_overlaps_donor(intron_end - 6, intron_end - 2) {
            consequences |= Consequence::SpliceDonorRegionVariant;
        }

        // Check conditions for splice_donor_5th_base_variant
        // (A sequence variant that causes a change at the 5th base pair after the start of the intron in the orientation of the transcript.)
        if strand == Strand::Plus && var_overlaps_donor(intron_start + 4, intron_start + 5) {
            consequences |= Consequence::SpliceDonorFifthBaseVariant;
        }
        if strand == Strand::Minus && var_overlaps_donor(intron_end - 5, intron_end - 4) {
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
            } => i32::from(edit.inner().is_ins()),
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

            // An insertion at an exon edge (`c.10_10+1ins`, `c.11-1_11ins`) adds its bases to
            // the exon, or to the intron if they keep the splice site there, see
            // `edge_bases_in_intron()`. In front of the first CDS base (`c.1-1_1ins`), the exon
            // bases join the 5' UTR. Behind the last one (`c.N_N+1ins`), they join the 3' UTR.
            let at_edge = matches!(edit, NaEdit::Ins { .. })
                && start_base == end_base
                && start_cds_from == end_cds_from
                && matches!((loc_start_offset, loc_end_offset), (0, 1) | (-1, 0));
            let at_exon_edge = at_edge && is_exonic;
            let is_intronic_or_utr = is_intronic_or_utr || at_edge && !is_exonic;
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

            // As in VEP, only a partial last codon counts. It is partial if the CDS length is
            // not a multiple of 3, and it starts at the last codon boundary.
            if incomplete_3p
                && let Some(cds_len) = available_cds_len
                && cds_len % 3 != 0
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

            // A variant with both ends in introns covers the exon bases between them, if any. The
            // first lies behind the intron of its start, the last in front of that of its end.
            let between_introns = is_intronic_or_utr && is_exonic;
            let first_exon_base_in_5_prime_utr = start_cds_from == CdsFrom::Start
                && if loc_start_offset < 0 {
                    start_base < 0
                } else {
                    start_base < -1
                };
            let last_exon_base_in_3_prime_utr =
                end_cds_from == CdsFrom::End && (loc_end_offset > 0 || end_base > 1);

            // Detect variants affecting the 5'/3' UTRs.
            if (starts_left_of_start && start_base < 0) || in_5_prime_utr || ins_in_front_of_start {
                if is_intronic_or_utr {
                    consequences |= Consequence::FivePrimeUtrIntronVariant;
                } else if is_exonic {
                    consequences |= Consequence::FivePrimeUtrExonVariant;
                }
                if between_introns && first_exon_base_in_5_prime_utr {
                    consequences |= Consequence::FivePrimeUtrExonVariant;
                }
            }
            if ends_right_of_stop {
                if is_intronic_or_utr {
                    consequences |= Consequence::ThreePrimeUtrIntronVariant;
                } else if is_exonic {
                    consequences |= Consequence::ThreePrimeUtrExonVariant;
                }
                if between_introns && last_exon_base_in_3_prime_utr {
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
            // Between two exon positions, a deletion or delins replaces the CDS bases from one to
            // the other. Its reference also holds the bases of an intron in between.
            let cds_bases = (loc_start_offset == 0 && loc_end_offset == 0)
                .then(|| start_base.abs_diff(end_base) as usize + 1);

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
                        let reference_len = cds_bases.unwrap_or(reference.len());
                        if reference_len.abs_diff(alternative.len()) % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else {
                            // Check for inframe insertions/deletions (that are not delins)
                            match (
                                reference_len.cmp(&alternative.len()),
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
                        if cds_bases.unwrap_or(reference.len()) % 3 != 0 {
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
                        if cds_bases.map_or(count % 3, |bases| bases as i32 % 3) != 0 {
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

    /// Terms read off the p. variant `var_p`, see `consequences_fix_special_cases` for those
    /// that the protein decides. `deletes_codons` tells whether the c. variant deletes whole
    /// codons.
    fn analyze_protein_variant(
        var_p: &HgvsVariant,
        conservative: bool,
        deletes_codons: bool,
        incomplete_3p: bool,
    ) -> Consequences {
        let mut consequences: Consequences = Consequences::empty();

        // TODO move to hgvs-rs library as method of `ProtPos` or similar
        // `X` (`Xaa`) is an unknown residue, e.g., of a codon with an `N`, not a stop.
        fn is_stop(s: &str) -> bool {
            s == "Ter" || s == "*"
        }

        match var_p {
            HgvsVariant::ProtVariant { loc_edit, .. } => match loc_edit {
                ProtLocEdit::Ordinary { loc, edit } => {
                    let loc = loc.inner();
                    match edit.inner() {
                        ProteinEdit::Fs { .. } => {
                            consequences |= Consequence::FrameshiftVariant;
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
                                    // A deletion of whole codons keeps the others, so a stop
                                    // in place of a residue means that the deletion ends in
                                    // front of the stop codon. It only shortens the protein,
                                    // e.g. 1:43450470:GCCT:G, ENST00000634258.3:
                                    // c.10294_10296del, p.Leu3432Ter.
                                    if deletes_codons {
                                        consequences |= Consequence::ConservativeInframeDeletion;
                                    } else {
                                        consequences |= Consequence::StopGained;
                                    }
                                }
                            } else if alternative == "X" {
                                // An unknown residue can be the same, another or a stop, so
                                // only the common SO parent applies.
                                consequences |= Consequence::CodingSequenceVariant;
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
                        // not a multiple of 3, as HGVS prescribes. The frame still shifts then,
                        // so the CDS terms tell `frameshift_variant` or the in-frame insertion.
                        ProteinEdit::Ins { alternative } if has_stop(alternative) => {
                            consequences |= Consequence::StopGained;
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

    /// Whether the annotations of a variant with ALT `alternative` get HGVS descriptions. VEP
    /// 116.3 writes none for an allele with bases other than A, C, G and T.
    fn writes_hgvs(&self, alternative: &str) -> bool {
        !self.config.vep_compatibility_mode || is_unambiguous(alternative)
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
        // Run each input variant through single-variant normalization first. Like in `predict`,
        // a VCF indel then lies behind its first REF base, as the VCF specification defines.
        let normalized_variants: Vec<VcfVariant> =
            variants.iter().map(|v| self.normalize_variant(v)).collect();

        // A variant without change changes nothing in the group.
        let mut paired_vars: Vec<_> = (variants.iter().cloned().zip(normalized_variants))
            .filter(|(_, norm)| !(norm.reference.is_empty() && norm.alternative.is_empty()))
            .collect();
        // An insertion lies in front of the base at its position, so it sorts before a variant
        // that starts at that base.
        paired_vars.sort_by_key(|(_, norm)| (norm.position, norm.reference.len()));

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
            // As in `predict`, `hgvs_g` is 3'-shifted on the strand of the transcript.
            let (fwd, rev) = self.normalize_g(&compound_var_g, chrom_acc)?;
            hgvs_g = Some((self.format_hgvs_g(&fwd), self.format_hgvs_g(&rev)));
        }

        let txs = self.get_transcripts_for_variant_group(&sorted_vars, chrom_acc)?;
        if txs.is_empty() {
            return Ok(None);
        }

        let multi_anns = self.select_annotations(txs, |tx_record| {
            let hgvs_g = hgvs_g.as_ref().map(|(fwd, rev)| {
                if tx_record.alt_strand == -1 {
                    rev.clone()
                } else {
                    fwd.clone()
                }
            });
            self.predict_multiple_for_transcript(
                &sorted_vars,
                &sorted_originals,
                tx_record,
                chrom_acc,
                hgvs_g,
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
        // As in `predict_multiple`, an insertion sorts before a variant at the same position.
        sorted_vars.sort_by_key(|v| (v.position, v.reference.len()));

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
            let (tx_loc, tx_csqs) = self
                .determine_transcript_context(alignment, strand, &var_g, var_start, var_end, false);

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

            // The n. edit stays at the input position. Normalized on its own, it could move
            // past another variant of the group. `predict_multiple_for_transcript` normalizes
            // the combined change.
            let mapper = self.mapper.variant_mapper();
            let var_n = match mapper.g_to_n(&var_g, &tx.id, ALT_ALN_METHOD) {
                Ok(var_n) => var_n,
                Err(Error::NonAdjacentExons(_, _, _, _)) => {
                    tracing::warn!("{}, {}: NonAdjacentExons, skipping", &tx.id, var_g);
                    return Ok(None);
                }
                Err(e) => return Err(SeqvarsError::from(e)),
            };
            let var_c = mapper
                .n_to_c(&var_n)
                .map_err(|e| SeqvarsError::HgvsProjection(e.to_string()))?;
            let proj = HgvsProjectionContext {
                n: Some(var_n),
                c: Some(var_c),
                p: None,
                region: None,
                unshifted_g: None,
            };
            if !proj.is_within_cds_bounds() {
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
            match self.get_group_projections(sorted_vars, tx, chrom_acc)? {
                Some((p, csqs)) => (p, csqs),
                None => return Ok(None),
            };

        if projections.is_empty() {
            // All variants in this group were purely intronic and thus skipped.
            return Ok(None);
        }

        let Ok(ref_data) = self.ref_data(&tx.id) else {
            return Ok(None);
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

        // Like `predict`, normalize the n. variant and project it to c.
        let compound_var_n = self.mapper.maybe_normalize(&compound_var_n)?;
        let compound_var_c = match self.mapper.variant_mapper().n_to_c(&compound_var_n) {
            Ok(c) => c,
            Err(e) => {
                tracing::debug!("Failed to map compound n_loc to c_loc: {}", e);
                return Ok(None);
            }
        };

        let (compound_var_p, region) = self.safe_project_c_to_p(&compound_var_c)?.unzip();
        let region = region.flatten();

        let mut compound_proj = HgvsProjectionContext {
            n: Some(compound_var_n),
            c: Some(compound_var_c),
            p: compound_var_p,
            region,
            unshifted_g: None,
        };
        self.mark_unknown_protein(&mut compound_proj, &tx.id);

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
            None,
            None,
            tx_len,
            transcript_biotype,
        )?;

        let mut consequences = c_ctx.cds_consequences | base_group_consequences;
        self.consequences_fix_special_cases(
            &mut consequences,
            c_ctx.cds_consequences,
            c_ctx.protein_consequences,
            c_ctx.sequence_consequences,
            c_ctx.vep_codon_consequences,
            &compound_proj,
            tx,
        );

        if self.vep_consequence_terms() {
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

        // As in `build_ann_field`, see `writes_hgvs`.
        let hgvs = (sorted_originals.iter()).all(|var| self.writes_hgvs(&var.alternative));
        let [hgvs_g, hgvs_n, hgvs_c, hgvs_p] =
            [hgvs_g, hgvs_n, hgvs_c, hgvs_p].map(|desc| desc.filter(|_| hgvs));

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

        let ins_at_boundary_is_intronic = (self.edge_bases_in_intron(&compound_proj))
            .unwrap_or(self.mapper.config.ins_at_boundary_is_intronic);
        Ok(Some(TxAnnotation {
            context: self
                .tx_context(Some(tx.id.clone()), compound_proj.c, compound_proj.n)
                .with_ref_data(ref_data)
                .with_alt_inputs(None, false, ins_at_boundary_is_intronic),
            ann,
        }))
    }

    /// Whether the c. variant of `projection` leaves no complete codon of the CDS behind it, so
    /// that none reads in another frame. A change of the length from the CDS through its last
    /// base `c.N` does so, see `crosses_cds_end`. Outside VEP compatibility mode, so does a
    /// frameshift that ends inside the stop codon, `c.(N-2)` to `c.N`: an insertion or the copy
    /// of a dup that goes inside it, or a deletion or delins from the CDS that ends in it. An
    /// in-frame change there keeps the frame, and its p. terms stay. `false` for a CDS with an
    /// incomplete 3' end, as it has no stop codon.
    fn leaves_no_codon_behind(&self, projection: &HgvsProjectionContext, tx: &Transcript) -> bool {
        let Some(var_c @ HgvsVariant::CdsVariant { loc_edit, .. }) = &projection.c else {
            return false;
        };
        let Some(cds_len) = (tx.start_codon.zip(tx.stop_codon))
            .filter(|_| !tx.is_incomplete_3p())
            .map(|(start, stop)| stop - start)
        else {
            return false;
        };
        let loc = loc_edit.loc.inner();
        let edit = loc_edit.edit.inner();
        if changes_length(edit) && crosses_cds_end(var_c, Some(cds_len)) {
            return true;
        }
        if self.config.vep_compatibility_mode
            || loc.start.cds_from != CdsFrom::Start
            || loc.start.base < 1
            || loc.start.offset.unwrap_or_default() != 0
        {
            return false;
        }
        let in_stop_codon = |pos: &CdsPos| {
            pos.cds_from == CdsFrom::Start
                && pos.base >= cds_len - 2
                && pos.offset.unwrap_or_default() == 0
        };
        // The CDS bases from `start` to `end`.
        let cds_bases = usize::try_from(loc.end.base - loc.start.base + 1).unwrap_or_default();
        // The bases go behind `start`, the copy behind `end`.
        let (inside, inserted, deleted) = match edit {
            NaEdit::Ins { alternative } => (
                in_stop_codon(&loc.start) && loc.start.base < cds_len,
                alternative.len(),
                0,
            ),
            NaEdit::Dup { .. } => (
                in_stop_codon(&loc.end) && loc.end.base < cds_len,
                cds_bases,
                0,
            ),
            NaEdit::DelRef { .. } | NaEdit::DelNum { .. } => {
                (in_stop_codon(&loc.end), 0, cds_bases)
            }
            NaEdit::RefAlt { alternative, .. } | NaEdit::NumAlt { alternative, .. } => {
                (in_stop_codon(&loc.end), alternative.len(), cds_bases)
            }
            NaEdit::InvRef { .. } | NaEdit::InvNum { .. } => (false, 0, 0),
        };
        inside && inserted.abs_diff(deleted) % 3 != 0
    }

    /// Terms of the protein that `var_c` produces, compared to the reference protein of
    /// `ref_data`. The first stop codon of the new protein decides:
    /// - in front of the reference stop codon: `stop_gained`;
    /// - at its position: `missense_variant` if residues change, and `stop_retained_variant` if
    ///   the variant changes the bases there;
    /// - behind it: `stop_lost` and `feature_elongation`;
    /// - none: `stop_lost`, with `feature_elongation` if the protein gets longer.
    ///
    /// Without a reference stop codon, the reference end takes its place, with neither
    /// `stop_lost` nor `stop_retained_variant`. A new protein without a stop codon then gets no
    /// term, as the residues behind that end are unknown.
    ///
    /// Empty if the variant changes no base up to the reference stop codon, e.g., a deletion that
    /// can go behind it, and without a protein, as for a UTR variant.
    fn sequence_consequences(
        &self,
        var_c: &HgvsVariant,
        ref_data: &RefTranscriptData,
    ) -> Consequences {
        // `build_altseq` builds no protein outside the CDS or for a deletion of the whole CDS.
        let Some(alt_data) = self
            .mapper
            .alt_seq_builder(var_c, ref_data)
            .and_then(|builder| builder.build_altseq())
            .ok()
            .and_then(|alt_data| alt_data.into_iter().next())
            .filter(|alt_data| alt_data.variant_start_aa.is_some())
        else {
            return Consequences::empty();
        };
        let ref_seq = ref_data.transcript_sequence.as_bytes();
        let alt_seq = alt_data.transcript_sequence.as_bytes();
        let stop_end = usize::try_from(ref_data.cds_stop).unwrap_or_default();
        let prefix = ref_seq
            .iter()
            .zip(alt_seq)
            .take_while(|(r, a)| r == a)
            .count();
        if prefix >= stop_end {
            return Consequences::empty();
        }
        // Behind a change of the length, the bases at the stop codon are other ones.
        let changes_stop_codon = ref_seq.len() != alt_seq.len() || {
            let suffix = (ref_seq.iter().rev().zip(alt_seq.iter().rev()))
                .take_while(|(r, a)| r == a)
                .count();
            ref_seq.len() - suffix > stop_end.saturating_sub(3)
        };

        // Count the residues before the stops. `Xaa` counts as a residue.
        let reference = ref_data.aa_sequence.trim_end_matches('*');
        let has_stop = reference.len() < ref_data.aa_sequence.len();
        let altered = alt_data.aa_sequence.as_str();
        match altered.find('*') {
            Some(pos) => match pos.cmp(&reference.len()) {
                Ordering::Less => Consequence::StopGained.into(),
                Ordering::Equal => {
                    let mut consequences = Consequences::empty();
                    if altered[..pos] != *reference {
                        consequences |= Consequence::MissenseVariant;
                    }
                    if has_stop && changes_stop_codon {
                        consequences |= Consequence::StopRetainedVariant;
                    }
                    consequences
                }
                Ordering::Greater if has_stop => {
                    Consequence::StopLost | Consequence::FeatureElongation
                }
                Ordering::Greater => Consequence::FeatureElongation.into(),
            },
            None if !has_stop => Consequences::empty(),
            None if altered.len() > reference.len() => {
                Consequence::StopLost | Consequence::FeatureElongation
            }
            None => Consequence::StopLost.into(),
        }
    }

    /// Make p. of `projection` on transcript `tx_ac` `p.?` if the protein reads a codon with an
    /// unknown base, which can be any amino acid or a stop. VEP compatibility mode keeps p., as
    /// VEP's rules for such bases decide the terms there, see `adjust_vep_terms`.
    fn mark_unknown_protein(&self, projection: &mut HgvsProjectionContext, tx_ac: &str) {
        if self.config.vep_compatibility_mode
            || projection.p.is_none()
            || !projection.has_ambiguous_alt()
        {
            return;
        }
        if let Some(var_c) = &projection.c
            && self.reads_unknown_codon(var_c, tx_ac)
        {
            projection.p = Some(HgvsVariant::ProtVariant {
                accession: var_c.accession().clone(),
                gene_symbol: var_c.gene_symbol().clone(),
                loc_edit: ProtLocEdit::Unknown,
            });
        }
    }

    /// Whether the protein of `var_c` on transcript `tx_ac` reads a codon with an unknown base,
    /// e.g. `N`, from the first changed residue up to a stop. hgvs-rs translates such a codon to
    /// `X` unless all its readings give the same amino acid.
    fn reads_unknown_codon(&self, var_c: &HgvsVariant, tx_ac: &str) -> bool {
        let Ok(ref_data) = self.ref_data(tx_ac) else {
            return false;
        };
        let Some((alt_data, start)) = self
            .mapper
            .alt_seq_builder(var_c, &ref_data)
            .and_then(|builder| builder.build_altseq())
            .ok()
            .and_then(|alt_data| alt_data.into_iter().next())
            .and_then(|alt_data| {
                let start = usize::try_from(alt_data.variant_start_aa?).ok()?;
                Some((alt_data, start))
            })
        else {
            return false;
        };
        (alt_data.aa_sequence.get(start.saturating_sub(1)..))
            .and_then(|changed| changed.split('*').next())
            .is_some_and(|read| read.contains('X'))
    }

    /// The reference data of transcript `tx_ac`, shared with the cache of hgvs-rs.
    fn ref_data(&self, tx_ac: &str) -> Result<Arc<RefTranscriptData>, Error> {
        let provider: Arc<dyn Provider + Send + Sync> = self.provider.clone();
        ref_transcript_data_shared(&provider, tx_ac, None)
    }

    /// The region of the transcript that `var_c` changes, as `c_to_p` of `self.mapper` uses it.
    fn variant_location(&self, var_c: &HgvsVariant) -> Option<VariantLocation> {
        let ref_data = self.ref_data(&var_c.accession().value).ok()?;
        self.mapper
            .alt_seq_builder(var_c, &ref_data)
            .and_then(|builder| builder.variant_location())
            .ok()
    }

    /// Whether `c_to_p` puts the bases of an insertion at an exon edge, or of the copy of a dup
    /// that can go there, into the intron, as `projection.region` tells. `None` if hgvs-rs does
    /// not decide this by sequence, as in VEP compatibility mode, and in the UTRs, where the
    /// region is that of the UTR.
    fn edge_bases_in_intron(&self, projection: &HgvsProjectionContext) -> Option<bool> {
        let Some(HgvsVariant::CdsVariant { loc_edit, .. }) = &projection.c else {
            return None;
        };
        let end = &loc_edit.loc.inner().end;
        if !self.mapper.config.ins_at_boundary_by_sequence
            || end.cds_from == CdsFrom::End
            || end.base < 0
        {
            return None;
        }
        projection
            .region
            .map(|region| region == VariantLocation::Intron)
    }

    /// Safely projects a CDS variant to a Protein variant, gracefully catching
    /// and swallowing expected incomplete-transcript errors as `None`. With the protein comes the
    /// region of the transcript that `c_to_p` built it for. For a dup whose intron bases are
    /// unknown (`DuplicationReferenceMissing`), `c_to_p` builds nothing: the protein is `p.?` and
    /// the region comes from the builder, see `variant_location`.
    fn safe_project_c_to_p(
        &self,
        var_c: &HgvsVariant,
    ) -> Result<Option<(HgvsVariant, Option<VariantLocation>)>, SeqvarsError> {
        self.mapper.c_to_p_with_region(var_c).map_or_else(
            |e| {
                if matches!(e, Error::DuplicationReferenceMissing(_)) {
                    // The intron bases of the dup are unknown, so is its protein change.
                    let var_p = HgvsVariant::ProtVariant {
                        accession: var_c.accession().clone(),
                        gene_symbol: var_c.gene_symbol().clone(),
                        loc_edit: ProtLocEdit::Unknown,
                    };
                    return Ok(Some((var_p, self.variant_location(var_c))));
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
            |(var_p, region)| Ok(Some((var_p, Some(region)))),
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

/// Whether `var_c` changes a base of the start codon c.1 to c.3: an insertion between two of
/// them, a dup whose copy goes there, behind its last base, or another edit of one of them.
fn changes_start_codon(var_c: &HgvsVariant) -> bool {
    let HgvsVariant::CdsVariant { loc_edit, .. } = var_c else {
        return false;
    };
    let loc = loc_edit.loc.inner();
    let in_exon = |pos: &CdsPos| pos.cds_from == CdsFrom::Start && pos.offset.unwrap_or(0) == 0;
    if !(in_exon(&loc.start) && in_exon(&loc.end)) {
        return false;
    }
    match loc_edit.edit.inner() {
        NaEdit::Ins { .. } => {
            (1..=2).contains(&loc.start.base) && loc.end.base == loc.start.base + 1
        }
        NaEdit::Dup { .. } => (1..=2).contains(&loc.end.base),
        _ => loc.start.base <= 3 && loc.end.base >= 1,
    }
}

/// The edit of the p. variant of `projection`, if it has one.
fn protein_edit(projection: &HgvsProjectionContext) -> Option<&ProteinEdit> {
    match &projection.p {
        Some(HgvsVariant::ProtVariant {
            loc_edit: ProtLocEdit::Ordinary { edit, .. },
            ..
        }) => Some(edit.inner()),
        _ => None,
    }
}

/// Whether the insertion or dup `var_c` keeps the frame and puts its bases in front of the stop
/// codon of a CDS of length `cds_len`, so that the stop codon stays.
fn inserts_in_frame_before_stop_codon(var_c: &HgvsVariant, cds_len: i32) -> bool {
    let HgvsVariant::CdsVariant { loc_edit, .. } = var_c else {
        return false;
    };
    let loc = loc_edit.loc.inner();
    if [&loc.start, &loc.end]
        .iter()
        .any(|pos| pos.cds_from != CdsFrom::Start || pos.offset.unwrap_or_default() != 0)
    {
        return false;
    }
    // The bases go behind `start`, the copy behind `end`.
    let (behind, len) = match loc_edit.edit.inner() {
        NaEdit::Ins { alternative } => (loc.start.base, alternative.len()),
        NaEdit::Dup { .. } => (
            loc.end.base,
            usize::try_from(loc.end.base - loc.start.base + 1).unwrap_or_default(),
        ),
        _ => return false,
    };
    behind <= cds_len - 3 && len % 3 == 0
}

/// Whether `var_c` changes bases of the stop codon of a CDS of length `cds_len`, i.e., starts in
/// the CDS and reaches its last 3 bases or behind.
fn touches_stop_codon(var_c: &HgvsVariant, cds_len: Option<i32>) -> bool {
    let (HgvsVariant::CdsVariant { loc_edit, .. }, Some(cds_len)) = (var_c, cds_len) else {
        return false;
    };
    let loc = loc_edit.loc.inner();
    let in_cds =
        |pos: &CdsPos| pos.cds_from == CdsFrom::Start && pos.offset.unwrap_or_default() == 0;
    in_cds(&loc.start)
        && loc.start.base >= 1
        && (loc.end.cds_from == CdsFrom::End || (in_cds(&loc.end) && loc.end.base > cds_len - 3))
}

/// VEP's `_overlaps_stop_codon_cil` (VEP 116): whether `var_g`, a variant at its input position,
/// overlaps the last 3 genomic bases of the CDS of `alignment`. Where an intron splits the stop
/// codon, these hold intron bases. VEP extends an insertion by the number of its bases, to the
/// right on the plus strand and to the left on the minus strand.
fn vep_overlaps_stop_codon(var_g: &HgvsVariant, alignment: &GenomeAlignment) -> bool {
    let HgvsVariant::GenomeVariant { loc_edit, .. } = var_g else {
        return false;
    };
    let loc = loc_edit.loc.inner();
    let (Some(start), Some(end), Some(cds_start), Some(cds_end)) =
        (loc.start, loc.end, alignment.cds_start, alignment.cds_end)
    else {
        return false;
    };
    let minus = Strand::try_from(alignment.strand) == Ok(Strand::Minus);
    // 1-based, from VEP's `coding_region_start` or to its `coding_region_end`.
    let (from, to) = if minus {
        (cds_start + 1, cds_start + 3)
    } else {
        (cds_end - 2, cds_end)
    };
    // VEP puts an insertion between `end`, the base in front, and `start = end + 1`.
    let (start, end) = match loc_edit.edit.inner() {
        NaEdit::Ins { alternative } => {
            let len = alternative.len() as i32;
            (if minus { end - len } else { end + len }, start)
        }
        _ => (start, end),
    };
    end >= from && start <= to
}

/// Whether `bases` hold only A, C, G and T in either case, as VEP's `seq_is_unambiguous_dna`
/// checks.
fn is_unambiguous(bases: &str) -> bool {
    bases
        .bytes()
        .all(|base| b"ACGT".contains(&base.to_ascii_uppercase()))
}

/// Whether the amino acids `s` of a p. edit hold a stop. hgvs-rs gives `X` (`Xaa`) for an
/// incomplete codon at the transcript end, not for a stop codon.
fn has_stop(s: &str) -> bool {
    s.contains('*') || s.contains("Ter")
}

/// Whether `edit` deletes bases or replaces them by a different number of bases. An edit
/// without its reference bases counts as such.
fn changes_length(edit: &NaEdit) -> bool {
    match edit {
        NaEdit::DelRef { .. } | NaEdit::DelNum { .. } => true,
        NaEdit::RefAlt {
            reference,
            alternative,
        } => reference.is_empty() || reference.len() != alternative.len(),
        NaEdit::NumAlt { count, alternative } => {
            usize::try_from(*count).ok() != Some(alternative.len())
        }
        _ => false,
    }
}

/// Whether `var_c` lies in one UTR: both ends in front of c.1, or both behind the CDS.
fn lies_in_one_utr(var_c: &HgvsVariant) -> bool {
    let HgvsVariant::CdsVariant { loc_edit, .. } = var_c else {
        return false;
    };
    let loc = loc_edit.loc.inner();
    let in_5p_utr = |pos: &CdsPos| pos.cds_from == CdsFrom::Start && pos.base < 1;
    (in_5p_utr(&loc.start) && in_5p_utr(&loc.end))
        || (loc.start.cds_from == CdsFrom::End && loc.end.cds_from == CdsFrom::End)
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
            // An insertion lies between two CDS bases, also at an exon edge, e.g. `c.9_9+1ins`
            // or `c.10-1_10ins` behind c.9. It is conservative if the first base ends a codon.
            if loc_edit.edit.inner().is_ins() {
                let in_front = if loc.start.offset.unwrap_or(0) < 0 {
                    end_base - 1
                } else {
                    start_base
                };
                return start_cds_from == CdsFrom::Start
                    && end_cds_from == CdsFrom::Start
                    && in_front % 3 == 0;
            }
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

/// Whether `var_c` starts in the CDS, not in an intron, and ends at its last base, `c.N` for a
/// CDS of length `cds_len`, or behind it, e.g. at `c.*N`.
fn crosses_cds_end(var_c: &HgvsVariant, cds_len: Option<i32>) -> bool {
    let HgvsVariant::CdsVariant { loc_edit, .. } = var_c else {
        return false;
    };
    let loc = loc_edit.loc.inner();
    let ends_at_cds_end = loc.end.cds_from == CdsFrom::End
        || (loc.end.cds_from == CdsFrom::Start && cds_len.is_some_and(|len| loc.end.base >= len));
    loc.start.cds_from == CdsFrom::Start
        && loc.start.base >= 1
        && loc.start.offset.unwrap_or_default() == 0
        && ends_at_cds_end
}

/// Whether deleting `tx_seq[del]` leaves the CDS that starts at `cds_start` intact.
///
/// `del` must end behind `cds_start`.  The CDS stays intact if it then starts
/// `del.len()` bases earlier, i.e. if the same bases could be deleted from the 5' UTR.
/// VEP's `_ins_del_start_altered` checks the same.
///
/// Only VEP compatibility mode uses this check. Outside it, the region from hgvs-rs decides.
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

/// The runs of bases that `edit`, starting at the 0-based position `start`, changes, as 0-based,
/// end-exclusive ranges. `None` unless `edit` replaces bases by as many others and keeps some of
/// them, e.g. the C of `GCT>ACC` or of the inversion `GCAT>ATGC`.
fn changed_runs(edit: &NaEdit, start: i32) -> Option<Vec<(i32, i32)>> {
    let (reference, alternative) = match edit {
        NaEdit::RefAlt {
            reference,
            alternative,
        } if reference.len() == alternative.len() => (reference, alternative.clone()),
        NaEdit::InvRef { reference } => (reference, revcomp(reference)),
        _ => return None,
    };
    let mut runs: Vec<(i32, i32)> = Vec::new();
    let bases = reference.bytes().zip(alternative.bytes());
    for (pos, (ref_base, alt_base)) in (start..).zip(bases) {
        if ref_base.eq_ignore_ascii_case(&alt_base) {
            continue;
        }
        match runs.last_mut() {
            Some((_, end)) if *end == pos => *end = pos + 1,
            _ => runs.push((pos, pos + 1)),
        }
    }
    let len = i32::try_from(reference.len()).ok()?;
    (runs != [(start, start + len)]).then_some(runs)
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
    use crate::annotate::seqvars::consequence::plugin::{
        NmdPlugin, SequencePlugin, TranscriptPlugin,
    };
    use crate::annotate::seqvars::consequence::terms::{
        ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_NMD, ANN_TX_SEQ_ALT, ANN_TX_SEQ_REF,
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

    /// Without a reference genome, mehari warns that HGVS descriptions may not be shifted and
    /// that reference alleles come from the transcripts. It warns once per assembly, also if two
    /// predictors share it, as in the server.
    #[tracing_test::traced_test]
    #[test]
    fn warn_without_reference() {
        let provider = Arc::new(MehariProvider::new(
            tx_db_with_indels_in_exon_cigars(),
            None::<PathBuf>,
            true,
        ));
        ConsequencePredictor::new(provider.clone(), Default::default());
        ConsequencePredictor::new(provider, Default::default());
        logs_assert(|lines: &[&str]| {
            let warnings = lines
                .iter()
                .filter(|line| line.contains("No reference genome given"))
                .collect::<Vec<_>>();
            match warnings[..] {
                [warning]
                    if warning.contains("No reference genome given for grch37:")
                        && warning.contains("transcripts that differ from the genome") =>
                {
                    Ok(())
                }
                _ => Err(format!("{warnings:?}")),
            }
        });
    }

    /// HGVS writes a predicted protein change in parentheses. c.4G>A changes GCC (Ala2) to ACC.
    #[rstest::rstest]
    #[case::predicted(HgvsParentheses::Predicted, "p.(Ala2Thr)")]
    #[case::never(HgvsParentheses::Never, "p.Ala2Thr")]
    fn predicted_protein_in_parentheses(
        #[case] hgvs_parentheses: HgvsParentheses,
        #[case] expected: &str,
    ) -> Result<(), anyhow::Error> {
        let provider = Arc::new(MehariProvider::new(
            tx_db_with_indels_in_exon_cigars(),
            None::<PathBuf>,
            true,
        ));
        let config = ConfigBuilder::default()
            .hgvs_parentheses(hgvs_parentheses)
            .build()?;
        let anns = ConsequencePredictor::new(provider, config)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1006,
                reference: "G".into(),
                alternative: "A".into(),
            })?
            .unwrap();

        assert_eq!(
            (anns[0].hgvs_c.as_deref(), anns[0].hgvs_p.as_deref()),
            (Some("c.4G>A"), Some(expected))
        );
        Ok(())
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
        assert_eq!(ann[0].hgvs_p.as_deref(), Some("p.(Lys8Ter)"));
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

    /// A deletion from the end of an exon through the intron into the next exon deletes the CDS
    /// bases on both sides of the intron. On OPA1 `NM_130837.3`, `c.1375_1380del` deletes 89
    /// genome bases, but 6 CDS bases, so it keeps the frame. The test database has no genome
    /// sequence, so the intron bases of the VCF REF allele are made up, as in
    /// `alt_sequences_deletion_covers_intron`.
    #[rstest::rstest]
    #[case(
        false,
        vec![
            Consequence::SpliceAcceptorVariant,
            Consequence::SpliceDonorVariant,
            Consequence::ConservativeInframeDeletion,
            Consequence::SpliceDonorFifthBaseVariant,
            Consequence::SpliceRegionVariant,
            Consequence::ExonicSpliceRegionVariant,
            Consequence::SpliceDonorRegionVariant,
            Consequence::SplicePolypyrimidineTractVariant
        ]
    )]
    // VEP counts the CDS bases for `frameshift_variant`, too, and builds codons for a deletion
    // between two exon positions
    #[case(
        true,
        vec![
            Consequence::SpliceAcceptorVariant,
            Consequence::SpliceDonorVariant,
            Consequence::InframeDeletion,
            Consequence::SpliceDonorFifthBaseVariant
        ]
    )]
    fn deletion_across_intron_keeps_frame(
        #[case] vep_compatibility_mode: bool,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .transcript_settings(pick_mane_by_gene())
                .vep_consequence_terms(vep_compatibility_mode)
                .vep_compatibility_mode(vep_compatibility_mode)
                .build()?,
        );

        let res = predictor
            .predict(&VcfVariant {
                chromosome: "3".into(),
                position: 193361230,
                reference: format!("TAAT{}ACT", "A".repeat(83)),
                alternative: "T".into(),
            })?
            .unwrap_or_default();
        let ann = res
            .iter()
            .find(|ann| ann.feature_id == "NM_130837.3")
            .ok_or_else(|| anyhow::anyhow!("no annotation on NM_130837.3: {res:?}"))?;
        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (Some("c.1375_1380del"), &expected_csqs)
        );

        Ok(())
    }

    /// With VEP compatibility, the terms of ClinVar variants on BRCA1 `NM_007294.4` as VEP 109
    /// reports them in `clinvar.excerpt.vep.brca1.tsv`. VEP 116 has the same predicates for them.
    #[rstest::rstest]
    // VEP builds codons for an exonic change in the splice region, so it reports no
    // `coding_sequence_variant`.
    #[case(
        "17:41197817:T:A",
        "c.5470A>T",
        "missense_variant&splice_region_variant"
    )]
    #[case(
        "17:41197818:T:C",
        "c.5469A>G",
        "splice_region_variant&synonymous_variant"
    )]
    #[case("17:41201209:G:A", "c.5335C>T", "splice_region_variant&stop_gained")]
    #[case(
        "17:41201137:CT:C",
        "c.5406del",
        "frameshift_variant&splice_region_variant"
    )]
    #[case(
        "17:41242960:CCTG:C",
        "c.4183_4185del",
        "inframe_deletion&splice_region_variant"
    )]
    // VEP puts an insertion at an exon edge into the exon.
    #[case(
        "17:41276033:C:CCAGATG",
        "c.80_80+1insCATCTG",
        "inframe_insertion&splice_region_variant"
    )]
    // VEP builds no codons for a deletion with an end in an intron.
    #[case(
        "17:41267741:ACTT:A",
        "c.133_134+1del",
        "coding_sequence_variant&splice_donor_variant"
    )]
    fn vep_terms_as_recorded(
        #[case] spdi: &str,
        #[case] expected_hgvs_c: &str,
        #[case] expected_csqs: &str,
    ) -> Result<(), anyhow::Error> {
        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .vep_consequence_terms(true)
                .vep_compatibility_mode(true)
                .build()?,
        );

        let res = predictor.predict(&spdi.parse()?)?.unwrap_or_default();
        let ann = res
            .iter()
            .find(|ann| ann.feature_id == "NM_007294.4")
            .ok_or_else(|| anyhow::anyhow!("no annotation on NM_007294.4: {res:?}"))?;
        let csqs = ann
            .consequences
            .iter()
            .map(|csq| csq.to_string())
            .sorted()
            .join("&");
        assert_eq!(
            (ann.hgvs_c.as_deref(), csqs.as_str()),
            (Some(expected_hgvs_c), expected_csqs)
        );

        Ok(())
    }

    /// Indels on Ensembl 108 chr22 transcripts.
    ///
    /// With VEP terms, the expected terms are those of VEP 116.3 with `--shift_3prime 1`. mehari
    /// adds `protein_altering_variant`. At the stop codon, mehari keeps `frameshift_variant` next
    /// to `stop_lost`, which VEP drops.
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
    #[case("22:22895417:CCT:C", "ENST00000531372", true, vec![Consequence::FrameshiftVariant, Consequence::StopLost])]
    // `p.Ter204TrpextTer73`
    #[case("22:42571197:T:TG", "ENST00000340239", true, vec![Consequence::FrameshiftVariant, Consequence::StopLost])]
    // `p.Ter704=`: VEP's changed codons end in a partial codon, and the bases at the stop codon
    // still read a stop codon
    #[case("22:45600444:TG:T", "ENST00000327858", true, vec![Consequence::StopRetainedVariant])]
    // `p.Ter204=`: VEP's changed codons end in a partial codon, and VEP does not count the
    // insertion behind the first base of the stop codon as overlapping it
    #[case("22:42571196:T:TTAAC", "ENST00000340239", true, vec![Consequence::FrameshiftVariant])]
    // `p.Ter512=`: the new frame starts with a stop codon
    #[case("22:17181484:C:CT", "ENST00000262607", true, vec![Consequence::StopRetainedVariant])]
    // `p.Ter704AspextTer1`: in-frame insertion in front of the stop codon
    #[case("22:45600443:C:CGAT", "ENST00000327858", true, vec![Consequence::InframeInsertion])]
    // `p.Cys203Ter` (VEP: `p.Cys203del`): in-frame deletion in front of the stop codon
    #[case("22:42571193:CTGT:C", "ENST00000340239", true, vec![Consequence::InframeDeletion])]
    // in-frame deletion across the stop codon
    #[case("22:22895418:CTCT:C", "ENST00000531372", true, vec![Consequence::StopLost, Consequence::InframeDeletion, Consequence::ProteinAlteringVariant])]
    // `c.*2_*4del`: VEP checks no stop codon for a deletion behind the CDS, also if the input
    // position overlaps the stop codon
    #[case("22:22895421:TGAG:T", "ENST00000531372", true, vec![Consequence::ThreePrimeUtrVariant])]
    // `c.1449_*18del`: VEP builds no altered sequence for a deletion beyond the transcript end,
    // so the stop codon counts as kept
    #[case("22:50525751:GTTTCGCGGCAAAGGAGCTT:G", "ENST00000487577", true, vec![Consequence::StopRetainedVariant, Consequence::ThreePrimeUtrVariant])]
    // `p.Ter704=`: in-frame insertion inside the stop codon that keeps it
    #[case("22:45600443:C:CTAA", "ENST00000327858", false, vec![Consequence::StopRetainedVariant])]
    #[case("22:45600443:C:CTAA", "ENST00000327858", true, vec![Consequence::InframeInsertion, Consequence::StopRetainedVariant])]
    // `p.Phe703=`: in-frame insertion that keeps codon 703 and forms a new stop codon in front
    // of the old one. VEP's changed codons translate to `F` and `F*H`, so VEP reports no
    // `synonymous_variant`.
    #[case("22:45600442:T:TTTAGCA", "ENST00000327858", false, vec![Consequence::SynonymousVariant])]
    #[case("22:45600442:T:TTTAGCA", "ENST00000327858", true, vec![Consequence::InframeInsertion, Consequence::StopRetainedVariant])]
    // The CDS of ENST00000430886 ends with a complete codon, so none of its cases gets
    // `incomplete_terminal_codon_variant`.
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
                .vep_compatibility_mode(vep_consequence_terms)
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

    /// Insertions at the end of a CDS without a stop codon (`cds_end_NF`), on transcripts of
    /// `CHR22_CDS_END_TXS`. With VEP terms, the expected stop terms are those of VEP 116.3 with
    /// `--shift_3prime 1`. Inside a partial last codon, VEP reports `coding_sequence_variant`
    /// instead of `frameshift_variant` or `inframe_insertion`, and mehari keeps those.
    #[rstest::rstest]
    // ENST00000452294 (minus strand) ends with the partial codon `tg`. `X/*` (`tg/tAg`)
    #[case("22:38750888:C:CT", "ENST00000452294", true, "c.433_434insA", vec![Consequence::StopGained, Consequence::FrameshiftVariant, Consequence::IncompleteTerminalCodonVariant])]
    // `X/*` (`tg/tAATg`)
    #[case("22:38750888:C:CATT", "ENST00000452294", true, "c.433_434insAAT", vec![Consequence::StopGained, Consequence::InframeInsertion, Consequence::IncompleteTerminalCodonVariant])]
    // `X/*L` (`tg/tAACTg`)
    #[case("22:38750888:C:CAGTT", "ENST00000452294", true, "c.433_434insAACT", vec![Consequence::StopGained, Consequence::FrameshiftVariant, Consequence::IncompleteTerminalCodonVariant])]
    #[case("22:38750888:C:CT", "ENST00000452294", false, "c.433_434insA", vec![Consequence::FrameshiftVariant, Consequence::IncompleteTerminalCodonVariant])]
    // ENST00000578984 (minus strand) ends with the partial codon `a`. `-/*` (`-/TAG`): the
    // insertion lies between the last complete codon and the partial codon.
    #[case("22:17819036:T:TCTA", "ENST00000578984", true, "c.567_568insTAG", vec![Consequence::StopGained, Consequence::InframeInsertion, Consequence::IncompleteTerminalCodonVariant])]
    fn annotate_incomplete_cds_end_csqs(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] vep_consequence_terms: bool,
        #[case] expected_hgvs_c: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let spdi = spdi.split(':').collect::<Vec<_>>();

        let tx_db = load_tx_db(CHR22_CDS_END_TXS)?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .vep_consequence_terms(vep_consequence_terms)
                .vep_compatibility_mode(vep_consequence_terms)
                .build()?,
        );

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].to_string(),
                position: spdi[1].parse()?,
                reference: spdi[2].to_string(),
                alternative: spdi[3].to_string(),
            })?
            .unwrap_or_default();
        let ann = res
            .iter()
            .find(|ann| ann.feature_id.starts_with(tx_id))
            .ok_or_else(|| anyhow::anyhow!("no annotation on {}", tx_id))?;
        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (Some(expected_hgvs_c), &expected_csqs),
            "spdi = {}",
            spdi.join(":")
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

    /// Ensembl 108 chr22 transcripts for `annotate_chr22_window`.
    const CHR22_ENSEMBL108_TXS: &str =
        "tests/data/annotate/db/grch38/GRCh38-ensembl.exon-edge-subset.txs.bin.zst";

    /// GENCODE 42 chr22 transcripts that share exon edges with those of `CHR22_ENSEMBL108_TXS`,
    /// and those of `annotate_dup_csqs`, `annotate_exon_edge_ins_stop_csqs` and
    /// `annotate_vep_ins_region_csqs`.
    const CHR22_GENCODE42_TXS: &str =
        "tests/data/annotate/db/grch38/GRCh38-gencode42.exon-edge-subset.txs.bin.zst";

    /// GENCODE 42 chr22 transcripts from `db create` with the GENCODE 42 chr22 GFF3. The CDS of
    /// a transcript without a stop codon runs to the transcript end, as in VEP.
    const CHR22_CDS_END_TXS: &str =
        "tests/data/annotate/db/grch38/GRCh38-gencode42.cds-end.txs.bin.zst";

    /// Annotation of `spdi` on the transcript `tx_id`, with Ensembl 108 chr22 transcripts and
    /// the chr22 reference around the variant, so that mehari shifts indels 3'.
    fn annotate_chr22_window(
        spdi: &str,
        tx_id: &str,
        vep_consequence_terms: bool,
    ) -> Result<AnnField, anyhow::Error> {
        let config = ConfigBuilder::default()
            .vep_consequence_terms(vep_consequence_terms)
            .vep_compatibility_mode(vep_consequence_terms)
            .build()?;
        annotate_chr22_window_with_config(spdi, tx_id, config, CHR22_ENSEMBL108_TXS)
    }

    /// `annotate_chr22_window` with the predictor configuration `config` and the transcripts
    /// in `tx_path`.
    fn annotate_chr22_window_with_config(
        spdi: &str,
        tx_id: &str,
        config: Config,
        tx_path: &str,
    ) -> Result<AnnField, anyhow::Error> {
        let spdi = spdi.split(':').collect::<Vec<_>>();

        let dir = tempfile::tempdir()?;
        let reference = write_chr22_from_windows(
            "tests/data/annotate/seqvars/exon-edge.chr22-windows.fa",
            dir.path(),
        )?;
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(tx_db, Some(reference), false));
        let predictor = ConsequencePredictor::new(provider, config);

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
    /// shifts them 3'. The terms and `p.` follow the 3'-shifted HGVS description, also where the
    /// same change at the other end of its repeat would keep a splice site, see
    /// `splice_site_intact_at_equivalent_position`. With VEP terms, as VEP 116.3 with
    /// `--shift_3prime 1`.
    #[rstest::rstest]
    // deleting the acceptor G equals deleting the first exon G
    #[case("22:37713208:AG:A", "ENST00000644935", false, "c.257del", "p.(Gly86AlafsTer33)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // deleting the donor G (`GCCG|GTGAGT`), which equals deleting the last exon G
    #[case("22:23772977:CG:C", "ENST00000215743", false, "c.108+1del", "p.?", vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant])]
    #[case("22:23772977:CG:C", "ENST00000215743", true, "c.108+1del", "p.?", vec![Consequence::SpliceDonorVariant])]
    // minus strand: deleting donor +1 to +3, which equals deleting the last three exon bases
    #[case("22:26481718:TCAC:T", "ENST00000336873", false, "c.41+1_41+3del", "p.?", vec![Consequence::SpliceDonorVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    // minus strand: deleting acceptor -3 and -2, which equals an intronic deletion that keeps
    // the AG
    #[case("22:50246093:CTG:C", "ENST00000216271", false, "c.1651-3_1651-2del", "p.?", vec![Consequence::SpliceAcceptorVariant, Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant])]
    // deletion of donor +3 to +9 that keeps the GT
    #[case("22:25763388:AGGTAAGT:A", "ENST00000335473", false, "c.198+3_198+9del", "p.?", vec![Consequence::SpliceDonorFifthBaseVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    // deleting one C of `CCCAG|GT`: the 3'-shifted C is one of the last three exon bases
    #[case("22:39101499:TC:T", "ENST00000442487", false, "c.416del", "p.(Pro139GlnfsTer25)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // the deletion removes donor +1 and +2 at every position
    #[case("22:28773775:ATGAG:A", "ENST00000249064", false, "c.238_239+2del", "p.?", vec![Consequence::SpliceDonorVariant, Consequence::ExonicSpliceRegionVariant])]
    // insertion between the last exon base and donor +1: the exon gains a base
    #[case("22:20996112:T:TA", "ENST00000646124", false, "c.2219_2219+1insA", "p.(Tyr741LeufsTer41)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    #[case("22:20996112:T:TA", "ENST00000646124", true, "c.2219_2219+1insA", "p.(Tyr741LeufsTer41)", vec![Consequence::FrameshiftVariant, Consequence::SpliceRegionVariant])]
    // minus strand: insertion between donor +1 and the last exon base
    #[case("22:32150991:C:CTT", "ENST00000382097", false, "c.493_493+1insAA", "p.(His165GlnfsTer19)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // insertion one base past the acceptor, between the first two exon bases
    #[case("22:50522343:G:GCTTTCTC", "ENST00000299821", false, "c.1234_1235insCTTTCTC", "p.(Val412AlafsTer5)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // insertion between acceptor -1 and the first CDS base: the 5' UTR gains a base
    #[case("22:28987075:G:GC", "ENST00000402174", false, "c.1-1_1insC", "p.?", vec![Consequence::ExonicSpliceRegionVariant, Consequence::FivePrimeUtrExonVariant])]
    // insertion between the last stop codon base and donor +1: the 3' UTR gains a base
    #[case("22:31879752:G:GA", "ENST00000646998", false, "c.3858_3858+1insA", "p.?", vec![Consequence::ExonicSpliceRegionVariant, Consequence::ThreePrimeUtrExonVariant])]
    // minus strand: dup of the last exon base, in front of donor +1 or inside the exon
    #[case("22:28710005:C:CA", "ENST00000348295", false, "c.846dup", "p.(Pro283SerfsTer7)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // dup inside the stop codon that equals one behind it (`TAA` to `TAAAA`)
    #[case("22:21469819:T:TAA", "ENST00000432134", false, "c.932_933dup", "p.?", vec![Consequence::ThreePrimeUtrExonVariant])]
    // frameshift inside the stop codon: every position changes the stop codon
    #[case("22:42571197:T:TG", "ENST00000340239", false, "c.611dup", "p.(Ter204TrpextTer73)", vec![Consequence::StopLost, Consequence::FeatureElongation])]
    // dup of donor +1 to +7: the copy lands behind +7 and leaves the donor intact. It starts
    // with the GT of the donor, so hgvs-rs keeps it in the intron, too.
    #[case("22:31604404:G:GGTAGGAA", "ENST00000400288", false, "c.1977+1_1977+7dup", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    #[case("22:31604404:G:GGTAGGAA", "ENST00000400288", true, "c.1977+1_1977+7dup", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::IntronVariant])]
    // dup of donor +5 and +6: the copy lands behind +6 and leaves donor +3 to +6 intact
    #[case("22:31616877:G:GGT", "ENST00000400288", false, "c.3433+5_3433+6dup", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    // dup of c.16 to c.20: the copy lands behind c.20, in the last three exon bases
    #[case("22:42554033:G:GCCGCA", "ENST00000340239", false, "c.16_20dup", "p.(Gly8HisfsTer4)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    fn annotate_exon_edge_indel_csqs(
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

    /// With VEP terms, a change of donor +5 is also an `intron_variant`, as VEP counts the intron
    /// bases from +3 on as intronic (`_intron_effects`).
    #[test]
    fn vep_fifth_base_is_intronic() -> Result<(), anyhow::Error> {
        let ann = annotate_chr22_window("22:25763388:AGGTAAGT:A", "ENST00000335473", true)?;
        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (
                Some("c.198+3_198+9del"),
                &vec![
                    Consequence::SpliceDonorFifthBaseVariant,
                    Consequence::IntronVariant
                ]
            )
        );
        Ok(())
    }

    /// `splice_polypyrimidine_tract_variant` covers acceptor -17 to -3. mehari counts an
    /// insertion as in the tract if its bases go there, and reports the tract also for a
    /// variant that reaches into the exon. With VEP terms, as VEP 116.3: an insertion is in the
    /// tract if the base on either side is (`_intron_effects`), and no variant that overlaps an
    /// exon is in it (`include => { exon => 0 }`), also next to `splice_acceptor_variant`.
    #[rstest::rstest]
    // plus strand: insertion between -18 and -17
    #[case("22:28773707:T:TC", "ENST00000249064", "c.186-18_186-17insC", [false, true])]
    // plus strand: deletion of -3 and -2
    #[case("22:28773721:ACA:A", "ENST00000249064", "c.186-3_186-2del", [true, true])]
    // plus strand: deletion of -3 to the first exon base
    #[case("22:28773721:ACAGT:A", "ENST00000249064", "c.186-3_186del", [true, false])]
    // minus strand: insertion between -18 and -17
    #[case("22:50246109:G:GA", "ENST00000216271", "c.1651-18_1651-17insT", [false, true])]
    // minus strand: deletion of -3 and -2
    #[case("22:50246093:CTG:C", "ENST00000216271", "c.1651-3_1651-2del", [true, true])]
    // minus strand: deletion of -3 to the first exon base
    #[case("22:50246091:ATCTG:A", "ENST00000216271", "c.1651-3_1651del", [true, false])]
    fn polypyrimidine_tract(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_hgvs_c: &str,
        #[case] expected_in_tract: [bool; 2],
    ) -> Result<(), anyhow::Error> {
        for (vep_consequence_terms, expected_in_tract) in
            [false, true].into_iter().zip(expected_in_tract)
        {
            let ann = annotate_chr22_window(spdi, tx_id, vep_consequence_terms)?;
            assert_eq!(
                (
                    ann.hgvs_c.as_deref(),
                    ann.consequences
                        .contains(&Consequence::SplicePolypyrimidineTractVariant)
                ),
                (Some(expected_hgvs_c), expected_in_tract),
                "spdi = {spdi}, VEP terms = {vep_consequence_terms}, {:?}",
                ann.consequences
            );
        }
        Ok(())
    }

    /// The 3' shift can put an indel onto an essential splice site that the same change at the
    /// other end of its repeat leaves intact. The terms follow the HGVS description, and a
    /// message says so.
    #[rstest::rstest]
    // c.108+1del equals c.108del
    #[case("22:23772977:CG:C", "ENST00000215743", true)]
    // minus strand: c.41+1_41+3del equals deleting the last three exon bases
    #[case("22:26481718:TCAC:T", "ENST00000336873", true)]
    // minus strand: c.1651-3_1651-2del equals an intronic deletion that keeps the AG
    #[case("22:50246093:CTG:C", "ENST00000216271", true)]
    // c.239+1dup joins the exon as c.239_239+1insG, which keeps the donor, so it has no splice
    // site term
    #[case("22:28773779:G:GG", "ENST00000249064", false)]
    // c.239_239+1dup breaks the donor, but equals c.238_239insAG, which keeps it
    #[case("22:28773779:G:GAG", "ENST00000249064", true)]
    // c.238_239+2del removes donor +1 and +2 at every position
    #[case("22:28773775:ATGAG:A", "ENST00000249064", false)]
    // c.416del has no splice site term
    #[case("22:39101499:TC:T", "ENST00000442487", false)]
    fn splice_site_intact_at_equivalent_position(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let ann = annotate_chr22_window(spdi, tx_id, false)?;
        let has_message = ann.messages.as_ref().is_some_and(|messages| {
            messages.contains(&Message::InfoSpliceSiteIntactAtEquivalentPosition)
        });
        assert_eq!(has_message, expected, "{ann:?}");

        Ok(())
    }

    /// The variants of `annotate_exon_edge_indel_csqs` at the same exon edge of other transcripts
    /// of the gene (GENCODE 42, i.e., Ensembl 108): in the CDS, in a UTR and on a non-coding
    /// transcript. The c. and the n. variants get the same terms apart from the region.
    #[rstest::rstest]
    // deleting the donor G (`GCCG|GTGAGT`), which equals deleting the last exon G
    #[case("22:23772977:CG:C", "ENST00000215743", false, "c.108+1del", vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant])]
    #[case("22:23772977:CG:C", "ENST00000428253", false, "n.133+1del", vec![Consequence::SpliceDonorVariant, Consequence::NonCodingTranscriptIntronVariant])]
    #[case("22:23772977:CG:C", "ENST00000428253", true, "n.133+1del", vec![Consequence::SpliceDonorVariant, Consequence::NonCodingTranscriptVariant])]
    // minus strand: deleting donor +1 to +3, which equals deleting the last three exon bases
    #[case("22:26481718:TCAC:T", "ENST00000336873", false, "c.41+1_41+3del", vec![Consequence::SpliceDonorVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    #[case("22:26481718:TCAC:T", "ENST00000483631", false, "c.-664+1_-664+3del", vec![Consequence::SpliceDonorVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::FivePrimeUtrIntronVariant, Consequence::CodingTranscriptIntronVariant])]
    #[case("22:26481718:TCAC:T", "ENST00000466781", false, "n.612+1_612+3del", vec![Consequence::SpliceDonorVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::NonCodingTranscriptIntronVariant])]
    // minus strand: deleting acceptor -3 and -2, which equals an intronic deletion that keeps
    // the AG
    #[case("22:50246093:CTG:C", "ENST00000216271", false, "c.1651-3_1651-2del", vec![Consequence::SpliceAcceptorVariant, Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant])]
    #[case("22:50246093:CTG:C", "ENST00000415993", false, "c.*1172-3_*1172-2del", vec![Consequence::SpliceAcceptorVariant, Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::ThreePrimeUtrIntronVariant, Consequence::CodingTranscriptIntronVariant])]
    #[case("22:50246093:CTG:C", "ENST00000477814", false, "n.1286-3_1286-2del", vec![Consequence::SpliceAcceptorVariant, Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::NonCodingTranscriptIntronVariant])]
    fn annotate_exon_edge_indel_csqs_on_other_transcripts(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] vep_consequence_terms: bool,
        #[case] expected_hgvs_c: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let config = ConfigBuilder::default()
            .vep_consequence_terms(vep_consequence_terms)
            .vep_compatibility_mode(vep_consequence_terms)
            .build()?;
        let ann = annotate_chr22_window_with_config(spdi, tx_id, config, CHR22_GENCODE42_TXS)?;
        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (Some(expected_hgvs_c), &expected_csqs),
            "spdi = {spdi}"
        );

        Ok(())
    }

    /// With `discard_utr_splice_variants`, a variant in a UTR loses its splice terms. The
    /// region of hgvs-rs decides, as for the UTR terms, not the c. description.
    #[rstest::rstest]
    // insertion between acceptor -1 and the first CDS base: the 5' UTR gains a base
    #[case("22:28987075:G:GC", "ENST00000402174", CHR22_ENSEMBL108_TXS, "c.1-1_1insC", vec![Consequence::FivePrimeUtrExonVariant])]
    // insertion between the last stop codon base and donor +1: the 3' UTR gains a base
    #[case("22:31879752:G:GA", "ENST00000646998", CHR22_ENSEMBL108_TXS, "c.3858_3858+1insA", vec![Consequence::ThreePrimeUtrExonVariant])]
    // dup of acceptor -1, the copy in front of the first CDS base: the 5' UTR gains a base
    #[case("22:28987075:G:GG", "ENST00000402174", CHR22_ENSEMBL108_TXS, "c.1-1dup", vec![Consequence::FivePrimeUtrExonVariant])]
    // dup of donor +1 behind the stop codon: the copy joins the 3' UTR, as the insertion does
    #[case("22:31879752:G:GG", "ENST00000646998", CHR22_ENSEMBL108_TXS, "c.3858+1dup", vec![Consequence::ThreePrimeUtrExonVariant])]
    // deleting the donor G behind the CDS: not in a UTR, so the splice term stays
    #[case("22:23772977:CG:C", "ENST00000215743", CHR22_ENSEMBL108_TXS, "c.108+1del", vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant])]
    // minus strand: deleting donor +1 to +3 of a 5' UTR intron
    #[case("22:26481718:TCAC:T", "ENST00000483631", CHR22_GENCODE42_TXS, "c.-664+1_-664+3del", vec![Consequence::FivePrimeUtrIntronVariant, Consequence::CodingTranscriptIntronVariant])]
    // minus strand: deletion in a 3' UTR intron that keeps the AG
    #[case("22:50246093:CTG:C", "ENST00000415993", CHR22_GENCODE42_TXS, "c.*1172-3_*1172-2del", vec![Consequence::ThreePrimeUtrIntronVariant, Consequence::CodingTranscriptIntronVariant])]
    fn annotate_discard_utr_splice_csqs(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] tx_path: &str,
        #[case] expected_hgvs_c: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let config = ConfigBuilder::default()
            .discard_utr_splice_variants(true)
            .build()?;
        let ann = annotate_chr22_window_with_config(spdi, tx_id, config, tx_path)?;
        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (Some(expected_hgvs_c), &expected_csqs),
            "spdi = {spdi}"
        );

        Ok(())
    }

    /// Dups on GENCODE 42 chr22 transcripts. The terms of a dup come from the insertion of its
    /// copy behind its bases. With VEP terms, the expected terms are those of VEP 116.3 with
    /// `--shift_3prime 1`.
    #[rstest::rstest]
    // ENST00000413051 has no 5' UTR, and the copy lands behind c.2, in the start codon. The copy
    // can also go in front of c.1, so hgvs-rs places the dup in the 5' UTR and gives `p.?`, and
    // the start codon remains. With VEP terms, the position of the copy counts, not its 5' UTR
    // region.
    #[case("22:37746290:C:CGCGGCGGCGGGCGGCCGAGAGGG", "ENST00000413051", false, "c.-21_2dup", "p.?", vec![Consequence::StartRetainedVariant, Consequence::FivePrimeUtrExonVariant])]
    #[case("22:37746290:C:CGCGGCGGCGGGCGGCCGAGAGGG", "ENST00000413051", true, "c.-21_2dup", "p.?", vec![Consequence::FrameshiftVariant])]
    // the copy lands in the last three exon bases
    #[case("22:50256598:C:CCTTCAGCTCACAGT", "ENST00000215659", false, "c.491_504dup", "p.(Lys168_Ile169insThrValSerTer)", vec![Consequence::StopGained, Consequence::FrameshiftVariant, Consequence::ExonicSpliceRegionVariant])]
    // VEP's changed codons at the 3'-shifted position hold the new stop codon
    #[case("22:50256598:C:CCTTCAGCTCACAGT", "ENST00000215659", true, "c.491_504dup", "p.(Lys168_Ile169insThrValSerTer)", vec![Consequence::StopGained, Consequence::FrameshiftVariant, Consequence::SpliceRegionVariant])]
    fn annotate_dup_csqs(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] vep_consequence_terms: bool,
        #[case] expected_hgvs_c: &str,
        #[case] expected_hgvs_p: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let config = ConfigBuilder::default()
            .vep_consequence_terms(vep_consequence_terms)
            .vep_compatibility_mode(vep_consequence_terms)
            .build()?;
        let ann = annotate_chr22_window_with_config(spdi, tx_id, config, CHR22_GENCODE42_TXS)?;
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

    /// `vep_consequence_terms` alone renames the terms but keeps mehari's rules; only
    /// `vep_compatibility_mode` follows VEP. The dup can go in front of c.1, so mehari puts it in
    /// the 5' UTR, while VEP counts its position in the CDS, see `annotate_dup_csqs`.
    #[test]
    fn vep_terms_without_compatibility_mode() -> Result<(), anyhow::Error> {
        let config = ConfigBuilder::default()
            .vep_consequence_terms(true)
            .build()?;
        let ann = annotate_chr22_window_with_config(
            "22:37746290:C:CGCGGCGGCGGGCGGCCGAGAGGG",
            "ENST00000413051",
            config,
            CHR22_GENCODE42_TXS,
        )?;
        assert_eq!(
            ann.consequences,
            [
                Consequence::StartRetainedVariant,
                Consequence::FivePrimeUtrVariant
            ]
        );

        Ok(())
    }

    /// Insertions at exon edges on GENCODE 42 chr22 transcripts, with VEP terms. VEP puts such an
    /// insertion into the exon and takes its stop terms and in-frame insertion terms from the
    /// changed codons there. The expected terms are those of VEP 116.3 with `--shift_3prime 1`.
    #[rstest::rstest]
    // ENST00000381982 (minus strand): an exon ends with c.1081, the first base of the stop codon
    // (`t|ag`). `*/Y*` (`tag/tACTag`): the stop codon stays behind the new codon.
    #[case("22:30491572:C:CAGT", "ENST00000381982", "c.1081_1081+1insACT", vec![Consequence::InframeInsertion, Consequence::SpliceRegionVariant])]
    // `*/*LX` (`tag/tAACTag`): the bases at the old stop codon still read a stop codon. On the
    // genome, an intron lies between the input position and the stop codon end, so VEP finds no
    // overlap with the stop codon.
    #[case("22:30491572:C:CAGTT", "ENST00000381982", "c.1081_1081+1insAACT", vec![Consequence::FrameshiftVariant, Consequence::SpliceRegionVariant])]
    // `*/YX` (`tag/tATag`)
    #[case("22:30491572:C:CAT", "ENST00000381982", "c.1081_1081+1insAT", vec![Consequence::FrameshiftVariant, Consequence::StopLost, Consequence::SpliceRegionVariant])]
    // `*/**` (`tag/tAATag`)
    #[case("22:30491572:C:CATT", "ENST00000381982", "c.1081_1081+1insAAT", vec![Consequence::InframeInsertion, Consequence::SpliceRegionVariant, Consequence::StopRetainedVariant])]
    // `-/K*` (`-/AAATAA`): in front of the exon that starts with c.520
    #[case("22:30494210:A:ATTATTT", "ENST00000381982", "c.520-1_520insAAATAA", vec![Consequence::StopGained, Consequence::InframeInsertion, Consequence::SpliceRegionVariant])]
    // `F/*F` (`ttt/tAGTtt`): the copy lands behind c.520
    #[case("22:30494209:A:AACT", "ENST00000381982", "c.520-2_520dup", vec![Consequence::StopGained, Consequence::SpliceRegionVariant])]
    // ENST00000398941 (plus strand): an exon ends with c.428, the next one starts with c.429.
    // `S/SX` (`agc/agCGc`)
    #[case("22:28744708:A:AGC", "ENST00000398941", "c.428_428+1insCG", vec![Consequence::FrameshiftVariant, Consequence::SpliceRegionVariant])]
    // `-/X` (`-/GC`): the copy lands behind c.429
    #[case("22:28751239:A:AGC", "ENST00000398941", "c.429-1_429dup", vec![Consequence::FrameshiftVariant, Consequence::SpliceRegionVariant])]
    // `K/N*` (`aag/aaCTAg`)
    #[case("22:37967169:A:AACT", "ENST00000442738", "c.293_293+1insCTA", vec![Consequence::StopGained, Consequence::ProteinAlteringVariant, Consequence::SpliceRegionVariant])]
    // minus strand: `-/P` (`-/CCG`)
    #[case("22:50775771:C:CCGG", "ENST00000395590", "c.297_297+1insCCG", vec![Consequence::InframeInsertion, Consequence::SpliceRegionVariant])]
    fn annotate_exon_edge_ins_stop_csqs(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_hgvs_c: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let config = ConfigBuilder::default()
            .vep_consequence_terms(true)
            .vep_compatibility_mode(true)
            .build()?;
        let ann = annotate_chr22_window_with_config(spdi, tx_id, config, CHR22_GENCODE42_TXS)?;
        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (Some(expected_hgvs_c), &expected_csqs),
            "spdi = {spdi}"
        );

        Ok(())
    }

    /// Insertions on GENCODE 42 chr22 transcripts, with VEP terms. VEP counts an insertion as
    /// inside a region only if the bases on both sides of it are inside. The expected terms are
    /// those of VEP 116.3 with `--shift_3prime 1`.
    #[rstest::rstest]
    // the copy lands between donor +2 and +3, so only one side is in donor +3 to +6
    #[case("22:15690709:G:GGT", "ENST00000343518", "c.632+1_632+2dup", vec![Consequence::SpliceRegionVariant, Consequence::IntronVariant])]
    // between donor +4 and +5: inside donor +3 to +6, but not at the fifth base
    #[case("22:15690713:A:AG", "ENST00000343518", "c.632+4_632+5insG", vec![Consequence::SpliceDonorRegionVariant, Consequence::IntronVariant])]
    // minus strand
    #[case("22:26481719:C:CG", "ENST00000466781", "n.612+2_612+3insC", vec![Consequence::SpliceRegionVariant, Consequence::NonCodingTranscriptVariant, Consequence::IntronVariant])]
    #[case("22:26481717:C:CA", "ENST00000466781", "n.612+4_612+5insT", vec![Consequence::SpliceDonorRegionVariant, Consequence::NonCodingTranscriptVariant, Consequence::IntronVariant])]
    // between an exon and an intron of a non-coding transcript: not in the exon
    #[case("22:20996112:T:TA", "ENST00000495142", "n.2571_2571+1insA", vec![Consequence::SpliceRegionVariant, Consequence::NonCodingTranscriptVariant])]
    #[case("22:28710005:C:CA", "ENST00000491919", "n.403dup", vec![Consequence::SpliceRegionVariant, Consequence::NonCodingTranscriptVariant])]
    #[case("22:26481721:C:CA", "ENST00000466781", "n.612_612+1insT", vec![Consequence::SpliceRegionVariant, Consequence::NonCodingTranscriptVariant])]
    #[case("22:26479355:C:CA", "ENST00000466781", "n.613-1_613insT", vec![Consequence::SpliceRegionVariant, Consequence::NonCodingTranscriptVariant])]
    // between the last two exon bases
    #[case("22:26481722:C:CG", "ENST00000466781", "n.611_612insC", vec![Consequence::SpliceRegionVariant, Consequence::NonCodingTranscriptExonVariant])]
    fn annotate_vep_ins_region_csqs(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_hgvs_c: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let config = ConfigBuilder::default()
            .vep_consequence_terms(true)
            .vep_compatibility_mode(true)
            .build()?;
        let ann = annotate_chr22_window_with_config(spdi, tx_id, config, CHR22_GENCODE42_TXS)?;
        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (Some(expected_hgvs_c), &expected_csqs),
            "spdi = {spdi}"
        );

        Ok(())
    }

    /// Indels at the end of the CDS on transcripts of `CHR22_CDS_END_TXS`, with VEP terms. The
    /// expected terms are those of VEP 116.3 with `--shift_3prime 1`.
    #[rstest::rstest]
    // The deletion shifts from `c.1637_*1del` to `c.1638_*1+1del` and keeps the stop codon
    // `TGA`. VEP builds no codons for it, as its end lies in the intron.
    #[case("22:15719775:TGAG:T", "ENST00000343518", Some("c.1638_*1+1del"), vec![Consequence::SpliceDonorVariant, Consequence::StopRetainedVariant, Consequence::ThreePrimeUtrVariant])]
    // ENST00000381982 (minus strand): an exon ends with c.1081, the first base of the stop codon
    // (`t|ag`). VEP takes the stop codon from the genome, so c.1081 lies outside it.
    #[case("22:30491570:TACAG:T", "ENST00000381982", Some("c.1080_1081+2del"), vec![Consequence::SpliceDonorVariant, Consequence::CodingSequenceVariant])]
    // ENST00000613577 and ENST00000540833 end with the stop codon `TAA`. The insertion shifts
    // behind the last base, so VEP puts it outside the transcript. The dup that describes it
    // stays inside.
    #[case("22:18530457:T:TA", "ENST00000613577", Some("c.1071dup"), vec![Consequence::DownstreamGeneVariant])]
    #[case("22:41767606:T:TA", "ENST00000540833", Some("c.1521dup"), vec![Consequence::DownstreamGeneVariant])]
    // an insertion behind the last base that does not shift
    #[case("22:18530459:A:AC", "ENST00000613577", None, vec![Consequence::DownstreamGeneVariant])]
    // also on the minus strand
    #[case("22:30489395:A:AC", "ENST00000381982", None, vec![Consequence::DownstreamGeneVariant])]
    fn annotate_vep_cds_end_csqs(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_hgvs_c: Option<&str>,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let config = ConfigBuilder::default()
            .vep_consequence_terms(true)
            .vep_compatibility_mode(true)
            .build()?;
        let ann = annotate_chr22_window_with_config(spdi, tx_id, config, CHR22_CDS_END_TXS)?;
        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (expected_hgvs_c, &expected_csqs),
            "spdi = {spdi}"
        );

        Ok(())
    }

    /// An insertion at an exon edge, or the copy of a dup that goes there, takes the location of
    /// the side its bases join, see `annotate_exon_boundary_dup_and_ins`: the exon number and a
    /// CDS position, or the intron number and the distance to the exon.
    #[rstest::rstest]
    #[case("22:32150991:C:CTT", "ENST00000382097", "c.493_493+1insAA", Rank { ord: 6, total: 9 }, 0, Some(Pos { ord: 493, total: Some(756) }))]
    #[case("22:28773779:G:GG", "ENST00000249064", "c.239+1dup", Rank { ord: 2, total: 5 }, 0, Some(Pos { ord: 239, total: Some(840) }))]
    #[case("22:28773724:G:GG", "ENST00000249064", "c.186-1dup", Rank { ord: 2, total: 5 }, 0, Some(Pos { ord: 186, total: Some(840) }))]
    #[case("22:28773724:G:GAAG", "ENST00000249064", "c.186-1_186insAAG", Rank { ord: 1, total: 4 }, -1, None)]
    #[case("22:28773724:G:GCAG", "ENST00000249064", "c.186-3_186-1dup", Rank { ord: 1, total: 4 }, -1, None)]
    fn annotate_ins_at_exon_edge_location(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_hgvs_c: &str,
        #[case] expected_rank: Rank,
        #[case] expected_distance: i32,
        #[case] expected_cds_pos: Option<Pos>,
    ) -> Result<(), anyhow::Error> {
        let ann = annotate_chr22_window(spdi, tx_id, false)?;
        assert_eq!(
            (ann.hgvs_c.as_deref(), ann.rank, ann.distance, ann.cds_pos),
            (
                Some(expected_hgvs_c),
                Some(expected_rank),
                Some(expected_distance),
                expected_cds_pos
            ),
            "spdi = {spdi}"
        );

        Ok(())
    }

    /// The alternative sequences follow the 3'-shifted HGVS description, as the terms do. Deleting
    /// the donor G of `GCCG|GTGAGT` is `c.108+1del`, which changes splicing, so the alternative
    /// transcript is unknown.
    #[test]
    fn alt_sequences_follow_hgvs() -> Result<(), anyhow::Error> {
        let config = ConfigBuilder::default()
            .plugins(vec![
                Arc::new(SequencePlugin::TxRef),
                Arc::new(SequencePlugin::TxAlt),
                Arc::new(SequencePlugin::AaRef),
                Arc::new(SequencePlugin::AaAlt),
            ])
            .build()?;
        let ann = annotate_chr22_window_with_config(
            "22:23772977:CG:C",
            "ENST00000215743",
            config,
            CHR22_ENSEMBL108_TXS,
        )?;
        assert_eq!(ann.hgvs_c.as_deref(), Some("c.108+1del"));
        assert_eq!(ann.hgvs_p.as_deref(), Some("p.?"));

        assert!(custom_field(&ann, ANN_TX_SEQ_REF).is_some());
        assert!(custom_field(&ann, ANN_AA_SEQ_REF).is_some());
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), None);
        assert_eq!(custom_field(&ann, ANN_AA_SEQ_ALT), None);

        Ok(())
    }

    /// The alternative sequences and NMD follow the reading of hgvs-rs at an exon edge, see
    /// `annotate_exon_boundary_dup_and_ins`. The copy of `c.239+1dup` joins the exon, so the
    /// transcript gains a G behind c.239 (n.394), and the frameshift triggers NMD. The copy of
    /// `c.186-3_186-1dup` stays in the intron, so the alternative transcript is unknown.
    #[test]
    fn alt_sequences_at_exon_edge() -> Result<(), anyhow::Error> {
        let annotate = |spdi| -> Result<AnnField, anyhow::Error> {
            let config = ConfigBuilder::default()
                .plugins(vec![
                    Arc::new(SequencePlugin::TxRef),
                    Arc::new(SequencePlugin::TxAlt),
                    Arc::new(NmdPlugin {
                        vep_compatible: false,
                    }),
                ])
                .build()?;
            annotate_chr22_window_with_config(spdi, "ENST00000249064", config, CHR22_ENSEMBL108_TXS)
        };

        let ann = annotate("22:28773779:G:GG")?;
        assert_eq!(ann.hgvs_c.as_deref(), Some("c.239+1dup"));
        let mut tx_alt = custom_field(&ann, ANN_TX_SEQ_REF).unwrap().to_string();
        tx_alt.insert(394, 'G');
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), Some(tx_alt.as_str()));
        assert_eq!(custom_field(&ann, ANN_NMD), Some("triggering"));

        let ann = annotate("22:28773724:G:GCAG")?;
        assert_eq!(ann.hgvs_c.as_deref(), Some("c.186-3_186-1dup"));
        assert!(custom_field(&ann, ANN_TX_SEQ_REF).is_some());
        assert_eq!(custom_field(&ann, ANN_TX_SEQ_ALT), None);

        Ok(())
    }

    /// Dups and insertions at the edges of a coding exon. A dup of the first or last exon base,
    /// or of intron bases next to the exon, equals an insertion of its copy at the exon edge. As
    /// in hgvs-rs, bases at the edge stay in the intron if they keep the splice site dinucleotide
    /// right next to the exon, with `p.?`. Else they join the exon, which keeps the splice site
    /// behind them, and the frame shifts. An insertion between the two donor site bases breaks
    /// the site.
    ///
    /// Exon 2 of ENST00000249064 (plus strand) spans c.186 to c.239 (`CAG|TGT`, `GATGA|GTAAGT`).
    /// On the minus strand, exon 7 of ENST00000348295 ends at c.846 (`CAT|GTAAGT`), and exon 18
    /// of ENST00000216271 starts at c.1651 (`CAG|ATG`).
    #[rstest::rstest]
    // codon 62 CG|T stays, codon 63 GTT (Val) becomes TGT (Cys)
    #[case("22:28773725:T:TT", "ENST00000249064", "c.186dup", "p.(Val63CysfsTer6)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 62 becomes CGG (Arg, as CGT), codon 63 TGT (Cys) as for c.186dup
    #[case("22:28773724:G:GG", "ENST00000249064", "c.186-1dup", "p.(Val63CysfsTer6)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 62 becomes CGC (Arg, as CGT), codon 63 TGT (Cys) as for c.186dup
    #[case("22:28773724:G:GC", "ENST00000249064", "c.186-1_186insC", "p.(Val63CysfsTer6)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // the inserted bases end with the AG of the acceptor, so they stay in the intron
    #[case("22:28773724:G:GAAG", "ENST00000249064", "c.186-1_186insAAG", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    // the copy ends with the AG of the acceptor
    #[case("22:28773724:G:GCAG", "ENST00000249064", "c.186-3_186-1dup", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    // codon 80 GA|T (Asp, its T in exon 3) becomes GAA (Glu)
    #[case("22:28773778:A:AA", "ENST00000249064", "c.239dup", "p.(Asp80GlufsTer11)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 80 becomes GAG (Glu), the frame as for c.239dup
    #[case("22:28773779:G:GG", "ENST00000249064", "c.239+1dup", "p.(Asp80GlufsTer11)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // the copy starts with the GT of the donor, so it stays in the intron, behind donor +2
    #[case("22:28773778:A:AGT", "ENST00000249064", "c.239+1_239+2dup", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    // codon 80 becomes GAC (Asp, as GAT), codon 81 TGT (Cys) becomes TTG (Leu)
    #[case("22:28773778:A:AC", "ENST00000249064", "c.239_239+1insC", "p.(Cys81LeufsTer10)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    #[case("22:28773779:G:GC", "ENST00000249064", "c.239+1_239+2insC", "p.?", vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant])]
    // minus strand: codon 282 CAT stays, codon 283 CCT (Pro) becomes TCC (Ser)
    #[case("22:28710005:C:CA", "ENST00000348295", "c.846dup", "p.(Pro283SerfsTer7)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 283 becomes GCC (Ala)
    #[case("22:28710005:C:CC", "ENST00000348295", "c.846+1dup", "p.(Pro283AlafsTer7)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // the copy starts with the GT of the donor
    #[case("22:28710005:C:CAC", "ENST00000348295", "c.846+1_846+2dup", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant])]
    // codon 283 becomes CCC (Pro, as CCT), codon 284 TGC (Cys) becomes TTG (Leu)
    #[case("22:28710005:C:CG", "ENST00000348295", "c.846_846+1insC", "p.(Cys284LeufsTer6)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    #[case("22:28710004:A:AG", "ENST00000348295", "c.846+1_846+2insC", "p.?", vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant])]
    // minus strand: codon 551 ATG (Met) becomes AAT (Asn)
    #[case("22:50246092:T:TT", "ENST00000216271", "c.1651dup", "p.(Met551AsnfsTer23)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 551 becomes GAT (Asp)
    #[case("22:50246093:C:CC", "ENST00000216271", "c.1651-1dup", "p.(Met551AspfsTer23)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // codon 551 becomes CAT (His)
    #[case("22:50246092:T:TG", "ENST00000216271", "c.1651-1_1651insC", "p.(Met551HisfsTer23)", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation, Consequence::ExonicSpliceRegionVariant])]
    // the inserted bases end with the AG of the acceptor
    #[case("22:50246092:T:TCTC", "ENST00000216271", "c.1651-1_1651insGAG", "p.?", vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant])]
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

    /// A variant with both ends in introns covers the exon bases between them. It is a UTR exon
    /// variant if one of these lies in a UTR, e.g. c.*1 to c.*5 for `c.112-59_*6-51del`. The
    /// first of them lies behind the intron of its start, the last in front of that of its end.
    #[rstest::rstest]
    #[case(
        "c.112-59_*6-51del",
        VariantLocation::Intron,
        true,
        Consequence::ThreePrimeUtrExonVariant | Consequence::ThreePrimeUtrIntronVariant
    )]
    #[case(
        "c.*160+197_*437-8del",
        VariantLocation::ThreePrimeUtr,
        true,
        Consequence::ThreePrimeUtrExonVariant | Consequence::ThreePrimeUtrIntronVariant
    )]
    #[case(
        "c.-89-3923_201+2793del",
        VariantLocation::Intron,
        true,
        Consequence::FivePrimeUtrExonVariant | Consequence::FivePrimeUtrIntronVariant
    )]
    #[case(
        "c.-19-2025_120+884delinsTTA",
        VariantLocation::Intron,
        true,
        Consequence::FivePrimeUtrExonVariant | Consequence::FivePrimeUtrIntronVariant
    )]
    #[case(
        "c.-90+5_-50-10del",
        VariantLocation::FivePrimeUtr,
        true,
        Consequence::FivePrimeUtrExonVariant | Consequence::FivePrimeUtrIntronVariant
    )]
    // The exons between start with c.1 and end with c.3857.
    #[case(
        "c.-1+5_100+5del",
        VariantLocation::Intron,
        true,
        Consequence::FivePrimeUtrIntronVariant.into()
    )]
    #[case(
        "c.110+5_*1-10del",
        VariantLocation::Intron,
        true,
        Consequence::ThreePrimeUtrIntronVariant.into()
    )]
    // within one intron
    #[case(
        "c.-10+5_-10+20del",
        VariantLocation::FivePrimeUtr,
        false,
        Consequence::FivePrimeUtrIntronVariant.into()
    )]
    fn analyze_cds_variant_utr_exon_between_introns(
        #[case] var_c: &str,
        #[case] region: VariantLocation,
        #[case] is_exonic: bool,
        #[case] expected: Consequences,
    ) -> Result<(), anyhow::Error> {
        let var_c = HgvsVariant::from_str(&format!("NM_000000.1:{var_c}"))?;
        let csqs = ConsequencePredictor::analyze_cds_variant(
            &var_c,
            Some(region),
            is_exonic,
            false,
            false,
            Some(3858),
        );
        let utr_terms = Consequence::FivePrimeUtrExonVariant
            | Consequence::FivePrimeUtrIntronVariant
            | Consequence::ThreePrimeUtrExonVariant
            | Consequence::ThreePrimeUtrIntronVariant;
        assert_eq!(csqs & utr_terms, expected, "{csqs:?}");
        Ok(())
    }

    /// An insertion or dup that can go before c.1 or after the last CDS base leaves the CDS
    /// intact, so it lies in a UTR. `c.1dup` inserts an A in front of the ATG, and `c.1_3dup` an
    /// ATG. The stop codon of NM_007294.4 ends at c.5592. The copy of `c.1dup` goes between c.1
    /// and c.2, so it changes the start codon, which remains: `start_retained_variant`.
    #[rstest::rstest]
    #[case(
        "3:193311166:G:GC",
        "NM_130837.3",
        "c.-1_1insC",
        vec![Consequence::FivePrimeUtrExonVariant]
    )] // OPA1, forward
    #[case(
        "17:41258543:T:TA",
        "NM_007297.4",
        "c.-1_1insT",
        vec![Consequence::FivePrimeUtrExonVariant]
    )] // BRCA1, reverse
    #[case(
        "3:193311167:A:AA",
        "NM_130837.3",
        "c.1dup",
        vec![Consequence::StartRetainedVariant, Consequence::FivePrimeUtrExonVariant]
    )]
    #[case(
        "17:41276112:A:AT",
        "NM_007294.4",
        "c.1dup",
        vec![Consequence::StartRetainedVariant, Consequence::FivePrimeUtrExonVariant]
    )]
    #[case(
        "3:193311169:G:GATG",
        "NM_130837.3",
        "c.1_3dup",
        vec![Consequence::FivePrimeUtrExonVariant]
    )]
    #[case(
        "17:41276110:C:CCAT",
        "NM_007294.4",
        "c.1_3dup",
        vec![Consequence::FivePrimeUtrExonVariant]
    )]
    #[case(
        "17:41197695:T:TT",
        "NM_007294.4",
        "c.5592dup",
        vec![Consequence::ThreePrimeUtrExonVariant]
    )]
    fn annotate_ins_keeping_cds(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_hgvs_c: &str,
        #[case] expected_csqs: Vec<Consequence>,
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
        assert_eq!(ann.consequences, expected_csqs);

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
                .vep_compatibility_mode(true)
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
    /// then gives `p.?`, and the default terms are UTR terms. A deletion of a start codon base
    /// keeps the start codon then: `start_retained_variant`. In NM_007294.4, c.-2, c.-1 and c.1
    /// are A bases. In NM_130837.3, c.3048 and c.*1 are A bases.
    ///
    /// VEP 116.3 with `--shift_3prime 1` keeps such a deletion at the start codon in the CDS and
    /// reports `start_retained_variant`. It shifts the one behind the stop codon into the 3' UTR.
    #[rstest::rstest]
    #[case("17:41276112:AT:A", "NM_007294.4", "c.1del", false, vec![Consequence::StartRetainedVariant, Consequence::FivePrimeUtrExonVariant])]
    #[case("17:41276112:AT:A", "NM_007294.4", "c.1del", true, vec![Consequence::FrameshiftVariant, Consequence::StartRetainedVariant])]
    #[case("17:41276112:ATT:A", "NM_007294.4", "c.-1_1del", false, vec![Consequence::StartRetainedVariant, Consequence::FivePrimeUtrExonVariant])]
    #[case("17:41276112:ATT:A", "NM_007294.4", "c.-1_1del", true, vec![Consequence::StartRetainedVariant, Consequence::FivePrimeUtrVariant])]
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
                .vep_compatibility_mode(vep_consequence_terms)
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

    /// An indel from the 5' UTR into the start codon changes a start codon base, by its 3'-shifted
    /// description. If the same change fits in front of c.1, the start codon remains in front of
    /// the CDS: `start_retained_variant`, not `start_lost`. VEP compatibility mode decides this
    /// for a deletion only, as VEP, and keeps the CDS terms of a dup. The 5' UTR of
    /// `one_exon_tx` ends with `CAGT`, and c.-3 and c.1 are A bases.
    #[rstest::rstest]
    #[case::del(
        (1018, "AGTA", "A"),
        false,
        "c.-2_1del",
        vec![Consequence::StartRetainedVariant, Consequence::FivePrimeUtrExonVariant]
    )]
    #[case::del_vep(
        (1018, "AGTA", "A"),
        true,
        "c.-2_1del",
        vec![Consequence::StartRetainedVariant, Consequence::FivePrimeUtrVariant]
    )]
    #[case::dup(
        (1021, "A", "AAGTA"),
        false,
        "c.-3_1dup",
        vec![Consequence::StartRetainedVariant, Consequence::FivePrimeUtrExonVariant]
    )]
    #[case::dup_vep(
        (1021, "A", "AAGTA"),
        true,
        "c.-3_1dup",
        vec![Consequence::FrameshiftVariant]
    )]
    fn start_retained_from_5_prime_utr(
        #[case] (position, reference, alternative): (i32, &str, &str),
        #[case] vep_compatibility_mode: bool,
        #[case] expected_c: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CGCATCA");
        let mut genome = NamedTempFile::new()?;
        write!(genome, ">1\n{}{seq}{}\n", "T".repeat(1000), "T".repeat(100))?;
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(vep_compatibility_mode)
            .build()?;
        let anns = predictor_with_reference(tx, seq, config, Some(genome.path()))
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (Some(expected_c), Some("p.?"), &expected_csqs)
        );
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
        let config = ConfigBuilder::default()
            .plugins(vec![
                Arc::new(SequencePlugin::TxRef),
                Arc::new(SequencePlugin::TxAlt),
            ])
            .build()
            .unwrap();
        predictor_with_config(tx, seq, config)
    }

    /// `predictor_for` with the predictor configuration `config`.
    fn predictor_with_config(tx: Transcript, seq: String, config: Config) -> ConsequencePredictor {
        predictor_with_reference(tx, seq, config, None)
    }

    /// `predictor_with_config` with the reference genome at `reference`.
    fn predictor_with_reference(
        tx: Transcript,
        seq: String,
        config: Config,
        reference: Option<&Path>,
    ) -> ConsequencePredictor {
        use crate::pbs::txs::{GeneToTxId, SequenceDb, SourceVersion, TranscriptDb, TxSeqDatabase};

        // hgvs-rs caches the reference protein by the database version, so each database
        // needs its own version.
        let version = format!("{tx:?} {seq}");
        let tx_seq_db = TxSeqDatabase {
            tx_db: Some(TranscriptDb {
                gene_to_tx: vec![GeneToTxId {
                    tx_ids: vec![tx.id.clone()],
                    filtered: Some(false),
                    filter_reason: None,
                    gene_id: tx.gene_id.clone(),
                    hgnc_id: tx.hgnc_id.clone(),
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
        let provider = Arc::new(MehariProvider::new(tx_seq_db, reference, true));
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

    /// A predictor for one selenoprotein transcript, see `selenoprotein_tx`.
    fn selenoprotein_predictor(tagged: bool, positions: Vec<u32>) -> ConsequencePredictor {
        let (tx, seq) = selenoprotein_tx(tagged, positions);
        predictor_for(tx, seq)
    }

    /// One selenoprotein transcript on the plus strand of chr1 (GRCh38).
    ///
    /// Its one exon spans chr1:1001-1073: 20 bases of 5' UTR, the CDS
    /// `ATG GCC AAG CTG TGG GAA CCA TGA CGC GTT TAA` (`MAKLWEPURV*`) at chr1:1021-1053, and 20
    /// bases of 3' UTR. The eighth codon is the Sec codon.
    fn selenoprotein_tx(tagged: bool, positions: Vec<u32>) -> (Transcript, String) {
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
        (tx, seq)
    }

    /// With Sec positions, UGA reads as selenocysteine only at them. Without, a transcript
    /// tagged as selenoprotein reads every UGA as selenocysteine.
    #[rstest::rstest]
    #[case::new_uga_is_a_stop(false, vec![8], "1035:G:A", "p.(Trp5Ter)", Consequence::StopGained)]
    #[case::new_uga_is_a_stop_despite_tag(
        true,
        vec![8],
        "1035:G:A",
        "p.(Trp5Ter)",
        Consequence::StopGained
    )]
    #[case::sec_codon_reads_sec(false, vec![8], "1044:A:G", "p.(Sec8Trp)", Consequence::MissenseVariant)]
    #[case::codon_after_sec_codon(
        false,
        vec![8],
        "1046:G:A",
        "p.(Arg9His)",
        Consequence::MissenseVariant
    )]
    #[case::sec_codon_moves_with_deletion(
        false,
        vec![8],
        "1035:GGAA:G",
        "p.(Glu6del)",
        Consequence::ConservativeInframeDeletion
    )]
    #[case::tag_without_positions(true, vec![], "1035:G:A", "p.(Trp5Sec)", Consequence::MissenseVariant)]
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

    /// What a plugin sees of the sequences of one transcript, through the accessors of
    /// `wit/plugin.wit`.
    #[derive(Debug, serde::Serialize)]
    struct SeenSequences {
        hgvs_c: Option<String>,
        reference_transcript: Option<String>,
        reference_protein: Option<String>,
        reference_cds_start: Option<i32>,
        reference_cds_end: Option<i32>,
        alternative_transcript: Option<String>,
        alternative_protein: Option<String>,
        alternative_cds_start: Option<i32>,
        frameshift: Option<bool>,
        first_changed_codon: Option<i32>,
    }

    /// A plugin that records what it sees of the sequences, as the WebAssembly host reads it.
    #[derive(Debug, Default)]
    struct SequenceProbe(std::sync::Mutex<Vec<SeenSequences>>);

    impl TranscriptPlugin for SequenceProbe {
        fn columns(&self) -> Vec<String> {
            Vec::new()
        }

        fn annotate(
            &self,
            _variants: &[VcfVariant],
            txs: &mut Vec<crate::annotate::seqvars::consequence::plugin::TxAnnotation>,
        ) -> Result<(), SeqvarsError> {
            let mut seen = self.0.lock().unwrap();
            for tx in txs {
                let ref_data = tx.context.ref_data();
                let alt_data = tx.context.alt_data();
                seen.push(SeenSequences {
                    hgvs_c: tx.ann.hgvs_c.clone(),
                    reference_transcript: ref_data.map(|d| d.transcript_sequence.to_string()),
                    reference_protein: ref_data.map(|d| d.aa_sequence.to_string()),
                    reference_cds_start: ref_data.map(|d| d.cds_start),
                    reference_cds_end: ref_data.map(|d| d.cds_stop),
                    alternative_transcript: alt_data.map(|d| d.transcript_sequence.clone()),
                    alternative_protein: alt_data.map(|d| d.aa_sequence.clone()),
                    alternative_cds_start: alt_data.map(|d| d.cds_start),
                    frameshift: alt_data.map(|d| d.is_frameshift),
                    first_changed_codon: alt_data.and_then(|d| d.variant_start_aa),
                });
            }
            Ok(())
        }
    }

    /// Plugins see the same reference and alternative sequences, whatever `c_to_p` computes for
    /// `p.`. The transcript has 20 bases of 5' UTR, the CDS
    /// `ATG GCC AAG CTG TGG GAA CCA CGC GTT TAA` (`MAKLWEPRV*`) at chr1:1021-1050, and 20 bases of
    /// 3' UTR; `seleno` is that of `selenoprotein_tx` with its Sec codon, and `mito` one on the
    /// mitochondrion whose CDS ends with a `T` that the poly(A) tail completes to a stop codon.
    #[rstest::rstest]
    #[case::synonymous("standard", "1026:C:T")]
    #[case::missense("standard", "1024:G:A")]
    #[case::stop_gain("standard", "1027:A:T")]
    #[case::stop_loss("standard", "1048:T:C")]
    #[case::start_loss("standard", "1021:A:G")]
    #[case::inframe_insertion("standard", "1026:C:CGGC")]
    #[case::inframe_deletion("standard", "1026:CAAG:C")]
    #[case::inframe_delins("standard", "1027:AAG:TTTCCC")]
    #[case::frameshift_insertion("standard", "1026:C:CA")]
    #[case::frameshift_deletion("standard", "1026:CA:C")]
    #[case::frameshift_into_utr3("standard", "1045:GT:G")]
    #[case::cds_end("standard", "1050:A:G")]
    #[case::n_allele("standard", "1027:AAG:NNN")]
    #[case::selenoprotein_missense("seleno", "1024:G:A")]
    #[case::selenoprotein_deletion("seleno", "1035:GGAA:G")]
    #[case::mitochondrial_poly_a_stop("mito", "1024:G:A")]
    fn plugins_see_sequences(#[case] tx: &str, #[case] var: &str) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{tx}-{}", var.replace(':', "-"));
        let utr = "CAGTCAGTCAGTCAGTCAGT";
        let (tx, seq, chromosome) = match tx {
            "seleno" => {
                let (tx, seq) = selenoprotein_tx(false, vec![8]);
                (tx, seq, "1")
            }
            "mito" => {
                let (mut tx, seq) = one_exon_tx(utr, "ATGGCCAAGCTGTGGT", "");
                tx.genome_alignments[0].contig = "NC_012920.1".into();
                (tx, seq, "MT")
            }
            _ => {
                let (tx, seq) = one_exon_tx(
                    utr,
                    "ATGGCCAAGCTGTGGGAACCACGCGTTTAA",
                    "CGCATCAGTCAGTCAGTCAG",
                );
                (tx, seq, "1")
            }
        };
        let probe = Arc::new(SequenceProbe::default());
        let config = ConfigBuilder::default()
            .plugins(vec![probe.clone() as Arc<dyn TranscriptPlugin>])
            .build()?;
        let var = var.split(':').collect::<Vec<_>>();
        predictor_with_config(tx, seq, config).predict(&VcfVariant {
            chromosome: chromosome.into(),
            position: var[0].parse()?,
            reference: var[1].into(),
            alternative: var[2].into(),
        })?;

        let seen = std::mem::take(&mut *probe.0.lock().unwrap());
        assert_yaml_snapshot!(seen);
        Ok(())
    }

    /// `Xaa` is an unknown residue, not a stop. A substitution by it can be synonymous, missense
    /// or a stop, so it is a `coding_sequence_variant`, their common SO parent. The CDS is
    /// `ATG GCC AAG CTG TGG TAA` at chr1:1021-1038, or `ATG GCN ...` with an unknown base.
    #[rstest::rstest]
    // `c.8_9delinsNC`: `AAG` to `ANC`, an `N` that the protein reads gives `p.?`
    #[case::n_in_alt("ATGGCCAAGCTGTGGTAA", (1028, "AG", "NC"), false, Some("p.?"))]
    // VEP writes no HGVS for such alleles
    #[case::n_in_alt_vep(
        "ATGGCCAAGCTGTGGTAA",
        (1028, "AG", "NC"),
        true,
        None
    )]
    // `c.5C>A`: `GCN` (Ala) to `GAN` (Asp or Glu)
    #[case::n_in_reference(
        "ATGGCNAAGCTGTGGTAA",
        (1025, "C", "A"),
        false,
        Some("p.(Ala2Xaa)")
    )]
    fn unknown_residue_is_no_stop(
        #[case] cds: &str,
        #[case] (position, reference, alternative): (i32, &str, &str),
        #[case] vep_compatibility_mode: bool,
        #[case] expected_hgvs_p: Option<&str>,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", cds, "CGCATCA");
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(vep_compatibility_mode)
            .build()?;
        let anns = predictor_with_config(tx, seq, config)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (ann.hgvs_p.as_deref(), &ann.consequences),
            (expected_hgvs_p, &vec![Consequence::CodingSequenceVariant])
        );
        Ok(())
    }

    /// SO: a conservative in-frame insertion "inserts one or more codons into the coding
    /// sequence between existing codons", a disruptive one "within an existing codon". A
    /// conservative in-frame deletion "deletes one or more entire codons from the coding sequence
    /// but does not change any remaining codons", a disruptive one "deletes bases from the coding
    /// sequence starting within an existing codon". The CDS is `ATG GCC AAG CTG TGG TAA` at
    /// chr1:1021-1038. VEP terms do not tell the two apart.
    #[rstest::rstest]
    #[case::ins_between_codons(
        "1026:C:CGGG",
        false,
        "c.6_7insGGG",
        Consequence::ConservativeInframeInsertion
    )]
    #[case::ins_behind_start_codon(
        "1023:G:GCCC",
        false,
        "c.3_4insCCC",
        Consequence::ConservativeInframeInsertion
    )]
    #[case::ins_in_front_of_stop_codon(
        "1035:G:GCCC",
        false,
        "c.15_16insCCC",
        Consequence::ConservativeInframeInsertion
    )]
    #[case::ins_within_codon(
        "1027:A:AGGG",
        false,
        "c.7_8insGGG",
        Consequence::DisruptiveInframeInsertion
    )]
    #[case::dup_ending_with_codon(
        "1026:C:CGCC",
        false,
        "c.4_6dup",
        Consequence::ConservativeInframeInsertion
    )]
    #[case::dup_within_codon(
        "1027:A:ACCA",
        false,
        "c.5_7dup",
        Consequence::DisruptiveInframeInsertion
    )]
    #[case::del_of_codon(
        "1026:CAAG:C",
        false,
        "c.7_9del",
        Consequence::ConservativeInframeDeletion
    )]
    #[case::del_within_codon(
        "1027:AAGC:A",
        false,
        "c.8_10del",
        Consequence::DisruptiveInframeDeletion
    )]
    #[case::vep_ins_between_codons(
        "1026:C:CGGG",
        true,
        "c.6_7insGGG",
        Consequence::InframeInsertion
    )]
    // `GCC` to `GCG CAC`: VEP's changed codons read A to AH, which starts with A
    #[case::vep_ins_within_codon("1025:C:CGCA", true, "c.5_6insGCA", Consequence::InframeInsertion)]
    fn conservative_or_disruptive_inframe_indel(
        #[case] var: &str,
        #[case] vep_compatibility_mode: bool,
        #[case] expected_hgvs_c: &str,
        #[case] expected: Consequence,
    ) -> Result<(), anyhow::Error> {
        let var = var.split(':').collect::<Vec<_>>();
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CGCATCA");
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(vep_compatibility_mode)
            .build()?;
        let anns = predictor_with_config(tx, seq, config)
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
        let inframe = Consequence::ConservativeInframeInsertion
            | Consequence::DisruptiveInframeInsertion
            | Consequence::ConservativeInframeDeletion
            | Consequence::DisruptiveInframeDeletion
            | Consequence::InframeInsertion
            | Consequence::InframeDeletion;
        let found = (ann.consequences.iter().copied())
            .filter(|csq| inframe.contains(*csq))
            .collect::<Vec<_>>();

        assert_eq!(
            (ann.hgvs_c.as_deref(), found),
            (Some(expected_hgvs_c), vec![expected]),
            "{:?}",
            ann.consequences
        );
        Ok(())
    }

    /// An insertion at an exon edge whose bases join the exon lies between the CDS bases on both
    /// sides of the edge, see `conservative_or_disruptive_inframe_indel`. Exon 1 ends with
    /// c.9, the end of a codon, and the CDS is `ATG GCC AAG | CTG TGG TAA`.
    #[rstest::rstest]
    #[case::behind_exon_end(1029, "G", "GCCC", "c.9_9+1insCCC")]
    #[case::in_front_of_exon_start(1053, "G", "GAAA", "c.10-1_10insAAA")]
    fn conservative_inframe_insertion_at_exon_edge(
        #[case] position: i32,
        #[case] reference: &str,
        #[case] alternative: &str,
        #[case] expected_hgvs_c: &str,
    ) -> Result<(), anyhow::Error> {
        let exon_1 = "CAGTCAGTCAGTCAGTCAGTATGGCCAAG";
        let intron = format!("GT{}AG", "C".repeat(20));
        let exon_2 = "CTGTGGTAACATCA";
        let exon = |ord, alt_start_i, alt_end_i, tx_start, tx_end| ExonAlignment {
            alt_start_i,
            alt_end_i,
            ord,
            alt_cds_start_i: Some(tx_start),
            alt_cds_end_i: Some(tx_end),
            cigar: format!("{}M", tx_end - tx_start + 1),
        };
        let tx = Transcript {
            id: "NM_000001.1".into(),
            gene_symbol: "GENE1".into(),
            biotype: TranscriptBiotype::Coding.into(),
            protein: Some("NP_000001.1".into()),
            start_codon: Some(20),
            stop_codon: Some(38),
            genome_alignments: vec![GenomeAlignment {
                genome_build: "grch38".into(),
                contig: "NC_000001.11".into(),
                cds_start: Some(1020),
                cds_end: Some(1062),
                strand: Strand::Plus.into(),
                exons: vec![exon(0, 1000, 1029, 1, 29), exon(1, 1053, 1067, 30, 43)],
                ..Default::default()
            }],
            filtered: Some(false),
            ..Default::default()
        };
        let mut genome = NamedTempFile::new()?;
        write!(
            genome,
            ">1\n{}{exon_1}{intron}{exon_2}{}\n",
            "T".repeat(1000),
            "T".repeat(100)
        )?;
        let predictor = predictor_with_reference(
            tx,
            format!("{exon_1}{exon_2}"),
            Config::default(),
            Some(genome.path()),
        );
        let anns = predictor
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.consequences
                    .contains(&Consequence::ConservativeInframeInsertion),
                ann.consequences
                    .contains(&Consequence::DisruptiveInframeInsertion)
            ),
            (Some(expected_hgvs_c), true, false),
            "{:?}",
            ann.consequences
        );
        Ok(())
    }

    /// `Gene_ID` is the gene ID of the transcript.
    #[test]
    fn ann_gene_id() -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGT", "ATGGCCTAA", "CAGT");
        let tx = Transcript {
            gene_id: "11".into(),
            hgnc_id: Some("HGNC:1".into()),
            ..tx
        };
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1008,
                reference: "G".into(),
                alternative: "A".into(),
            })?
            .unwrap_or_default();

        let gene_ids: Vec<_> = anns.iter().map(|ann| ann.gene_id.as_str()).collect();
        assert_eq!(gene_ids, ["11"]);

        Ok(())
    }

    /// Two transcripts and their gene entries. `merged` has a RefSeq and an Ensembl gene of one
    /// HGNC ID. `versions` has two versions of a gene without HGNC ID.
    fn two_transcripts_db(db: &str) -> (Vec<Transcript>, Vec<crate::pbs::txs::GeneToTxId>) {
        use crate::annotate::seqvars::provider::test_db;
        use crate::pbs::txs::GeneToTxId;

        let tx = |id: &str, gene_id: &str, hgnc_id: Option<&str>| Transcript {
            gene_id: gene_id.into(),
            hgnc_id: hgnc_id.map(String::from),
            ..test_db::transcript(id, 1000)
        };
        if db == "merged" {
            let refseq = tx("NM_000001.1", "672", Some("HGNC:1100"));
            let ensembl = tx("ENST00000000001.1", "ENSG00000012048.1", Some("HGNC:1100"));
            let genes = vec![
                test_db::gene(&refseq),
                GeneToTxId {
                    gene_id: "ENSG00000012048".into(),
                    ..test_db::gene(&ensembl)
                },
            ];
            (vec![refseq, ensembl], genes)
        } else {
            let old = tx("ENST00000000001.1", "ENSG00000000001.9", None);
            let new = tx("ENST00000000001.2", "ENSG00000000001.10", None);
            let gene = GeneToTxId {
                gene_id: "ENSG00000000001".into(),
                tx_ids: vec![old.id.clone(), new.id.clone()],
                ..test_db::gene(&new)
            };
            (vec![old, new], vec![gene])
        }
    }

    /// Picking groups the transcripts by gene entry, by length and by severity. A database with a
    /// RefSeq and an Ensembl gene of one HGNC ID thus picks one transcript of each gene, and a
    /// database with two versions of a gene picks one.
    #[rstest::rstest]
    #[case::refseq_and_ensembl_gene("merged", &["NM_000001.1", "ENST00000000001.1"])]
    #[case::two_versions_of_a_gene("versions", &["ENST00000000001.1"])]
    fn picks_per_gene(#[case] db: &str, #[case] expected: &[&str]) -> Result<(), anyhow::Error> {
        use crate::annotate::seqvars::provider::test_db;

        let (transcripts, genes) = two_transcripts_db(db);
        let txs: Vec<_> = transcripts
            .iter()
            .map(|tx| TxForRegionRecord {
                tx_ac: tx.id.clone(),
                alt_ac: "NC_000001.11".into(),
                alt_strand: 1,
                alt_aln_method: ALT_ALN_METHOD.into(),
                start_i: 1000,
                end_i: 1100,
            })
            .collect();
        // The first transcript has the more severe consequence.
        let consequences = transcripts
            .iter()
            .map(|tx| tx.id.clone())
            .zip([Consequence::MissenseVariant, Consequence::SynonymousVariant])
            .collect::<BTreeMap<_, _>>();
        for settings in [
            TranscriptSettings {
                pick: Some(PickUnit::Gene),
                pick_order: vec![PickCriterion::Length],
                ..Default::default()
            },
            TranscriptSettings::most_severe(PickUnit::Gene),
        ] {
            let provider = MehariProvider::new(
                test_db::tx_seq_db(transcripts.clone(), genes.clone()),
                None::<PathBuf>,
                true,
            );
            let predictor = ConsequencePredictor::new(
                Arc::new(provider),
                ConfigBuilder::default()
                    .transcript_settings(settings)
                    .build()?,
            );

            let reported: Vec<_> = predictor
                .select_annotations(txs.clone(), |record| {
                    Ok(Some(TxAnnotation {
                        context: TxContext::empty(),
                        ann: AnnField {
                            feature_id: record.tx_ac.clone(),
                            consequences: vec![consequences[&record.tx_ac]],
                            ..Default::default()
                        },
                    }))
                })?
                .into_iter()
                .map(|tx| tx.ann.feature_id)
                .collect();
            assert_eq!(reported, expected);
        }

        Ok(())
    }

    /// A selenoprotein without Sec positions reads UGA as selenocysteine inside its CDS, but a
    /// UGA at the CDS end is its stop codon, as in SEPHS2.
    ///
    /// The transcript is the one of `selenoprotein_predictor`, but its CDS
    /// `ATG GCC AAG CTG TGG GAA CCA TGA CGC GTT TGA` (`MAKLWEPURV*`) ends with UGA.
    #[rstest::rstest]
    #[case::stop_codon("1051:T:C", "p.(Ter11ArgextTer?)", Consequence::StopLost)]
    #[case::sec_codon("1044:A:G", "p.(Sec8Trp)", Consequence::MissenseVariant)]
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
    #[case::synonymous("1038:A:G", "p.(Glu6=)", &[Consequence::SynonymousVariant])]
    #[case::missense("1038:A:C", "p.(Glu6Asp)", &[Consequence::MissenseVariant])]
    // c.16_18del
    #[case::deletion(
        "1035:GGAA:G",
        "p.(Glu6del)",
        &[Consequence::ConservativeInframeDeletion]
    )]
    // c.17_18insC
    #[case::frameshift("1037:A:AC", "p.(Glu6AspfsTer?)", &[Consequence::FrameshiftVariant])]
    // c.18del
    #[case::frameshift_in_last_codon("1036:GA:G", "p.?", &[Consequence::FrameshiftVariant])]
    // c.16del
    #[case::frameshift_at_last_codon("1034:GG:G", "p.?", &[Consequence::FrameshiftVariant])]
    // c.16_17insGGGGGTAAGG: the new frame reads a stop codon behind the end of the CDS, which
    // is no `feature_elongation`, as the residues behind that end are unknown
    #[case::frameshift_with_stop_behind_end(
        "1035:G:GGGGGGGTAAG",
        "p.(Glu6delinsGlyGlyTer)",
        &[Consequence::FrameshiftVariant]
    )]
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
    #[case::stop_gained("1035:G:A", "p.(Trp5Ter)", &[Consequence::StopGained])]
    #[case::last_complete_codon("1038:A:G", "p.(Glu6=)", &[Consequence::SynonymousVariant])]
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
        "p.(Trp5SerfsTer2)",
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
    #[case::missense("1033:T:C", "p.(Trp5Arg)", Consequence::MissenseVariant)]
    #[case::stop_lost("1039:T:C", "p.(Ter7GlnextTer?)", Consequence::StopLost)]
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
    // VEP 116.3 reports `stop_lost` only, and mehari keeps `frameshift_variant`
    #[case(true, vec![Consequence::FrameshiftVariant, Consequence::StopLost])]
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
                .vep_compatibility_mode(vep_consequence_terms)
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
            (
                Some("c.699del"),
                Some("p.(Ter233XaaextTer?)"),
                &expected_csqs
            )
        );
        Ok(())
    }

    /// An indel crosses the end of a CDS of 30 bases if it starts in the CDS and ends at c.30 or
    /// behind it.
    #[rstest::rstest]
    #[case("c.10_30del", true)]
    #[case("c.10_29del", false)]
    #[case("c.10_*2del", true)]
    #[case("c.30_30+5del", true)]
    #[case("c.-2_30del", false)]
    #[case("c.10+5_*2del", false)]
    fn indel_crosses_cds_end(
        #[case] var_c: &str,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let var_c = HgvsVariant::from_str(&format!("NM_000001.1:{var_c}"))?;
        assert_eq!(crosses_cds_end(&var_c, Some(30)), expected);
        Ok(())
    }

    /// A deletion from the CDS through the stop codon leaves only CDS bases in front of it, so no
    /// CDS codon reads in a new frame. HGVS still describes it as a frameshift, or as `p.?` if it
    /// reaches past the transcript end. The terms come from the protein it produces.
    ///
    /// GLA `ENST00000218516.4` ends 6 bases behind its stop codon. `c.1286_*7del` deletes 12
    /// bases and reaches past the transcript end, so `p.?`, and the new frame has no stop codon.
    /// IRF7 `c.1507_*40del` forms a new stop codon at the position of the old one.
    #[rstest::rstest]
    #[case::gla(
        "X",
        101397801,
        "TAAACATTTTAAA",
        "ENST00000218516.4",
        false,
        ("c.1286_*7del", "p.?"),
        vec![Consequence::StopLost, Consequence::ThreePrimeUtrExonVariant]
    )]
    // VEP 116.3 gives the same terms
    #[case::gla_vep(
        "X",
        101397801,
        "TAAACATTTTAAA",
        "ENST00000218516.4",
        true,
        ("c.1286_*7del", "p.?"),
        vec![Consequence::StopRetainedVariant, Consequence::ThreePrimeUtrVariant]
    )]
    #[case::irf7(
        "11",
        612604,
        "TGGGCTGCTCCAGCTTTCTGGAGTTCTCATTAGACTGGGTTCTAGGC",
        "ENST00000525445.6",
        false,
        ("c.1507_*40del", "p.(Ala503ThrfsTer2)"),
        vec![
            Consequence::MissenseVariant,
            Consequence::StopRetainedVariant,
            Consequence::ThreePrimeUtrExonVariant
        ]
    )]
    // VEP 116.3 gives the same terms
    #[case::irf7_vep(
        "11",
        612604,
        "TGGGCTGCTCCAGCTTTCTGGAGTTCTCATTAGACTGGGTTCTAGGC",
        "ENST00000525445.6",
        true,
        ("c.1507_*40del", "p.(Ala503ThrfsTer2)"),
        vec![Consequence::StopRetainedVariant, Consequence::ThreePrimeUtrVariant]
    )]
    fn deletion_over_stop_codon_without_frameshift(
        #[case] chromosome: &str,
        #[case] position: i32,
        #[case] reference: &str,
        #[case] tx_id: &str,
        #[case] vep_consequence_terms: bool,
        #[case] expected_hgvs: (&str, &str),
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let tx_db = load_tx_db(
            "tests/data/annotate/db/grch38/GRCh38-ensembl.disagreement-subset.txs.bin.zst",
        )?;
        let provider = Arc::new(MehariProvider::new(tx_db, None::<PathBuf>, true));
        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .vep_consequence_terms(vep_consequence_terms)
                .vep_compatibility_mode(vep_consequence_terms)
                .build()?,
        );

        let anns = predictor
            .predict(&VcfVariant {
                chromosome: chromosome.into(),
                position,
                reference: reference.into(),
                alternative: reference[..1].into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == tx_id)
            .ok_or_else(|| anyhow::anyhow!("no annotation for {tx_id}: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (Some(expected_hgvs.0), Some(expected_hgvs.1), &expected_csqs)
        );
        Ok(())
    }

    /// As above, on a test transcript. Its 3' UTR decides where the new frame reads a stop codon.
    #[rstest::rstest]
    // `c.14_*7delinsCCC` replaces 12 bases by 3 and reaches past the transcript end, so `p.?`.
    // The new frame reads `Xaa` in place of the old stop codon.
    #[case::delins(
        "CGCATC",
        (1033, "TGGTAACGCATCA", "TCCC"),
        ("c.14_*7delinsCCC", "p.?"),
        vec![
            Consequence::StopLost,
            Consequence::FeatureElongation,
            Consequence::ThreePrimeUtrExonVariant
        ]
    )]
    // `c.14_*7delinsCC` replaces 12 bases by 2 and reaches past the transcript end, so `p.?`.
    // The new frame has no stop codon and no more residues than the reference.
    #[case::delins_in_frame_in_transcript(
        "CGCATC",
        (1033, "TGGTAACGCATCA", "TCC"),
        ("c.14_*7delinsCC", "p.?"),
        vec![Consequence::StopLost, Consequence::ThreePrimeUtrExonVariant]
    )]
    #[case::new_stop_in_front(
        "ACAGTAAGCAT",
        (1030, "CTGTGGTAAAC", "C"),
        ("c.11_*2del", "p.(Leu4GlnfsTer2)"),
        vec![Consequence::StopGained, Consequence::ThreePrimeUtrExonVariant]
    )]
    #[case::new_stop_behind(
        "ACCGCATTAAGC",
        (1033, "TGGTAAAC", "T"),
        ("c.14_*2del", "p.(Trp5SerfsTer3)"),
        vec![
            Consequence::StopLost,
            Consequence::FeatureElongation,
            Consequence::ThreePrimeUtrExonVariant
        ]
    )]
    // No new stop codon, but the new frame reads past the old one to the transcript end
    #[case::no_new_stop_reading_past_old_one(
        "ACCGCATCAGCGCC",
        (1033, "TGGTAAAC", "T"),
        ("c.14_*2del", "p.(Trp5SerfsTer?)"),
        vec![
            Consequence::StopLost,
            Consequence::FeatureElongation,
            Consequence::ThreePrimeUtrExonVariant
        ]
    )]
    fn indel_over_stop_codon_without_frameshift(
        #[case] utr3: &str,
        #[case] variant: (i32, &str, &str),
        #[case] expected_hgvs: (&str, &str),
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", utr3);
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: variant.0,
                reference: variant.1.into(),
                alternative: variant.2.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (Some(expected_hgvs.0), Some(expected_hgvs.1), &expected_csqs)
        );
        Ok(())
    }

    /// A substitution of the last stop codon base changes no length, so it keeps the terms of p.:
    /// TAA to TAC (Tyr) reads on into the 3' UTR, TAA to TAG keeps the stop codon.
    #[rstest::rstest]
    #[case("C", "c.18A>C", "p.(Ter6TyrextTer?)", vec![Consequence::StopLost, Consequence::FeatureElongation])]
    #[case("G", "c.18A>G", "p.(Ter6=)", vec![Consequence::StopRetainedVariant])]
    fn substitution_of_last_stop_codon_base(
        #[case] alternative: &str,
        #[case] expected_c: &str,
        #[case] expected_p: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CGCATCA");
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1038,
                reference: "A".into(),
                alternative: alternative.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;
        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (Some(expected_c), Some(expected_p), &expected_csqs)
        );
        Ok(())
    }

    /// Only deletions and edits that change the number of bases change the length.
    #[rstest::rstest]
    #[case("c.10_12del", true)]
    #[case("c.10_12delinsA", true)]
    #[case("c.10_12delinsAAA", false)]
    #[case("c.10A>C", false)]
    #[case("c.10_11insA", false)]
    fn edit_changes_length(
        #[case] var_c: &str,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let HgvsVariant::CdsVariant { loc_edit, .. } =
            HgvsVariant::from_str(&format!("NM_000001.1:{var_c}"))?
        else {
            anyhow::bail!("not a c. variant: {var_c}");
        };
        let edit = match loc_edit.edit.inner() {
            // The parser keeps the reference of `delins` empty; give it one as mehari's c. has.
            NaEdit::RefAlt {
                reference,
                alternative,
            } if reference.is_empty() => NaEdit::RefAlt {
                reference: "GCC".into(),
                alternative: alternative.clone(),
            },
            edit => edit.clone(),
        };
        assert_eq!(changes_length(&edit), expected);
        Ok(())
    }

    /// `c.7_10del` shifts the frame of codons 4 and 5 and of the stop codon, and a new stop codon
    /// from the 3' UTR takes the old position. The protein keeps its length, so the frameshift is
    /// also a missense change.
    #[test]
    fn frameshift_keeping_protein_length() -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CTAAGCATCA");
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1026,
                reference: "CAAGC".into(),
                alternative: "C".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (
                Some("c.7_10del"),
                Some("p.(Lys3CysfsTer4)"),
                &vec![
                    Consequence::FrameshiftVariant,
                    Consequence::MissenseVariant,
                    Consequence::StopRetainedVariant
                ]
            )
        );
        Ok(())
    }

    /// Where p. and the protein disagree, the protein decides the terms, and `frameshift_variant`
    /// requires a codon of the CDS behind the change that reads in another frame. On a test
    /// transcript with 20 bases of 5' UTR, the CDS `cds` at chr1:1021 and the 3' UTR `utr3`.
    #[rstest::rstest]
    // `c.11_18del` deletes through the stop codon, and the protein ends at a new stop codon in
    // front of the old position. No CDS codon is left to read in another frame.
    #[case::through_stop_codon(
        ("ATGGCCAAGCTGTGGTAA", "CATAAGCATCA"),
        (1030, "CTGTGGTAA", "C"),
        ("c.11_18del", "p.(Leu4ProfsTer2)"),
        vec![Consequence::StopGained]
    )]
    // `c.14_17del` ends inside the stop codon. A new stop codon from the 3' UTR takes the old
    // position, so only Trp5 changes.
    #[case::into_stop_codon(
        ("ATGGCCAAGCTGTGGTAA", "CTAAGCATCA"),
        (1033, "TGGTA", "T"),
        ("c.14_17del", "p.(Trp5TyrfsTer2)"),
        vec![Consequence::MissenseVariant, Consequence::StopRetainedVariant]
    )]
    // `c.16_17insC` splits the stop codon, and the protein reads on into the 3' UTR.
    #[case::ins_into_stop_codon(
        ("ATGGCCAAGCTGTGGTAA", "CTAAGCATCA"),
        (1036, "T", "TC"),
        ("c.16_17insC", "p.(Ter6SerextTer?)"),
        vec![Consequence::StopLost, Consequence::FeatureElongation]
    )]
    // `c.14_16dup` keeps the frame, so its p. terms stay. Its copy splits the stop codon.
    #[case::in_frame_dup_into_stop_codon(
        ("ATGGCCAAGCTGTGGTAA", "CTAAGCATCA"),
        (1036, "T", "TGGT"),
        ("c.14_16dup", "p.(Ter6TrpextTer1)"),
        vec![
            Consequence::StopLost,
            Consequence::FeatureElongation,
            Consequence::DisruptiveInframeInsertion
        ]
    )]
    // `c.15_16insGAT` adds Asp in front of the stop codon, which stays.
    #[case::in_frame_ins_before_stop_codon(
        ("ATGGCCAAGCTGTGGTAA", "CTAAGCATCA"),
        (1035, "G", "GGAT"),
        ("c.15_16insGAT", "p.(Ter6AspextTer1)"),
        vec![
            Consequence::FeatureElongation,
            Consequence::ConservativeInframeInsertion
        ]
    )]
    // `c.7_10del` shifts the frame of codons 4 and 5 and of the stop codon, and a new stop codon
    // from the 3' UTR takes the old position.
    #[case::fs_keeping_length(
        ("ATGGCCAAGCTGTGGTAA", "CTAAGCATCA"),
        (1026, "CAAGC", "C"),
        ("c.7_10del", "p.(Lys3CysfsTer4)"),
        vec![
            Consequence::FrameshiftVariant,
            Consequence::MissenseVariant,
            Consequence::StopRetainedVariant
        ]
    )]
    // `c.15del` keeps Ala5 (`GCA` to `GCT`), and the stop codon reads in another frame.
    #[case::stop_codon_in_new_frame(
        ("ATGGCCAAGCTGGCATAA", "GCATCAGTAAC"),
        (1034, "CA", "C"),
        ("c.15del", "p.(Ter6LysextTer3)"),
        vec![
            Consequence::FrameshiftVariant,
            Consequence::StopLost,
            Consequence::FeatureElongation
        ]
    )]
    // `c.2del` loses the start codon, so no reading frame of the protein shifts.
    #[case::start_lost(
        ("ATGGCCAAGCTGTGGTAA", "CGCATCA"),
        (1021, "AT", "A"),
        ("c.2del", "p.Met1?"),
        vec![Consequence::StartLost]
    )]
    // `c.2_5del` loses the start codon. Read from c.1, the new frame ends at the old position,
    // but the protein is unknown.
    #[case::start_lost_frame_ending_at_old_position(
        ("ATGGCCAAGCTGTGGTAA", "CTAAGCATCA"),
        (1021, "ATGGC", "A"),
        ("c.2_5del", "p.Met1?"),
        vec![Consequence::StartLost]
    )]
    // `c.9_11delinsCTA` keeps the number of residues in p., but the protein ends at the new stop
    // codon.
    #[case::delins_ending_with_stop(
        ("ATGGCCAAGCTGTGGTAA", "CGCATCA"),
        (1029, "GCT", "CTA"),
        ("c.9_11delinsCTA", "p.(Lys3_Leu4delinsAsnTer)"),
        vec![Consequence::StopGained]
    )]
    // `c.15_17delinsCTG` turns TGG (Trp) into TGC (Cys) and the stop codon TAA into TGA. p. omits
    // the unchanged stop, but its codon changed.
    #[case::delins_keeping_stop_position(
        ("ATGGCCAAGCTGTGGTAA", "CGCATCA"),
        (1035, "GTA", "CTG"),
        ("c.15_17delinsCTG", "p.(Trp5Cys)"),
        vec![Consequence::MissenseVariant, Consequence::StopRetainedVariant]
    )]
    // `c.15delinsCTAA` turns TGG (Trp) into TGC (Cys) and puts a new TAA at the old position.
    // The protein keeps its length.
    #[case::delins_with_stop_at_old_position(
        ("ATGGCCAAGCTGTGGTAA", "CGCATCA"),
        (1035, "G", "CTAA"),
        ("c.15delinsCTAA", "p.(Trp5delinsCysTer)"),
        vec![Consequence::MissenseVariant, Consequence::StopRetainedVariant]
    )]
    #[case::delins_of_two_bases_with_stop_at_old_position(
        ("ATGGCCAAGCTGTGGTAA", "CGCATCA"),
        (1034, "GG", "CCTAA"),
        ("c.14_15delinsCCTAA", "p.(Trp5delinsSerTer)"),
        vec![Consequence::MissenseVariant, Consequence::StopRetainedVariant]
    )]
    // The frameshifts `c.13_14insCCTAA` and `c.14delinsCCTAA` both give `MAKLS*`, with a new TAA
    // at the old position. p. describes one as a substitution, the other as `fs`.
    #[case::fs_ins_keeping_length(
        ("ATGGCCAAGCTGTGGTAA", "CGCATCA"),
        (1033, "T", "TCCTAA"),
        ("c.13_14insCCTAA", "p.(Trp5Ser)"),
        vec![
            Consequence::FrameshiftVariant,
            Consequence::MissenseVariant,
            Consequence::StopRetainedVariant
        ]
    )]
    #[case::fs_delins_keeping_length(
        ("ATGGCCAAGCTGTGGTAA", "CGCATCA"),
        (1034, "G", "CCTAA"),
        ("c.14delinsCCTAA", "p.(Trp5SerfsTer2)"),
        vec![
            Consequence::FrameshiftVariant,
            Consequence::MissenseVariant,
            Consequence::StopRetainedVariant
        ]
    )]
    fn protein_decides_terms(
        #[case] (cds, utr3): (&str, &str),
        #[case] variant: (i32, &str, &str),
        #[case] expected_hgvs: (&str, &str),
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", cds, utr3);
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: variant.0,
                reference: variant.1.into(),
                alternative: variant.2.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (Some(expected_hgvs.0), Some(expected_hgvs.1), &expected_csqs)
        );
        Ok(())
    }

    /// p. writes a stop in place of a residue as a substitution, e.g. `p.(Ala5Ter)`. A deletion
    /// of whole codons in front of the stop codon only shortens the protein, other changes gain
    /// a stop. The CDS is `ATG GCC AAG CTG GCA TAG`.
    #[rstest::rstest]
    #[case::codon_deletion(
        (1032, "GGCA", "G"),
        ("c.13_15del", "p.(Ala5Ter)"),
        vec![Consequence::ConservativeInframeDeletion]
    )]
    #[case::two_codon_deletion(
        (1029, "GCTGGCA", "G"),
        ("c.10_15del", "p.(Leu4Ter)"),
        vec![Consequence::ConservativeInframeDeletion]
    )]
    #[case::delins(
        (1033, "GCA", "TAAC"),
        ("c.13_15delinsTAAC", "p.(Ala5Ter)"),
        vec![Consequence::StopGained, Consequence::FrameshiftVariant]
    )]
    fn stop_in_place_of_residue(
        #[case] (position, reference, alternative): (i32, &str, &str),
        #[case] expected_hgvs: (&str, &str),
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx(
            "CAGTCAGTCAGTCAGTCAGT",
            "ATGGCCAAGCTGGCATAG",
            "GCATCAGTAACTGAGG",
        );
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (Some(expected_hgvs.0), Some(expected_hgvs.1), &expected_csqs)
        );
        Ok(())
    }

    /// With VEP compatibility, `missense_variant` does not go with a stop term, as in VEP's
    /// `missense_variant`. `c.9_11delinsCTA` changes `KL` to `N*`.
    #[test]
    fn vep_missense_without_stop_terms() -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CGCATCA");
        let config = ConfigBuilder::default()
            .vep_consequence_terms(true)
            .vep_compatibility_mode(true)
            .build()?;
        let anns = predictor_with_config(tx, seq, config)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1029,
                reference: "GCT".into(),
                alternative: "CTA".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (Some("c.9_11delinsCTA"), &vec![Consequence::StopGained])
        );
        Ok(())
    }

    /// In VEP compatibility mode, VEP shifts an insertion 3' also across the CDS end, where
    /// hgvs-rs stops. The A inserted behind c.16 is `c.18dup`. If the 3' UTR starts with an A,
    /// VEP shifts it there, outside the CDS. Else it stays behind c.18, and its input position
    /// in the stop codon gives `stop_retained_variant`.
    #[rstest::rstest]
    #[case::into_utr("AAGCATCA", vec![Consequence::ThreePrimeUtrVariant])]
    #[case::at_cds_end(
        "CAGCATCA",
        vec![Consequence::StopRetainedVariant, Consequence::ThreePrimeUtrVariant]
    )]
    fn vep_shifts_insertion_across_cds_end(
        #[case] utr3: &str,
        #[case] expected: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", utr3);
        let config = ConfigBuilder::default()
            .vep_consequence_terms(true)
            .vep_compatibility_mode(true)
            .build()?;
        let anns = predictor_with_config(tx, seq, config)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1036,
                reference: "T".into(),
                alternative: "TA".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (Some("c.18dup"), &expected)
        );
        Ok(())
    }

    /// In VEP compatibility mode, VEP applies the translation exceptions to the reference
    /// codons. Deleting the Sec codon `TGA` (codon 8) is no `stop_lost`.
    #[test]
    fn vep_codons_read_translation_exceptions() -> Result<(), anyhow::Error> {
        use crate::pbs::txs::TranslationException;

        let utr = "CAGTCAGTCAGTCAGTCAGT";
        let (tx, seq) = one_exon_tx(utr, "ATGGCCAAGCTGTGGGAACCATGACGCGTTTAA", utr);
        let tx = Transcript {
            translation_exceptions: vec![TranslationException {
                position: 8,
                amino_acid: "U".into(),
            }],
            ..tx
        };
        let config = ConfigBuilder::default()
            .vep_consequence_terms(true)
            .vep_compatibility_mode(true)
            .build()?;
        let anns = predictor_with_config(tx, seq, config)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1041,
                reference: "ATGA".into(),
                alternative: "A".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (
                Some("c.22_24del"),
                Some("p.(Sec8del)"),
                &vec![Consequence::InframeDeletion]
            )
        );
        Ok(())
    }

    /// In VEP compatibility mode, an in-frame insertion is `inframe_insertion` only if the
    /// peptide of the changed codons, cut behind its first stop, starts or ends with the
    /// reference peptide, as in VEP's `inframe_insertion`. Else it is `protein_altering_variant`,
    /// unless the peptide starts with a stop. The CDS is `ATG GCC AAG CTG TGG TAA`.
    #[rstest::rstest]
    // `AAG` to `ACC CAG`: K to TQ
    #[case::protein_altering(
        (1027, "A", "ACCC"),
        true,
        "p.(Lys3delinsThrGln)",
        vec![Consequence::ProteinAlteringVariant]
    )]
    #[case::protein_altering_default(
        (1027, "A", "ACCC"),
        false,
        "p.(Lys3delinsThrGln)",
        vec![Consequence::DisruptiveInframeInsertion, Consequence::ProteinAlteringVariant]
    )]
    // `TGG` to `TAG AGG`: W to *R
    #[case::stop_first((1033, "T", "TAGA"), true, "p.(Trp5Ter)", vec![Consequence::StopGained])]
    // `CTG` to `CAA CTG`: L to QL
    #[case::ends_with_reference(
        (1030, "C", "CAAC"),
        true,
        "p.(Lys3_Leu4insGln)",
        vec![Consequence::InframeInsertion]
    )]
    // a stop between two codons: VEP's reference peptide is empty
    #[case::stop_between_codons(
        (1029, "G", "GTAG"),
        true,
        "p.(Leu4Ter)",
        vec![Consequence::StopGained, Consequence::InframeInsertion]
    )]
    fn vep_inframe_insertion_terms_from_codons(
        #[case] (position, reference, alternative): (i32, &str, &str),
        #[case] vep_compatibility_mode: bool,
        #[case] expected_hgvs_p: &str,
        #[case] expected: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CGCATCA");
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(vep_compatibility_mode)
            .build()?;
        let anns = predictor_with_config(tx, seq, config)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (ann.hgvs_p.as_deref(), &ann.consequences),
            (Some(expected_hgvs_p), &expected)
        );
        Ok(())
    }

    /// An in-frame deletion of the stop codon can read to the transcript end without a new stop
    /// codon. If the last bases there do not fill a codon, hgvs-rs ends the protein with `Xaa`.
    /// `Xaa` is an unknown amino acid, not a stop. HGVS describes the deletion as a frameshift,
    /// and the protein gets longer.
    #[test]
    fn in_frame_deletion_of_stop_codon_ending_in_xaa() -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CGCATCA");
        let anns = predictor_for(tx, seq)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1032,
                reference: "GTGGTAA".into(),
                alternative: "G".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (
                Some("c.13_18del"),
                Some("p.(Trp5ArgfsTer?)"),
                &vec![
                    Consequence::StopLost,
                    Consequence::FeatureElongation,
                    Consequence::ConservativeInframeDeletion
                ]
            )
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

    /// The reference genome for `tx_db_with_indels_in_exon_cigars`: the genome around the exons,
    /// with the two genome-only bases of exon 1 and without the three transcript-only bases of
    /// exon 2. The transcript ends at n.33, before the `A` bases that the stored sequence
    /// appends.
    fn reference_with_indels_in_exon_cigars() -> Result<NamedTempFile, anyhow::Error> {
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
        Ok(reference)
    }

    /// A phased insertion in front of an SNV on a plus strand transcript. Both edits start at
    /// the same transcript index, and the inserted bases stay in front of the changed base.
    #[test]
    fn annotate_multiple_insertion_before_snv() -> Result<(), anyhow::Error> {
        let reference = reference_with_indels_in_exon_cigars()?;
        let provider = Arc::new(MehariProvider::new(
            tx_db_with_indels_in_exon_cigars(),
            Some(reference.path()),
            true,
        ));
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

        assert_eq!(res[0].hgvs_g.as_deref(), Some("g.2006delinsGGGA"));
        assert_eq!(res[0].hgvs_p.as_deref(), Some("p.(Pro5delinsArgAsp)"));

        Ok(())
    }

    /// A phased insertion in a repeat in front of an SNV in the same repeat. The phased path
    /// annotates the combined change like `predict` does as one variant.
    #[test]
    fn annotate_multiple_insertion_in_repeat_before_snv() -> Result<(), anyhow::Error> {
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
        // A C in front of c.14 and c.14C>A turn CCC (Pro5) into CCAC.
        let phased = predictor
            .predict_multiple(&[var(2005, "C", "CC"), var(2006, "C", "A")])?
            .unwrap();
        // The same change as one variant.
        let single = predictor.predict(&var(2006, "C", "CA"))?.unwrap();

        assert_eq!(phased[0].hgvs_c.as_deref(), Some("c.14_15insA"));
        assert_eq!(phased[0].hgvs_p.as_deref(), Some("p.(Phe6LeufsTer?)"));
        let fields = |ann: &AnnField| {
            (
                ann.hgvs_n.clone(),
                ann.hgvs_c.clone(),
                ann.hgvs_p.clone(),
                ann.consequences.clone(),
                ann.cdna_pos.clone(),
                ann.cds_pos.clone(),
                ann.protein_pos.clone(),
            )
        };
        assert_eq!(fields(&phased[0]), fields(&single[0]));

        Ok(())
    }

    /// Like `predict`, the phased path reads a VCF insertion as behind the first REF base. An
    /// insertion and an SNV at the same VCF position therefore insert behind the changed base.
    #[test]
    fn annotate_multiple_snv_before_insertion_in_repeat() -> Result<(), anyhow::Error> {
        let reference = reference_with_indels_in_exon_cigars()?;
        let provider = Arc::new(MehariProvider::new(
            tx_db_with_indels_in_exon_cigars(),
            Some(reference.path()),
            true,
        ));
        let predictor = ConsequencePredictor::new(provider, Default::default());

        let var = |position, reference: &str, alternative: &str| VcfVariant {
            chromosome: "1".into(),
            position,
            reference: reference.into(),
            alternative: alternative.into(),
        };
        // c.14C>A and a C behind c.14 turn CCC (Pro5) into CACC.
        let phased = predictor
            .predict_multiple(&[var(2006, "C", "CC"), var(2006, "C", "A")])?
            .unwrap();
        // The same change as one variant.
        let single = predictor.predict(&var(2006, "C", "AC"))?.unwrap();

        assert_eq!(phased[0].hgvs_g.as_deref(), Some("g.2005_2006insA"));
        assert_eq!(phased[0].hgvs_c.as_deref(), Some("c.13_14insA"));
        assert_eq!(phased[0].hgvs_p.as_deref(), Some("p.(Pro5HisfsTer?)"));
        let fields = |ann: &AnnField| {
            (
                ann.hgvs_g.clone(),
                ann.hgvs_n.clone(),
                ann.hgvs_c.clone(),
                ann.hgvs_p.clone(),
                ann.consequences.clone(),
                ann.cdna_pos.clone(),
                ann.cds_pos.clone(),
                ann.protein_pos.clone(),
            )
        };
        assert_eq!(fields(&phased[0]), fields(&single[0]));

        Ok(())
    }

    /// An insertion lies between two bases, so it does not overlap an SNV next to it. The
    /// order of the input variants does not change the result.
    #[rstest::rstest]
    #[case(2005)] // insertion in front of the SNV base
    #[case(2006)] // insertion behind the SNV base
    fn annotate_multiple_insertion_and_snv_in_any_order(
        #[case] ins_position: i32,
    ) -> Result<(), anyhow::Error> {
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
        let ins = var(ins_position, "C", "CC");
        let snv = var(2006, "C", "A");
        let ins_first = predictor.predict_multiple(&[ins.clone(), snv.clone()])?;
        let snv_first = predictor.predict_multiple(&[snv, ins])?;

        assert!(ins_first.is_some());
        assert_eq!(snv_first, ins_first);

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
            Some("p.(Pro1859_His1860delinsLeuGln)")
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
        assert_eq!(ann.hgvs_p.as_deref(), Some("p.(His1862Gln)"));

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
        assert_eq!(ann.hgvs_p.as_deref(), Some("p.(His1860ThrfsTer62)"));

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

    /// A variant over a whole exon that reaches past it changes the splice site of the exon on
    /// that side. On the genome, the left side of an exon holds its acceptor on the plus strand
    /// and its donor on the minus strand. The first exon in transcript order has no acceptor, the
    /// last no donor. The exon spans chr1:1001-1100.
    #[rstest::rstest]
    #[case::plus_left_first(Strand::Plus, (990, 1100), Rank { ord: 1, total: 3 }, None)]
    #[case::plus_left_middle(
        Strand::Plus, (990, 1100), Rank { ord: 2, total: 3 }, Some(Consequence::SpliceAcceptorVariant)
    )]
    #[case::plus_left_last(
        Strand::Plus, (990, 1100), Rank { ord: 3, total: 3 }, Some(Consequence::SpliceAcceptorVariant)
    )]
    #[case::plus_right_first(
        Strand::Plus, (1000, 1110), Rank { ord: 1, total: 3 }, Some(Consequence::SpliceDonorVariant)
    )]
    #[case::plus_right_middle(
        Strand::Plus, (1000, 1110), Rank { ord: 2, total: 3 }, Some(Consequence::SpliceDonorVariant)
    )]
    #[case::plus_right_last(Strand::Plus, (1000, 1110), Rank { ord: 3, total: 3 }, None)]
    #[case::minus_left_first(
        Strand::Minus, (990, 1100), Rank { ord: 1, total: 3 }, Some(Consequence::SpliceDonorVariant)
    )]
    #[case::minus_left_middle(
        Strand::Minus, (990, 1100), Rank { ord: 2, total: 3 }, Some(Consequence::SpliceDonorVariant)
    )]
    #[case::minus_left_last(Strand::Minus, (990, 1100), Rank { ord: 3, total: 3 }, None)]
    #[case::minus_right_first(Strand::Minus, (1000, 1110), Rank { ord: 1, total: 3 }, None)]
    #[case::minus_right_middle(
        Strand::Minus,
        (1000, 1110),
        Rank { ord: 2, total: 3 },
        Some(Consequence::SpliceAcceptorVariant)
    )]
    #[case::minus_right_last(
        Strand::Minus,
        (1000, 1110),
        Rank { ord: 3, total: 3 },
        Some(Consequence::SpliceAcceptorVariant)
    )]
    #[case::plus_single(Strand::Plus, (990, 1110), Rank { ord: 1, total: 1 }, None)]
    #[case::minus_single(Strand::Minus, (990, 1110), Rank { ord: 1, total: 1 }, None)]
    fn splice_site_of_covered_exon(
        #[case] strand: Strand,
        #[case] (var_start, var_end): (i32, i32),
        #[case] rank: Rank,
        #[case] expected: Option<Consequence>,
    ) {
        let csqs = ConsequencePredictor::analyze_exonic_variant(
            0, strand, var_start, var_end, 1000, 1100, &rank, false,
        );
        assert_eq!(
            csqs & SPLICE_SITES,
            expected.map_or_else(Consequences::empty, Consequences::from)
        );
    }

    /// The deletion of the whole first exon of a minus-strand transcript and of the bases in
    /// front of it changes the donor of the exon, but no acceptor, as the exon has none. That of
    /// the last exon and of the bases behind it changes the acceptor of the exon. These are
    /// `two_exon_tx(20)` mirrored onto the minus strand, with exon 1 at chr1:1161-1220, exon 2 at
    /// chr1:1001-1060 and the CDS at chr1:1041-1200.
    #[rstest::rstest]
    // chr1:1150-1230, the end of the intron, exon 1 and 10 bases in front of the transcript
    #[case::first_exon(1149, 82, Consequence::SpliceDonorVariant)]
    // chr1:995-1070, 6 bases behind the transcript, exon 2 and the start of the intron
    #[case::last_exon(994, 77, Consequence::SpliceAcceptorVariant)]
    fn deletion_of_minus_strand_terminal_exon(
        #[case] position: i32,
        #[case] ref_len: usize,
        #[case] expected: Consequence,
    ) -> Result<(), anyhow::Error> {
        let (mut tx, seq) = two_exon_tx(20);
        let alignment = &mut tx.genome_alignments[0];
        alignment.strand = Strand::Minus.into();
        alignment.cds_start = Some(1040);
        alignment.cds_end = Some(1200);
        for exon in &mut alignment.exons {
            exon.ord = 1 - exon.ord;
            exon.alt_cds_start_i = Some(1 + 60 * exon.ord);
            exon.alt_cds_end_i = Some(60 + 60 * exon.ord);
        }
        let anns = predictor_with_config(tx, seq, Config::default())
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: "A".repeat(ref_len),
                alternative: "A".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        let found = (ann.consequences.iter())
            .filter(|csq| SPLICE_SITES.contains(**csq))
            .collect_vec();
        assert_eq!(found, [&expected], "{ann:?}");
        Ok(())
    }

    /// An MNV, also an inversion, changes only the bases where REF and ALT differ, and only
    /// these get splice terms. The genome holds `two_exon_tx(20)` with the intron `GTACTA`,
    /// `C` bases and `AACTAG`. Its donor +1 lies at chr1:1061, its acceptor -1 at chr1:1160.
    #[rstest::rstest]
    // donor +5 stays T
    #[case::donor_5th_base_kept(
        (1064, "CTA", "GTT"),
        "c.40+4_40+6delinsGTT",
        vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant],
        vec![Consequence::SpliceDonorRegionVariant]
    )]
    // acceptor -2 and -1 stay AG
    #[case::acceptor_kept(
        (1158, "TAGCCG", "CAGCCA"),
        "c.41-3_43delinsCAGCCA",
        vec![
            Consequence::SpliceRegionVariant,
            Consequence::ExonicSpliceRegionVariant,
            Consequence::SplicePolypyrimidineTractVariant
        ],
        vec![Consequence::SpliceRegionVariant]
    )]
    // the inversion keeps the palindrome `CTAG` from -4 to -1
    #[case::acceptor_kept_by_inversion(
        (1155, "AACTAGCC", "GGCTAGTT"),
        "c.41-6_42inv",
        vec![
            Consequence::SpliceRegionVariant,
            Consequence::ExonicSpliceRegionVariant,
            Consequence::SplicePolypyrimidineTractVariant
        ],
        vec![Consequence::SpliceRegionVariant]
    )]
    fn splice_terms_of_changed_mnv_bases(
        #[case] (position, reference, alternative): (i32, &str, &str),
        #[case] expected_c: &str,
        #[case] expected: Vec<Consequence>,
        #[case] expected_vep: Vec<Consequence>,
        #[values(false, true)] vep_compatibility_mode: bool,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = two_exon_tx(20);
        let intron = format!("GTACTA{}AACTAG", "C".repeat(88));
        let mut genome = NamedTempFile::new()?;
        write!(
            genome,
            ">1\n{}{}{intron}{}{}\n",
            "T".repeat(1000),
            &seq[..60],
            &seq[60..],
            "T".repeat(100)
        )?;
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(vep_compatibility_mode)
            .build()?;
        let anns = predictor_with_reference(tx, seq, config, Some(genome.path()))
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        let found = (ann.consequences.iter())
            .filter(|csq| SPLICE_TERMS.contains(**csq))
            .cloned()
            .collect_vec();
        let expected = if vep_compatibility_mode {
            expected_vep
        } else {
            expected
        };
        assert_eq!(
            (ann.hgvs_c.as_deref(), found),
            (Some(expected_c), expected),
            "{ann:?}"
        );
        Ok(())
    }

    /// VEP terms alone only rename terms. VEP's exclusions of overlapping splice and intron
    /// terms, and its terms of a frameshift, need `vep_compatibility_mode`.
    #[rstest::rstest]
    // donor +1 of `two_exon_tx(20)`
    #[case::donor_1(
        true, (1061, "C", "T"),
        vec![Consequence::SpliceDonorVariant, Consequence::IntronVariant],
        vec![Consequence::SpliceDonorVariant]
    )]
    // donor +5 of `two_exon_tx(20)`
    #[case::donor_5(
        true, (1065, "C", "T"),
        vec![
            Consequence::SpliceDonorFifthBaseVariant,
            Consequence::SpliceRegionVariant,
            Consequence::SpliceDonorRegionVariant,
            Consequence::IntronVariant
        ],
        vec![Consequence::SpliceDonorFifthBaseVariant, Consequence::IntronVariant]
    )]
    // `c.7_10del`, `p.(Lys3CysfsTer4)`: a new stop codon takes the old position, see
    // `frameshift_keeping_protein_length`.
    #[case::fs_keeping_length(
        false, (1026, "CAAGC", "C"),
        vec![
            Consequence::FrameshiftVariant,
            Consequence::MissenseVariant,
            Consequence::StopRetainedVariant
        ],
        vec![Consequence::FrameshiftVariant]
    )]
    fn vep_terms_rename_only(
        #[case] two_exons: bool,
        #[case] (position, reference, alternative): (i32, &str, &str),
        #[case] expected_terms_only: Vec<Consequence>,
        #[case] expected_compatibility: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = if two_exons {
            two_exon_tx(20)
        } else {
            one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CTAAGCATCA")
        };
        let mut found = Vec::new();
        for compatibility in [false, true] {
            let config = ConfigBuilder::default()
                .vep_consequence_terms(true)
                .vep_compatibility_mode(compatibility)
                .build()?;
            let anns = predictor_with_config(tx.clone(), seq.clone(), config)
                .predict(&VcfVariant {
                    chromosome: "1".into(),
                    position,
                    reference: reference.into(),
                    alternative: alternative.into(),
                })?
                .unwrap_or_default();
            let ann = anns
                .iter()
                .find(|ann| ann.feature_id == "NM_000001.1")
                .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;
            found.push(ann.consequences.clone());
        }
        assert_eq!(found, [expected_terms_only, expected_compatibility]);
        Ok(())
    }

    /// VEP compatibility mode implies VEP terms, also without `vep_consequence_terms`. Donor +1
    /// of `two_exon_tx(20)`, see `vep_terms_rename_only`.
    #[test]
    fn vep_compatibility_implies_vep_terms() -> Result<(), anyhow::Error> {
        let (tx, seq) = two_exon_tx(20);
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(true)
            .build()?;
        let anns = predictor_with_config(tx, seq, config)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1061,
                reference: "C".into(),
                alternative: "T".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(ann.consequences, vec![Consequence::SpliceDonorVariant]);
        Ok(())
    }

    /// VEP translates no bases other than A, C, G and T, e.g. the `N` bases of an unknown
    /// stretch that ClinVar inserts. So in VEP compatibility mode, an indel that puts in such
    /// bases gets no terms of the changed amino acids, but `coding_sequence_variant` unless it
    /// shifts the frame. By default, a stop codon in front of the `N` bases ends the protein, as
    /// p. has it. The CDS is `ATG GCC AAG CTG TGG GAA CCA TAA` at chr1:1021-1044.
    #[rstest::rstest]
    // `c.6_7insTAGNNNN`
    #[case::frameshift_with_stop(
        "1026:C:CTAGNNNN",
        Some(Consequence::StopGained),
        vec![Consequence::FrameshiftVariant]
    )]
    // `c.5_6insNNN`
    #[case::inframe(
        "1025:C:CNNN",
        Some(Consequence::DisruptiveInframeInsertion),
        vec![Consequence::CodingSequenceVariant]
    )]
    // `c.6_7insTAGNNN`
    #[case::inframe_with_stop(
        "1026:C:CTAGNNN",
        Some(Consequence::StopGained),
        vec![Consequence::CodingSequenceVariant]
    )]
    // `c.7delinsNNNC`, whose default terms depend on how p. reads the `N` bases
    #[case::delins("1027:A:NNNC", None, vec![Consequence::CodingSequenceVariant])]
    // `c.-1_1insNNN`, whose bases go into the 5' UTR
    #[case::utr(
        "1020:T:TNNN",
        Some(Consequence::FivePrimeUtrExonVariant),
        vec![Consequence::FivePrimeUtrVariant]
    )]
    // `AA>ct` at c.7: VEP reads lower-case bases as A, C, G and T
    #[case::lower_case(
        "1027:AA:ct",
        Some(Consequence::MissenseVariant),
        vec![Consequence::MissenseVariant]
    )]
    fn vep_terms_of_ambiguous_bases(
        #[case] var: &str,
        #[case] expected_default: Option<Consequence>,
        #[case] expected_compatibility: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        let var = var.split(':').collect::<Vec<_>>();
        let mut found = Vec::new();
        for vep_compatibility_mode in [false, true] {
            let utr = "CAGTCAGTCAGTCAGTCAGT";
            let (tx, seq) = one_exon_tx(utr, "ATGGCCAAGCTGTGGGAACCATAA", utr);
            let config = ConfigBuilder::default()
                .vep_compatibility_mode(vep_compatibility_mode)
                .build()?;
            let anns = predictor_with_config(tx, seq, config)
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
            found.push(ann.consequences.clone());
        }

        if let Some(expected_default) = expected_default {
            assert!(found[0].contains(&expected_default), "{:?}", found[0]);
        }
        assert_eq!(found[1], expected_compatibility);
        Ok(())
    }

    /// `hgvs_c`, `hgvs_p` and the terms of an annotation.
    type HgvsAndTerms = (Option<&'static str>, Option<&'static str>, Vec<Consequence>);

    /// An `N` in ALT is one unknown base, also an ALT of `N` alone. By default, the terms of the
    /// position and length stay, c. shows the `N`, and p. is `p.?` where the unknown bases can
    /// change the protein. In VEP compatibility mode, there is no HGVS, as in VEP. A variant
    /// without change gets no annotation. The genome holds the transcript at chr1:1001-1064,
    /// with the CDS `ATG GCC AAG CTG TGG GAA CCA TAA` at chr1:1021-1044, and `T` bases behind
    /// it.
    #[rstest::rstest]
    // no change: hgvs-rs would put it at both ends of the contig
    #[case::no_change("1027:A:A", None, None)]
    // `N` alone: one unknown base in place of REF, also of a longer one
    #[case::alt_n(
        "1027:A:N",
        Some((Some("c.7A>N"), Some("p.?"), vec![Consequence::CodingSequenceVariant])),
        Some((None, None, vec![Consequence::CodingSequenceVariant]))
    )]
    #[case::alt_n_longer_ref(
        "1027:AAG:N",
        Some((Some("c.7_9delinsN"), Some("p.?"), vec![Consequence::FrameshiftVariant])),
        Some((None, None, vec![Consequence::FrameshiftVariant]))
    )]
    #[case::alt_n_inframe(
        "1027:AAGC:N",
        Some((Some("c.7_10delinsN"), Some("p.?"), vec![Consequence::CodingSequenceVariant])),
        Some((None, None, vec![Consequence::CodingSequenceVariant]))
    )]
    // one unknown base behind c.6: the frameshift reads it at once
    #[case::frameshift(
        "1026:C:CN",
        Some((Some("c.6_7insN"), Some("p.?"), vec![Consequence::FrameshiftVariant])),
        Some((None, None, vec![Consequence::FrameshiftVariant]))
    )]
    #[case::inframe(
        "1025:C:CNNN",
        Some((Some("c.5_6insNNN"), Some("p.?"), vec![Consequence::DisruptiveInframeInsertion])),
        Some((None, None, vec![Consequence::CodingSequenceVariant]))
    )]
    // `TAG` ends the protein in front of the unknown bases
    #[case::stop_in_front(
        "1026:C:CTAGNNNN",
        Some((
            Some("c.6_7insTAGNNNN"),
            Some("p.(Lys3Ter)"),
            vec![Consequence::StopGained, Consequence::FrameshiftVariant]
        )),
        Some((None, None, vec![Consequence::FrameshiftVariant]))
    )]
    #[case::substitution(
        "1027:AA:NA",
        Some((Some("c.7A>N"), Some("p.?"), vec![Consequence::CodingSequenceVariant])),
        Some((None, None, vec![Consequence::CodingSequenceVariant]))
    )]
    #[case::utr(
        "1020:T:TNNN",
        Some((Some("c.-1_1insNNN"), Some("p.?"), vec![Consequence::FivePrimeUtrExonVariant])),
        Some((None, None, vec![Consequence::FivePrimeUtrVariant]))
    )]
    fn unknown_bases_in_alt(
        #[case] var: &str,
        #[case] expected_default: Option<HgvsAndTerms>,
        #[case] expected_compatibility: Option<HgvsAndTerms>,
    ) -> Result<(), anyhow::Error> {
        let var = var.split(':').collect::<Vec<_>>();
        let utr = "CAGTCAGTCAGTCAGTCAGT";
        let (tx, seq) = one_exon_tx(utr, "ATGGCCAAGCTGTGGGAACCATAA", utr);
        let mut reference = NamedTempFile::new()?;
        write!(
            reference,
            ">1\n{}{seq}{}\n",
            "T".repeat(1000),
            "T".repeat(100)
        )?;
        for (vep_compatibility_mode, expected) in
            [(false, expected_default), (true, expected_compatibility)]
        {
            let config = ConfigBuilder::default()
                .vep_compatibility_mode(vep_compatibility_mode)
                .build()?;
            let anns =
                predictor_with_reference(tx.clone(), seq.clone(), config, Some(reference.path()))
                    .predict(&VcfVariant {
                        chromosome: "1".into(),
                        position: var[0].parse()?,
                        reference: var[1].into(),
                        alternative: var[2].into(),
                    })?
                    .unwrap_or_default();
            let found = anns
                .iter()
                .find(|ann| ann.feature_id == "NM_000001.1")
                .map(|ann| {
                    (
                        ann.hgvs_c.as_deref(),
                        ann.hgvs_p.as_deref(),
                        ann.consequences.clone(),
                    )
                });
            assert_eq!(
                found, expected,
                "compatibility mode {vep_compatibility_mode}: {anns:?}"
            );
            // Without HGVS, there is neither n. nor g.
            if let Some(ann) = anns.iter().find(|ann| ann.feature_id == "NM_000001.1")
                && ann.hgvs_c.is_none()
            {
                assert_eq!((&ann.hgvs_n, &ann.hgvs_g), (&None, &None), "{ann:?}");
            }
        }
        Ok(())
    }

    /// An unknown base in place of the last base of exon 1 of `two_exon_tx(20)`, c.40, keeps
    /// the term of its position and gets `coding_sequence_variant` in place of the terms that
    /// need the amino acid.
    #[test]
    fn unknown_base_at_exon_edge() -> Result<(), anyhow::Error> {
        let (tx, seq) = two_exon_tx(20);
        let anns = predictor_with_config(tx, seq, Config::default())
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1059,
                reference: "CG".into(),
                alternative: "CN".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (
                Some("c.40G>N"),
                Some("p.?"),
                &vec![
                    Consequence::ExonicSpliceRegionVariant,
                    Consequence::CodingSequenceVariant
                ]
            )
        );
        Ok(())
    }

    /// The bases of `c.-1_1ins` go into the 5' UTR, also in front of a transcript without 5'
    /// UTR, where hgvs-rs counts them to the CDS. Unknown bases there leave the CDS as it is, so
    /// they give no `coding_sequence_variant`, e.g. those that ClinVar 2031467 inserts in front
    /// of ENST00000552914.5. The CDS `ATG GCC AAG CTG TGG GAA CCA TAA` starts the transcript at
    /// chr1:1001.
    #[test]
    fn unknown_bases_in_front_of_transcript() -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("", "ATGGCCAAGCTGTGGGAACCATAA", "CAGTCAGTCAGTCAGTCAGT");
        let mut reference = NamedTempFile::new()?;
        write!(
            reference,
            ">1\n{}{seq}{}\n",
            "T".repeat(1000),
            "T".repeat(100)
        )?;
        let config = ConfigBuilder::default().build()?;
        let anns = predictor_with_reference(tx, seq, config, Some(reference.path()))
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1000,
                reference: "T".into(),
                alternative: "TNNNN".into(),
            })?
            .unwrap_or_default();
        let found = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .map(|ann| {
                (
                    ann.hgvs_c.as_deref(),
                    ann.hgvs_p.as_deref(),
                    ann.consequences.clone(),
                )
            });

        assert_eq!(
            found,
            Some((
                Some("c.-1_1insNNNN"),
                Some("p.?"),
                vec![Consequence::FivePrimeUtrExonVariant]
            )),
            "{anns:?}"
        );
        Ok(())
    }

    /// An intergenic variant with `keep_intergenic` follows the rules of `unknown_bases_in_alt`
    /// for its g. description: no HGVS for an unknown base in VEP compatibility mode. The
    /// transcript at chr1:1001-1064 lies farther than `PADDING` from chr1:15000.
    #[rstest::rstest]
    #[case::insertion("TN", Some(Some("g.15000_15001insN")), Some(None))]
    #[case::alt_n("N", Some(Some("g.15000T>N")), Some(None))]
    fn intergenic_unknown_bases(
        #[case] alternative: &str,
        #[case] expected_default: Option<Option<&str>>,
        #[case] expected_compatibility: Option<Option<&str>>,
    ) -> Result<(), anyhow::Error> {
        let utr = "CAGTCAGTCAGTCAGTCAGT";
        let (tx, seq) = one_exon_tx(utr, "ATGGCCAAGCTGTGGGAACCATAA", utr);
        let mut reference = NamedTempFile::new()?;
        write!(
            reference,
            ">1\n{}{seq}{}\n",
            "T".repeat(1000),
            "T".repeat(20000)
        )?;
        for (vep_compatibility_mode, expected) in
            [(false, expected_default), (true, expected_compatibility)]
        {
            let config = ConfigBuilder::default()
                .vep_compatibility_mode(vep_compatibility_mode)
                .keep_intergenic(true)
                .build()?;
            let anns =
                predictor_with_reference(tx.clone(), seq.clone(), config, Some(reference.path()))
                    .predict(&VcfVariant {
                        chromosome: "1".into(),
                        position: 15000,
                        reference: "T".into(),
                        alternative: alternative.into(),
                    })?
                    .unwrap_or_default();
            let found = anns
                .iter()
                .find(|ann| ann.consequences == [Consequence::IntergenicVariant])
                .map(|ann| ann.hgvs_g.as_deref());
            assert_eq!(
                found, expected,
                "compatibility mode {vep_compatibility_mode}: {anns:?}"
            );
        }
        Ok(())
    }

    /// Phased variants follow the rules of `unknown_bases_in_alt` for the combined change, see
    /// there for the transcript.
    #[rstest::rstest]
    // one unknown base behind c.6, and c.10C>T
    #[case::frameshift(
        ("1026:C:CN", "1030:C:T"),
        (Some("c.7_10delinsNAAGT"), Some("p.?"), vec![Consequence::FrameshiftVariant]),
        (None, None, vec![Consequence::FrameshiftVariant])
    )]
    // `N` alone at c.7, and c.10C>T
    #[case::alt_n(
        ("1027:A:N", "1030:C:T"),
        (Some("c.7_10delinsNAGT"), Some("p.?"), vec![Consequence::CodingSequenceVariant]),
        (None, None, vec![Consequence::CodingSequenceVariant])
    )]
    fn phased_unknown_bases(
        #[case] vars: (&str, &str),
        #[case] expected_default: HgvsAndTerms,
        #[case] expected_compatibility: HgvsAndTerms,
    ) -> Result<(), anyhow::Error> {
        let utr = "CAGTCAGTCAGTCAGTCAGT";
        let (tx, seq) = one_exon_tx(utr, "ATGGCCAAGCTGTGGGAACCATAA", utr);
        let mut reference = NamedTempFile::new()?;
        write!(
            reference,
            ">1\n{}{seq}{}\n",
            "T".repeat(1000),
            "T".repeat(100)
        )?;
        let vars = [vars.0, vars.1].map(|var| {
            let var = var.split(':').collect::<Vec<_>>();
            VcfVariant {
                chromosome: "1".into(),
                position: var[0].parse().unwrap_or_default(),
                reference: var[1].into(),
                alternative: var[2].into(),
            }
        });
        for (vep_compatibility_mode, expected) in
            [(false, expected_default), (true, expected_compatibility)]
        {
            let config = ConfigBuilder::default()
                .vep_compatibility_mode(vep_compatibility_mode)
                .build()?;
            let anns =
                predictor_with_reference(tx.clone(), seq.clone(), config, Some(reference.path()))
                    .predict_multiple(&vars)?
                    .unwrap_or_default();
            let ann = anns
                .iter()
                .find(|ann| ann.feature_id == "NM_000001.1")
                .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;
            assert_eq!(
                (
                    ann.hgvs_c.as_deref(),
                    ann.hgvs_p.as_deref(),
                    ann.consequences.clone()
                ),
                expected,
                "compatibility mode {vep_compatibility_mode}: {ann:?}"
            );
            if ann.hgvs_c.is_none() {
                assert_eq!((&ann.hgvs_n, &ann.hgvs_g), (&None, &None), "{ann:?}");
            }
        }
        Ok(())
    }

    /// A dup of one base covers that base, also the first transcript base, `n.1dup`. Its rank and
    /// positions are those of the base.
    #[test]
    fn dup_of_first_transcript_base() -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CGCATCA");
        // The genome holds the transcript at chr1:1001, so that mehari writes the C behind the
        // first base as `g.1001dup`.
        let mut reference = NamedTempFile::new()?;
        write!(
            reference,
            ">1\n{}{seq}{}\n",
            "T".repeat(1000),
            "T".repeat(100)
        )?;
        let predictor =
            predictor_with_reference(tx, seq, Config::default(), Some(reference.path()));
        let anns = predictor
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1001,
                reference: "C".into(),
                alternative: "CC".into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_g.as_deref(),
                ann.hgvs_c.as_deref(),
                &ann.consequences,
                &ann.rank,
                ann.distance,
                &ann.cdna_pos,
            ),
            (
                Some("g.1001dup"),
                Some("c.-20dup"),
                &vec![Consequence::FivePrimeUtrExonVariant],
                &Some(Rank { ord: 1, total: 1 }),
                Some(0),
                &Some(Pos {
                    ord: 1,
                    total: Some(45)
                }),
            )
        );
        Ok(())
    }

    /// The transcripts come from `PADDING` around the variant at both ends of its repeat. A
    /// transcript farther than `PADDING` from the 3'-shifted bases gets no annotation. With
    /// `keep_intergenic`, a variant without a transcript in reach is intergenic. The genome holds
    /// the transcript at chr1:1001-1017, `T` bases behind it and then `AGGG`.
    #[rstest::rstest]
    // `g.6020_6021insAAGG`, 5003 bases behind the transcript; 5000 bases before the 3' shift
    #[case::ins_shifted_beyond(5000, 6017, "TAGGA", None)]
    // `g.6017_6018insAAGG`, 5000 bases behind the transcript
    #[case::ins_shifted_to_limit(4997, 6014, "TAGGA", Some(vec![Consequence::DownstreamGeneVariant]))]
    // `g.6015_6018dup` starts 4998 bases behind the transcript, but its copy goes 5001 bases
    // behind it
    #[case::dup_copy_beyond(4997, 6014, "TAGGG", None)]
    fn no_annotation_beyond_padding(
        #[case] t_bases: usize,
        #[case] position: i32,
        #[case] alternative: &str,
        #[case] expected: Option<Vec<Consequence>>,
        #[values(false, true)] vep_compatibility_mode: bool,
        #[values(false, true)] keep_intergenic: bool,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGT", "ATGGCCTAA", "CAGT");
        let mut reference = NamedTempFile::new()?;
        write!(
            reference,
            ">1\n{}{seq}{}AGGG{}\n",
            "T".repeat(1000),
            "T".repeat(t_bases),
            "C".repeat(100)
        )?;
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(vep_compatibility_mode)
            .keep_intergenic(keep_intergenic)
            .build()?;
        let anns = predictor_with_reference(tx, seq, config, Some(reference.path()))
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: "T".into(),
                alternative: alternative.into(),
            })?
            .unwrap_or_default();
        let found = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .map(|ann| ann.consequences.clone());
        let intergenic =
            (anns.iter()).any(|ann| ann.consequences == [Consequence::IntergenicVariant]);

        assert_eq!(found, expected, "{anns:?}");
        assert_eq!(
            intergenic,
            keep_intergenic && expected.is_none(),
            "{anns:?}"
        );
        Ok(())
    }

    /// A dup whose bases overlap the transcript keeps its annotation, also if its copy goes
    /// farther than `PADDING` behind the transcript. The dup of chr1:1011-6020 covers the stop
    /// codon, the 3' UTR and the 5003 `T` bases behind the transcript.
    #[rstest::rstest]
    fn dup_over_transcript_end_keeps_annotation(
        #[values(false, true)] vep_compatibility_mode: bool,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = one_exon_tx("CAGT", "ATGGCCTAA", "CAGT");
        let t_bases = "T".repeat(5003);
        let mut reference = NamedTempFile::new()?;
        write!(
            reference,
            ">1\n{}{seq}{t_bases}AGGG{}\n",
            "T".repeat(1000),
            "C".repeat(100)
        )?;
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(vep_compatibility_mode)
            .build()?;
        let anns = predictor_with_reference(tx, seq, config, Some(reference.path()))
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position: 1010,
                reference: "C".into(),
                alternative: format!("CTAACAGT{t_bases}"),
            })?
            .unwrap_or_default();
        let found = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .map(|ann| (ann.hgvs_c.as_deref(), ann.consequences.clone()));

        assert_eq!(
            found,
            Some((
                Some("c.7_*5007dup"),
                vec![Consequence::CodingSequenceVariant]
            )),
            "{anns:?}"
        );
        Ok(())
    }

    /// A dup whose bases span an exon edge equals the insertion of its rotated copy at that
    /// edge, e.g. `c.186-2_186dup` of `AGT` and `c.186-1_186insTAG`. hgvs-rs puts it in the
    /// intron, so its terms, rank and positions are those of the intron, in sync with `p.?`. On
    /// the donor side, the copy behind the dup lies in the intron, too, and keeps its terms. The
    /// distance to the exon counts on the genome, so only its size is checked.
    #[rstest::rstest]
    // acceptor side, plus strand
    #[case(
        "22:28773722:C:CAGT",
        "ENST00000249064",
        "c.186-2_186dup",
        vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant],
        Rank { ord: 1, total: 4 }
    )]
    // acceptor side, minus strand
    #[case(
        "22:50246091:A:ATCT",
        "ENST00000216271",
        "c.1651-2_1651dup",
        vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant],
        Rank { ord: 17, total: 19 }
    )]
    // donor side, plus strand
    #[case(
        "22:28773776:T:TGAG",
        "ENST00000249064",
        "c.238_239+1dup",
        vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant],
        Rank { ord: 2, total: 4 }
    )]
    // donor side, minus strand
    #[case(
        "22:26481720:A:ACC",
        "ENST00000336873",
        "c.41_41+1dup",
        vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant],
        Rank { ord: 2, total: 13 }
    )]
    fn dup_across_exon_edge(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_c: &str,
        #[case] expected_csqs: Vec<Consequence>,
        #[case] expected_rank: Rank,
    ) -> Result<(), anyhow::Error> {
        let ann = annotate_chr22_window(spdi, tx_id, false)?;
        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences,
                &ann.rank,
                ann.distance.map(i32::abs),
                &ann.cdna_pos,
                &ann.cds_pos,
            ),
            (
                Some(expected_c),
                Some("p.?"),
                &expected_csqs,
                &Some(expected_rank),
                Some(1),
                &None,
                &None,
            )
        );
        Ok(())
    }

    /// An insertion at either end of a transcript joins its terminal exon on both strands, as
    /// one at an exon edge where no splice site decides. Its HGVS description is reported, in
    /// sync with the terms. In VEP compatibility mode, it lies outside, as in VEP 116.3.
    ///
    /// `one_exon_tx` spans chr1:1001-1045, also mirrored on the minus strand.
    #[rstest::rstest]
    #[case::plus_3p_end(
        Strand::Plus,
        (1045, "A", "AC"),
        false,
        Some("c.*7_*8insC"),
        Consequence::ThreePrimeUtrExonVariant
    )]
    #[case::plus_5p_end(
        Strand::Plus,
        (1000, "T", "TG"),
        false,
        Some("c.-21_-20insG"),
        Consequence::FivePrimeUtrExonVariant
    )]
    #[case::minus_3p_end(
        Strand::Minus,
        (1000, "T", "TC"),
        false,
        Some("c.*7_*8insG"),
        Consequence::ThreePrimeUtrExonVariant
    )]
    #[case::minus_5p_end(
        Strand::Minus,
        (1045, "G", "GA"),
        false,
        Some("c.-21_-20insT"),
        Consequence::FivePrimeUtrExonVariant
    )]
    #[case::plus_3p_end_vep(
        Strand::Plus,
        (1045, "A", "AC"),
        true,
        None,
        Consequence::DownstreamGeneVariant
    )]
    #[case::plus_5p_end_vep(
        Strand::Plus,
        (1000, "T", "TG"),
        true,
        None,
        Consequence::UpstreamGeneVariant
    )]
    #[case::minus_3p_end_vep(
        Strand::Minus,
        (1000, "T", "TC"),
        true,
        None,
        Consequence::DownstreamGeneVariant
    )]
    #[case::minus_5p_end_vep(
        Strand::Minus,
        (1045, "G", "GA"),
        true,
        None,
        Consequence::UpstreamGeneVariant
    )]
    fn insertion_at_transcript_end(
        #[case] strand: Strand,
        #[case] (position, reference, alternative): (i32, &str, &str),
        #[case] vep_compatibility_mode: bool,
        #[case] expected_c: Option<&str>,
        #[case] expected_csq: Consequence,
    ) -> Result<(), anyhow::Error> {
        let (mut tx, seq) = one_exon_tx("CAGTCAGTCAGTCAGTCAGT", "ATGGCCAAGCTGTGGTAA", "CGCATCA");
        if strand == Strand::Minus {
            // n.1 lies at chr1:1045, the CDS at chr1:1008-1025.
            let alignment = &mut tx.genome_alignments[0];
            alignment.strand = Strand::Minus.into();
            alignment.cds_start = Some(1007);
            alignment.cds_end = Some(1025);
        }
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(vep_compatibility_mode)
            .build()?;
        let anns = predictor_with_config(tx, seq, config)
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (ann.hgvs_c.as_deref(), &ann.consequences),
            (expected_c, &vec![expected_csq])
        );
        Ok(())
    }

    /// The insertions behind the last base of ENST00000613577 (plus strand) and of
    /// ENST00000381982 (minus strand) join the last exon, see `insertion_at_transcript_end`.
    #[rstest::rstest]
    #[case("22:18530459:A:AC", "ENST00000613577", "c.1071_*1insC")]
    #[case("22:30489395:A:AC", "ENST00000381982", "c.*481_*482insG")]
    fn insertion_behind_transcript_end(
        #[case] spdi: &str,
        #[case] tx_id: &str,
        #[case] expected_c: &str,
    ) -> Result<(), anyhow::Error> {
        let ann =
            annotate_chr22_window_with_config(spdi, tx_id, Config::default(), CHR22_CDS_END_TXS)?;
        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.hgvs_p.as_deref(),
                &ann.consequences
            ),
            (
                Some(expected_c),
                Some("p.?"),
                &vec![Consequence::ThreePrimeUtrExonVariant]
            )
        );
        Ok(())
    }

    /// A dup of one base at an exon edge covers that base. On ENST00000491919 (minus strand,
    /// non-coding), `n.403dup` is a dup of the last base of exon 2, its first base on the genome.
    /// It lies in exon 2, not in the intron behind it.
    #[test]
    fn dup_of_exon_edge_base() -> Result<(), anyhow::Error> {
        let ann = annotate_chr22_window_with_config(
            "22:28710005:C:CA",
            "ENST00000491919",
            Config::default(),
            CHR22_GENCODE42_TXS,
        )?;

        assert_eq!(
            (
                ann.hgvs_n.as_deref(),
                &ann.consequences,
                &ann.rank,
                ann.distance,
                ann.cdna_pos.map(|pos| pos.ord),
            ),
            (
                Some("n.403dup"),
                &vec![
                    Consequence::ExonicSpliceRegionVariant,
                    Consequence::NonCodingTranscriptExonVariant
                ],
                &Some(Rank { ord: 2, total: 5 }),
                Some(0),
                Some(403),
            )
        );
        Ok(())
    }

    /// VEP checks whether an indel overlaps the stop codon at its input position on the genome,
    /// in the last 3 genomic bases of the CDS (`_overlaps_stop_codon_cil`). An intron splits the
    /// stop codon `TA|G` of this transcript, so these are the last 2 intron bases and c.18.
    /// Deleting c.17 lies outside them. Deleting the intron G in front of c.18 lies inside, also
    /// though it is `c.18del` 3'-shifted. Both leave a stop codon at the old position.
    ///
    /// Exon 1 spans chr1:1001-1037 (20 bases of 5' UTR, c.1 to c.17), the intron
    /// chr1:1038-1061 and exon 2 chr1:1062-1068 (c.18 and 6 bases of 3' UTR).
    #[rstest::rstest]
    #[case::cds_base_in_front_of_intron((1036, "TA", "T"), "c.17del", false)]
    #[case::intron_base_in_window((1060, "AG", "A"), "c.18del", true)]
    fn vep_stop_codon_window_on_genome(
        #[case] (position, reference, alternative): (i32, &str, &str),
        #[case] expected_c: &str,
        #[case] expected_stop_retained: bool,
    ) -> Result<(), anyhow::Error> {
        let exon_1 = "CAGTCAGTCAGTCAGTCAGTATGGCCAAGCTGTGGTA";
        let intron = format!("GT{}AG", "C".repeat(20));
        let exon_2 = "GACATCA";
        let exon = |ord, alt_start_i, alt_end_i, tx_start, tx_end| ExonAlignment {
            alt_start_i,
            alt_end_i,
            ord,
            alt_cds_start_i: Some(tx_start),
            alt_cds_end_i: Some(tx_end),
            cigar: format!("{}M", tx_end - tx_start + 1),
        };
        let tx = Transcript {
            id: "NM_000001.1".into(),
            gene_symbol: "GENE1".into(),
            biotype: TranscriptBiotype::Coding.into(),
            protein: Some("NP_000001.1".into()),
            start_codon: Some(20),
            stop_codon: Some(38),
            genome_alignments: vec![GenomeAlignment {
                genome_build: "grch38".into(),
                contig: "NC_000001.11".into(),
                cds_start: Some(1020),
                cds_end: Some(1062),
                strand: Strand::Plus.into(),
                exons: vec![exon(0, 1000, 1037, 1, 37), exon(1, 1061, 1068, 38, 44)],
                ..Default::default()
            }],
            filtered: Some(false),
            ..Default::default()
        };
        let mut genome = NamedTempFile::new()?;
        write!(
            genome,
            ">1\n{}{exon_1}{intron}{exon_2}{}\n",
            "T".repeat(1000),
            "T".repeat(100)
        )?;
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(true)
            .build()?;
        let predictor =
            predictor_with_reference(tx, format!("{exon_1}{exon_2}"), config, Some(genome.path()));
        let anns = predictor
            .predict(&VcfVariant {
                chromosome: "1".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            })?
            .unwrap_or_default();
        let ann = anns
            .iter()
            .find(|ann| ann.feature_id == "NM_000001.1")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_000001.1: {anns:?}"))?;

        assert_eq!(
            (
                ann.hgvs_c.as_deref(),
                ann.consequences.contains(&Consequence::StopRetainedVariant)
            ),
            (Some(expected_c), expected_stop_retained),
            "{:?}",
            ann.consequences
        );
        Ok(())
    }

    /// VEP's stop codon window on the genome, for a CDS at chr1:1021-1038. On the plus strand,
    /// it is chr1:1036-1038, and an insertion reaches its length to the right. On the minus
    /// strand, it is chr1:1021-1023, and an insertion reaches its length to the left.
    #[rstest::rstest]
    #[case::plus_del_in_front(Strand::Plus, "g.1035del", false)]
    #[case::plus_del_first_stop_base(Strand::Plus, "g.1036del", true)]
    #[case::plus_del_behind(Strand::Plus, "g.1039del", false)]
    #[case::plus_ins_in_front(Strand::Plus, "g.1035_1036insA", false)]
    #[case::plus_ins_fitting(Strand::Plus, "g.1036_1037insA", true)]
    #[case::plus_ins_too_long(Strand::Plus, "g.1036_1037insAA", false)]
    #[case::plus_ins_in_front_of_last_base(Strand::Plus, "g.1037_1038insA", false)]
    #[case::minus_ins_in_front_of_last_base(Strand::Minus, "g.1021_1022insA", true)]
    #[case::minus_ins_reaching_stop_codon(Strand::Minus, "g.1023_1024insAAAA", true)]
    #[case::minus_ins_in_front(Strand::Minus, "g.1024_1025insA", false)]
    #[case::minus_ins_behind(Strand::Minus, "g.1020_1021insA", false)]
    fn vep_stop_codon_window(
        #[case] strand: Strand,
        #[case] var_g: &str,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let alignment = GenomeAlignment {
            cds_start: Some(1020),
            cds_end: Some(1038),
            strand: strand.into(),
            ..Default::default()
        };
        let var_g = HgvsVariant::from_str(&format!("NC_000001.11:{var_g}"))?;
        assert_eq!(vep_overlaps_stop_codon(&var_g, &alignment), expected);
        Ok(())
    }

    /// hgvs-rs takes the intron bases of a dup next to an exon from the variant. Without them,
    /// the protein change is unknown. The region then follows `ins_at_boundary_is_intronic`,
    /// with or without the decision by sequence that VEP compatibility mode switches off.
    /// mehari sets it to `false` in both modes, so the copy joins the exon.
    #[rstest::rstest]
    fn safe_project_c_to_p_dup_without_intron_bases(
        #[values(false, true)] vep_compatibility_mode: bool,
    ) -> Result<(), anyhow::Error> {
        let (tx, seq) = two_exon_tx(20);
        let config = ConfigBuilder::default()
            .vep_compatibility_mode(vep_compatibility_mode)
            .build()?;
        let var_c = HgvsVariant::from_str("NM_000001.1:c.40+1dup")?;
        let projection = predictor_with_config(tx, seq, config).safe_project_c_to_p(&var_c)?;
        assert_eq!(
            projection.map(|(var_p, region)| (FormattedLoc(&var_p).to_string(), region)),
            Some(("p.?".to_string(), Some(VariantLocation::Cds)))
        );
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

    /// Without a region, a variant counts as UTR variant only if both ends of c. lie in one UTR.
    #[rstest::rstest]
    #[case("c.-5_-1del", true)]
    #[case("c.-12+1_-12+5del", true)]
    #[case("c.*1_*5del", true)]
    #[case("c.-5_3del", false)]
    #[case("c.10_*2del", false)]
    #[case("c.-5_*2del", false)]
    fn variant_lies_in_one_utr(
        #[case] var_c: &str,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let var_c = HgvsVariant::from_str(&format!("NM_000001.1:{var_c}"))?;
        assert_eq!(lies_in_one_utr(&var_c), expected);
        Ok(())
    }
}
