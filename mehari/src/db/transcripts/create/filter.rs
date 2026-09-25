use crate::common::progress::{Progress, Unit};
use crate::db::transcripts::create::models::{
    Fix, GeneId, Identifier, Reason, TranscriptExt, TranscriptId, TranscriptLoader,
};
use crate::db::transcripts::create::reference::SequenceProvider;
use anyhow::Error;
use derive_new::new;
use enumflags2::{BitFlag, BitFlags};
use hgvs::data::cdot::json::models::{BioType, Gene, GenomeAlignment, Strand, Tag, Transcript};
use hgvs::sequences::{TranslationTable, translate_cds};
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use once_cell::sync::Lazy;
use rayon::iter::Either;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use seqrepo::AliasOrSeqId;
use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::time::Instant;

pub(crate) fn filter_initial_gene_id_entries(loader: &mut TranscriptLoader) -> Result<(), Error> {
    if loader.disable_filters {
        return Ok(());
    }
    tracing::info!("Filtering Gene Id entries …");
    for (gene_id, tx_ids) in &loader.gene_id_to_transcript_ids {
        if tx_ids.is_empty() {
            *loader
                .discards
                .entry(Identifier::Gene(gene_id.clone()))
                .or_default() |= Reason::NoTranscripts;
        }
        if let Some(gene) = loader.gene_id_to_gene.get(gene_id) {
            if gene
                .biotype
                .as_ref()
                .is_some_and(|bt| bt.iter().any(|b| DISCARD_BIOTYPES_GENES.contains(b)))
            {
                *loader
                    .discards
                    .entry(Identifier::Gene(gene_id.clone()))
                    .or_default() |= Reason::Biotype;
            }
        } else {
            *loader
                .discards
                .entry(Identifier::Gene(gene_id.clone()))
                .or_default() |= Reason::MissingGene;
        }
    }
    // loader.discard()?;
    tracing::info!("… done filter Gene Id entries");
    Ok(())
}

pub(crate) fn filter_genes(loader: &mut TranscriptLoader) -> Result<(), Error> {
    if loader.disable_filters {
        return Ok(());
    }
    tracing::info!("Filtering genes …");
    // Check whether the gene has a gene symbol.
    let missing_symbol = |_gene_id: GeneId, gene: &Gene, _txs: &[&Transcript]| -> bool {
        gene.gene_symbol.is_none() || gene.gene_symbol.as_ref().unwrap().is_empty()
    };
    // Check whether the gene has a biotype that should be discarded.
    let biotype = |_gene_id: GeneId, gene: &Gene, _txs: &[&Transcript]| -> bool {
        gene.biotype
            .as_ref()
            .map(|bt| bt.iter().any(|bt| DISCARD_BIOTYPES_GENES.contains(bt)))
            .unwrap_or(false)
    };
    // Check whether all transcripts for a gene are predicted.
    let predicted_only = |_gene_id: GeneId, _gene: &Gene, txs: &[&Transcript]| -> bool {
        !txs.is_empty() && txs.iter().all(|tx| tx.id.starts_with('X'))
    };
    type Filter = fn(GeneId, &Gene, &[&Transcript]) -> bool;
    let filters: [(Filter, Reason); 3] = [
        (missing_symbol, Reason::MissingGeneSymbol),
        (biotype, Reason::Biotype),
        (predicted_only, Reason::PredictedTranscriptsOnly),
    ];

    let discarded_genes = loader
        .gene_id_to_gene
        .iter()
        .filter_map(|(gene_id, gene)| {
            let txs = loader
                .gene_id_to_transcript_ids
                .get(gene_id)
                .map(|tx_ids| {
                    tx_ids
                        .iter()
                        .filter_map(|tx_id| loader.transcript_id_to_transcript.get(tx_id))
                        .collect_vec()
                })
                .unwrap_or_default();
            let reason = filters
                .iter()
                .filter_map(|(f, r)| f(gene_id.clone(), gene, &txs).then_some(*r))
                .fold(BitFlags::<Reason>::default(), |a, b| a | b);
            let empty = reason.is_empty();
            (!empty).then_some((Identifier::Gene(gene_id.clone()), reason))
        })
        .collect_vec();

    for (id, reason) in discarded_genes {
        loader.mark_discarded(&id, reason)?;
    }
    // loader.discard()?;
    tracing::info!("… done filtering genes");
    Ok(())
}

/// Filter transcripts for gene.
///
/// We employ the following rules:
///
/// - Remove redundant transcripts with the same identifier and pick only the
///   transcripts that have the highest version number for one assembly.
/// - Do not pick any `XM_`/`XR_` (NCBI predicted only) transcripts.
/// - Do not pick any `NR_` transcripts when there are coding `NM_` transcripts.
pub(crate) fn filter_transcripts(loader: &mut TranscriptLoader) -> Result<(), Error> {
    if loader.disable_filters {
        return Ok(());
    }
    tracing::info!("Filtering transcripts …");

    /// Params used as args for filter functions defined below
    #[derive(new)]
    struct Params<'a> {
        tx: &'a Transcript,
        tx_id: &'a str,
    }

    // Check whether the transcript has an HGNC id.
    let missing_hgnc = |p: &Params| -> bool {
        p.tx.hgnc.is_none()
            || p.tx.hgnc.as_ref().unwrap().is_empty()
            || p.tx.hgnc.as_ref().unwrap().starts_with("GENE:")
    };

    // Check whether the transcript has any genome builds.
    let empty_genome_builds = |p: &Params| -> bool { p.tx.genome_builds.is_empty() };

    // Check whether the transcript has tags cds_start_NF or cds_end_NF.
    let cds_start_end_nf = |p: &Params| -> bool {
        p.tx.genome_builds.values().any(|gb| {
            if let Some(tag) = &gb.tag {
                tag.contains(&Tag::Other("cds_start_NF".into()))
                    || tag.contains(&Tag::Other("cds_end_NF".into()))
            } else {
                false
            }
        })
    };

    // Check whether the transcript is marked as partial.
    // Ignore "partial" tag for NR transcripts.
    // For protein coding transcripts which are marked as partial,
    // check whether there's at least 1bp of 5' and 3' UTR
    // (and then discard the 'partial' tag).
    let partial = |p: &Params| -> bool {
        if !matches!(p.tx.partial, Some(1)) {
            return false;
        }
        if p.tx_id.starts_with("NR") {
            false
        } else if p.tx.protein_coding() {
            // Safely handle build lookup without panicking
            if let Some(alignment) = p.tx.genome_builds.values().next() {
                // Flag partial when either end is truncated, or when there are no exons
                truncated_ends(alignment)
                    .is_none_or(|(five_prime, three_prime)| five_prime || three_prime)
            } else {
                // No builds available, treat as non-partial
                false
            }
        } else {
            true
        }
    };

    // Check whether the transcript is predicted.
    let predicted = |p: &Params| -> bool { p.tx_id.starts_with('X') };

    // Check whether the transcript has a biotype that should be discarded.
    let biotype = |p: &Params| -> bool {
        p.tx.biotype
            .as_ref()
            .map(|bt| {
                bt.iter()
                    .any(|bt| DISCARD_BIOTYPES_TRANSCRIPTS.contains(bt))
            })
            .unwrap_or(false)
    };

    // We have two sets of filters here:
    // `tx_filters`: for every transcript known to us at this point (so every transcript from cdot)
    // `gene_id_grouped_tx_filters`: all transcripts associated to a gene id
    // The latter set relies on the grouping of transcripts by gene id,
    // e.g. to discard transcripts with an older version within the same gene id group
    type Filter = fn(&Params) -> bool;
    let tx_filters: [(Filter, Reason); 6] = [
        (missing_hgnc, Reason::MissingHgncId),
        (empty_genome_builds, Reason::EmptyGenomeBuilds),
        (partial, Reason::OnlyPartialAlignmentInRefSeq),
        (predicted, Reason::PredictedTranscript),
        (biotype, Reason::Biotype),
        (cds_start_end_nf, Reason::CdsStartOrEndNotConfirmed),
    ];

    #[derive(new)]
    struct GroupParams<'a> {
        idx: usize,
        tx_ids: &'a [TranscriptId],
        txs: &'a [&'a Transcript],
    }
    type GroupFilter = fn(&GroupParams) -> bool;

    // Check whether we have a newer version of the same transcript accession available.
    let old_version = |p: &GroupParams| -> bool {
        let versioned = group_transcripts_by_release_and_version(p.txs);
        let (tx_ac, tx_version) = p.tx_ids[p.idx].split_version();
        for release in p.txs[p.idx].genome_builds.keys() {
            if let Some(other) = versioned.get(&(release.into(), tx_ac.to_string()))
                && other.iter().any(|(version, _tx)| *version > tx_version)
            {
                return true;
            }
        }
        false
    };

    // Check whether we have already seen an NM transcript for the gene.
    let nm_instead_of_nr = |p: &GroupParams| -> bool {
        p.tx_ids[p.idx].starts_with("NR_") && p.tx_ids.iter().any(|tx_id| tx_id.starts_with("NM_"))
    };

    let gene_id_grouped_tx_filters: [(GroupFilter, Reason); 2] = [
        (old_version, Reason::OldVersion),
        (nm_instead_of_nr, Reason::UseNmTranscriptInsteadOfNr),
    ];

    let start = Instant::now();

    // Apply first set of filters (which do not depend on gene id grouping)
    let discarded = loader
        .transcript_id_to_transcript
        .iter()
        .filter_map(|(tx_id, tx)| {
            let p = Params::new(tx, tx_id);
            let reason = tx_filters
                .iter()
                .filter_map(|(f, r)| f(&p).then_some(*r))
                .fold(BitFlags::<Reason>::default(), |a, b| a | b);
            let empty = reason.is_empty();
            (!empty).then(|| (Identifier::Transcript(tx_id.clone()), reason))
        })
        .collect_vec();

    for (id, reason) in discarded {
        loader.mark_discarded(&id, reason)?;
    }

    // Apply second set of filters (which depend on gene id grouping)
    let discarded = loader
        .gene_id_to_transcript_ids
        .values()
        .filter_map(|tx_ids| {
            let txs = tx_ids
                .iter()
                .map(|tx_id| loader.transcript_id_to_transcript.get(tx_id).unwrap())
                .collect_vec();
            let tx_reasons = tx_ids
                .iter()
                .enumerate()
                .filter_map(|(idx, tx)| {
                    let p = GroupParams::new(idx, tx_ids, &txs);
                    let reason = gene_id_grouped_tx_filters
                        .iter()
                        .filter_map(|(f, r)| f(&p).then_some(*r))
                        .fold(BitFlags::<Reason>::default(), |a, b| a | b);
                    let empty = reason.is_empty();
                    if !empty {
                        Some((Identifier::Transcript(tx.clone()), reason))
                    } else {
                        None
                    }
                })
                .collect_vec();
            if tx_reasons.is_empty() {
                None
            } else {
                Some(tx_reasons)
            }
        })
        .flatten()
        .collect_vec();

    for (id, reason) in discarded {
        loader.mark_discarded(&id, reason)?;
    }
    // loader.discard()?;

    tracing::info!("… done filtering transcripts in {:?}", start.elapsed());
    Ok(())
}

pub(crate) fn filter_transcripts_with_sequence(
    loader: &mut TranscriptLoader,
    seq_provider: &mut SequenceProvider,
    progress: &dyn Progress,
) -> Result<HashMap<TranscriptId, String>, Error> {
    tracing::info!("Filtering transcripts with sequences …");
    let start = Instant::now();

    let five_prime_truncated = |tx: &Transcript| -> bool {
        let is_mt = MITOCHONDRIAL_ACCESSIONS.iter().any(|a| tx.is_on_contig(a));
        if tx.protein_coding() && !is_mt {
            tx.genome_builds
                .values()
                .any(|alignment| truncated_ends(alignment).is_none_or(|(five_prime, _)| five_prime))
        } else {
            false
        }
    };
    let three_prime_truncated = |tx: &Transcript| -> bool {
        let is_mt = MITOCHONDRIAL_ACCESSIONS.iter().any(|a| tx.is_on_contig(a));
        if tx.protein_coding() && !is_mt {
            tx.genome_builds.values().any(|alignment| {
                truncated_ends(alignment).is_none_or(|(_, three_prime)| three_prime)
            })
        } else {
            false
        }
    };

    // Check transcript's CDS length for being a multiple of 3.
    //
    // Note that the chrMT transcripts have been fixed earlier already to
    // accommodate for how they are fixed by poly-A tailing.
    let invalid_cds_length = |tx: &Transcript| -> bool {
        if tx.protein_coding() {
            tx.cds_length().is_none_or(|l| l % 3 != 0)
        } else {
            false
        }
    };

    let bar = progress.bar(
        "Fetching transcript sequences",
        loader.transcript_id_to_transcript.len() as u64,
        Unit::Transcripts,
    );
    let (discards, keeps_with_reasons): (Vec<_>, Vec<_>) = loader
        .transcript_id_to_transcript
        .par_iter()
        .progress_with(bar.clone())
        .partition_map(|(tx_id, tx)| {
            if let Some(d) = loader.discards.get(&Identifier::Transcript(tx_id.clone()))
                && d.intersects(Reason::hard())
            {
                return Either::Left((Identifier::Transcript(tx_id.clone()), *d));
            }

            let has_invalid_cds_length = if loader.disable_filters {
                false
            } else {
                invalid_cds_length(tx)
            };

            let namespace: String = if tx_id.starts_with("ENST") {
                String::from("Ensembl")
            } else {
                String::from("NCBI")
            };
            let seq = seq_provider.fetch_sequence(&AliasOrSeqId::Alias {
                value: (*tx_id).to_string(),
                namespace: Some(namespace.clone()),
            });

            if let Ok(seq) = seq {
                let is_mt = MITOCHONDRIAL_ACCESSIONS.iter().any(|a| tx.is_on_contig(a));
                let is_seleno = tx
                    .biotype
                    .as_ref()
                    .map(|bt| bt.contains(&BioType::Selenoprotein))
                    .unwrap_or(false);
                if seq.is_empty() {
                    return Either::Left((
                        Identifier::Transcript(tx_id.clone()),
                        Reason::MissingSequence.into(),
                    ));
                };
                // Skip transcript if it is coding and the translated CDS does not have a stop codon.
                let cds = tx
                    .start_codon
                    .and_then(|start| tx.stop_codon.map(|stop| (start as usize, stop as usize)));
                let cds = if tx.protein_coding() { cds } else { None };
                if let Some((cds_start, cds_end)) = cds {
                    // `fix_cds` completes a stop codon at the transcript end with 1 or 2
                    // bases after the last exon. If the sequence ends before the exons or
                    // the CDS, pad it with `A` bases.
                    let tx_end = tx
                        .genome_builds
                        .values()
                        .flat_map(|gb| gb.exons.iter())
                        .filter_map(|exon| usize::try_from(exon.alt_cds_end_i).ok())
                        .max()
                        .unwrap_or_default();
                    let delta = tx_end.max(cds_end).saturating_sub(seq.len());
                    let seq = append_poly_a(seq, delta);

                    let safe_end = cds_end.min(seq.len());
                    let safe_start = cds_start.min(safe_end);
                    // An incomplete CDS end can leave a partial last codon. It is no stop codon.
                    let safe_end = safe_end - (safe_end - safe_start) % 3;

                    let tx_seq_to_translate = &seq[safe_start..safe_end];

                    let aa_sequence_result = translate_cds(
                        tx_seq_to_translate,
                        true,
                        "*",
                        if is_mt {
                            TranslationTable::VertebrateMitochondrial
                        } else if is_seleno {
                            TranslationTable::Selenocysteine
                        } else {
                            TranslationTable::Standard
                        },
                    );

                    let mut reason = Reason::empty();

                    if let Ok(aa_sequence) = aa_sequence_result {
                        let has_missing_stop_codon = (!is_mt && !aa_sequence.ends_with('*'))
                            || (is_mt && !aa_sequence.contains('*'));
                        // Annotation reads this flag, so set it also without filters.
                        if has_missing_stop_codon {
                            reason |= Reason::MissingStopCodon;
                        }
                        if has_missing_stop_codon && !loader.disable_filters {
                            if five_prime_truncated(tx) {
                                reason |= Reason::FivePrimeEndTruncated;
                            }
                            if three_prime_truncated(tx) {
                                reason |= Reason::ThreePrimeEndTruncated;
                            }
                        }
                    }

                    if has_invalid_cds_length
                        || loader
                            .fixes
                            .get(&Identifier::Transcript(tx_id.clone()))
                            .is_some_and(|f| f.contains(Fix::Cds))
                    {
                        reason |= Reason::InvalidCdsLength;
                    }

                    if reason.intersects(Reason::hard()) {
                        Either::Left((Identifier::Transcript(tx_id.clone()), reason))
                    } else {
                        Either::Right((tx_id.clone(), seq, reason))
                    }
                } else {
                    Either::Right((tx_id.clone(), seq, Reason::empty()))
                }
            } else {
                Either::Left((
                    Identifier::Transcript(tx_id.clone()),
                    Reason::MissingSequence.into(),
                ))
            }
        });
    bar.finish();

    for (id, reason) in discards {
        loader.mark_discarded(&id, reason)?;
    }
    // loader.discard()?;

    let mut sequence_map = HashMap::new();
    for (tx_id, seq, reason) in keeps_with_reasons {
        if !reason.is_empty() {
            loader.mark_discarded(&Identifier::Transcript(tx_id.clone()), reason)?;
        }
        sequence_map.insert(tx_id, seq);
    }

    tracing::info!(
        "… done filtering transcripts with sequence provider in {:?}",
        start.elapsed()
    );

    Ok(sequence_map)
}

/// Remove all gene ids for which there are no transcripts left.
pub(crate) fn filter_empty_gene_id_mappings(loader: &mut TranscriptLoader) -> Result<(), Error> {
    tracing::info!("Removing empty Gene Id mappings …");
    let gene_ids_gene: HashSet<_> = loader.gene_id_to_gene.keys().collect();
    let gene_ids_tx: HashSet<_> = loader.gene_id_to_transcript_ids.keys().collect();
    let gene_id_but_no_tx_id = &gene_ids_gene - &gene_ids_tx;
    let tx_id_but_no_gene_id = &gene_ids_tx - &gene_ids_gene;
    for gene_id in gene_id_but_no_tx_id {
        *loader
            .discards
            .entry(Identifier::Gene(gene_id.clone()))
            .or_default() |= Reason::NoTranscripts;
    }
    for gene_id in tx_id_but_no_gene_id {
        *loader
            .discards
            .entry(Identifier::Gene(gene_id.clone()))
            .or_default() |= Reason::MissingGene;
    }

    for (gene_id, txs) in loader.gene_id_to_transcript_ids.iter() {
        if txs.is_empty() {
            *loader
                .discards
                .entry(Identifier::Gene(gene_id.clone()))
                .or_default() |= Reason::NoTranscriptLeft;
        }
    }
    // loader.discard()?;
    tracing::info!("… done removing empty HGNC mappings");
    Ok(())
}

/// Mitochondrial accessions.
pub const MITOCHONDRIAL_ACCESSIONS: &[&str] = &["NC_012920.1", "NC_001807.4"];
static DISCARD_BIOTYPES_TRANSCRIPTS: Lazy<HashSet<BioType>> = Lazy::new(|| {
    HashSet::from([
        // BioType::PseudogenicTranscript,
        BioType::IgCGene,
        BioType::IgDGene,
        BioType::IgJGene,
        BioType::IgVGene,
        BioType::TrCGene,
        BioType::TrDGene,
        BioType::TrJGene,
        BioType::TrVGene,
        BioType::AberrantProcessedTranscript,
        BioType::UnconfirmedTranscript,
        BioType::NmdTranscriptVariant,
    ])
});
// this used to contain `BioType::Pseudogene`,
// but we keep pseudogenes (as there are cases where a protein still is produced in certain individuals)
static DISCARD_BIOTYPES_GENES: Lazy<HashSet<BioType>> = Lazy::new(|| {
    HashSet::from([
        BioType::CGeneSegment,
        BioType::DGeneSegment,
        BioType::JGeneSegment,
        BioType::VGeneSegment,
    ])
});

/// Group transcripts by release and version (sorted descending).
fn group_transcripts_by_release_and_version<'a>(
    txs: &[&'a Transcript],
) -> HashMap<(String, String), Vec<(u32, &'a Transcript)>> {
    let mut versioned = txs
        .iter()
        .flat_map(|tx| {
            let releases = tx.genome_builds.keys().cloned().collect_vec();
            let tx_id = TranscriptId::try_new(&tx.id).unwrap();
            let (ac, version) = tx_id.split_version();
            let ac = ac.to_string();
            releases
                .into_iter()
                .map(move |r| ((r, ac.clone()), (version, *tx)))
        })
        .into_group_map();
    versioned.values_mut().for_each(|txs| {
        txs.sort_unstable_by_key(|(version, _)| Reverse(*version));
    });
    versioned
}

/// Returns `(5' truncated, 3' truncated)`: whether the CDS reaches the 5' end and the 3' end
/// of the transcript, i.e. whether the 5' UTR and the 3' UTR are missing.
/// Returns `None` if the alignment has no exons.
///
/// `cds_start`, `cds_end` and the exon coordinates are genomic on both strands.
/// The 5' end of a minus-strand transcript is therefore its highest coordinate.
fn truncated_ends(alignment: &GenomeAlignment) -> Option<(bool, bool)> {
    let tx_start = alignment.exons.iter().map(|e| e.alt_start_i).min()?;
    let tx_end = alignment.exons.iter().map(|e| e.alt_end_i).max()?;
    let low_end = alignment.cds_start == Some(tx_start);
    let high_end = alignment.cds_end == Some(tx_end);
    Some(match alignment.strand {
        Strand::Plus => (low_end, high_end),
        Strand::Minus => (high_end, low_end),
    })
}

fn append_poly_a(seq: String, length: usize) -> String {
    let mut seq = seq.into_bytes();
    seq.extend_from_slice(b"A".repeat(length).as_slice());
    String::from_utf8(seq).expect("must be valid UTF-8")
}

/// For each transcript that has been discarded for whatever reason, propagate the reason to its
/// parent gene id entry *if* the gene id entry already exists and has a non-empty reason
/// (to avoid erroneously discarding gene id entries for a non-important reason).
pub(crate) fn propagate_discard_reasons(loader: &mut TranscriptLoader) -> Result<(), Error> {
    // First check whether all transcripts of a gene have been marked as discarded.
    for (gene_id, _) in loader.gene_id_to_gene.iter() {
        let tx_ids = loader
            .gene_id_to_transcript_ids
            .get(gene_id)
            .map(|v| v.as_slice())
            .unwrap_or_default();
        if !tx_ids.is_empty()
            && tx_ids.iter().all(|tx_id| {
                loader
                    .discards
                    .get(&Identifier::Transcript(tx_id.clone()))
                    .is_some_and(|d| d.intersects(Reason::hard()))
            })
        {
            *loader
                .discards
                .entry(Identifier::Gene(gene_id.clone()))
                .or_default() |= Reason::NoTranscriptLeft;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::progress::NoProgress;
    use hgvs::data::cdot::json::models::Exon;
    use indexmap::IndexMap;

    const TX_ID: &str = "NM_000001.1";

    /// Loader with one partial, coding transcript whose exons are at genomic
    /// `[100, 200)` and `[300, 400)`. Its all-`C` sequence has no stop codon, so
    /// `filter_transcripts_with_sequence` runs its 5'/3' truncation checks.
    fn loader_with_tx(strand: Strand, cds_start: i32, cds_end: i32) -> TranscriptLoader {
        let (low_ord, high_ord) = match strand {
            Strand::Plus => (0, 1),
            Strand::Minus => (1, 0),
        };
        let exon = |ord: i32, alt_start_i: i32| Exon {
            alt_start_i,
            alt_end_i: alt_start_i + 100,
            ord,
            alt_cds_start_i: ord * 100 + 1,
            alt_cds_end_i: ord * 100 + 100,
            cigar: "100M".to_string(),
        };
        let alignment = GenomeAlignment {
            cds_start: Some(cds_start),
            cds_end: Some(cds_end),
            contig: "NC_000017.11".to_string(),
            exons: vec![exon(low_ord, 100), exon(high_ord, 300)],
            strand,
            tag: None,
            note: None,
        };
        let tx = Transcript {
            biotype: None,
            gene_name: None,
            gene_version: String::new(),
            genome_builds: IndexMap::from([("GRCh38".to_string(), alignment)]),
            hgnc: None,
            id: TX_ID.to_string(),
            partial: Some(1),
            protein: Some("NP_000001.1".to_string()),
            start_codon: Some(0),
            stop_codon: Some(99),
            transl_except: None,
            transl_table: None,
        };
        let mut loader = TranscriptLoader::new("GRCh38".to_string(), false);
        loader
            .transcript_id_to_transcript
            .insert(TranscriptId::try_new(TX_ID).unwrap(), tx);
        loader
    }

    /// `cds_start` and `cds_end` are genomic on both strands. The 5' end of a
    /// minus-strand transcript is therefore its highest coordinate.
    #[rstest::rstest]
    #[case::plus_no_5p_utr(Strand::Plus, 100, 350, true, false)]
    #[case::plus_no_3p_utr(Strand::Plus, 150, 400, false, true)]
    #[case::plus_both_utrs(Strand::Plus, 150, 350, false, false)]
    #[case::minus_no_5p_utr(Strand::Minus, 150, 400, true, false)]
    #[case::minus_no_3p_utr(Strand::Minus, 100, 350, false, true)]
    #[case::minus_both_utrs(Strand::Minus, 150, 350, false, false)]
    fn truncation_is_checked_at_strand_aware_ends(
        #[case] strand: Strand,
        #[case] cds_start: i32,
        #[case] cds_end: i32,
        #[case] five_prime: bool,
        #[case] three_prime: bool,
    ) -> Result<(), Error> {
        let mut loader = loader_with_tx(strand, cds_start, cds_end);
        let mut seq_provider =
            SequenceProvider::FastaMap(HashMap::from([(TX_ID.to_string(), "C".repeat(200))]));

        filter_transcripts(&mut loader)?;
        filter_transcripts_with_sequence(&mut loader, &mut seq_provider, &NoProgress)?;

        let reason = loader.discards[&Identifier::Transcript(TranscriptId::try_new(TX_ID)?)];
        assert_eq!(
            (
                reason.contains(Reason::OnlyPartialAlignmentInRefSeq),
                reason.contains(Reason::FivePrimeEndTruncated),
                reason.contains(Reason::ThreePrimeEndTruncated),
            ),
            (five_prime || three_prime, five_prime, three_prime),
        );
        Ok(())
    }

    /// A coding transcript with one exon at genomic `[1000, 1000 + len)` and the CDS
    /// `0..stop_codon`.
    fn one_exon_tx(contig: &str, len: usize, stop_codon: i32) -> Result<Transcript, Error> {
        let len = i32::try_from(len)?;
        let alignment = GenomeAlignment {
            cds_start: Some(1000),
            cds_end: Some(1000 + stop_codon),
            contig: contig.to_string(),
            exons: vec![Exon {
                alt_start_i: 1000,
                alt_end_i: 1000 + len,
                ord: 0,
                alt_cds_start_i: 1,
                alt_cds_end_i: len,
                cigar: format!("{len}M"),
            }],
            strand: Strand::Plus,
            tag: None,
            note: None,
        };
        Ok(Transcript {
            biotype: None,
            gene_name: None,
            gene_version: String::new(),
            genome_builds: IndexMap::from([("GRCh38".to_string(), alignment)]),
            hgnc: None,
            id: TX_ID.to_string(),
            partial: None,
            protein: Some("NP_000001.1".to_string()),
            start_codon: Some(0),
            stop_codon: Some(stop_codon),
            transl_except: None,
            transl_table: None,
        })
    }

    /// Apply `fix_cds` and `filter_transcripts_with_sequence` to `tx` with the sequence
    /// `seq`. Return the stored sequence, the discard reasons and the stop codon.
    fn stored_sequence(
        tx: Transcript,
        seq: &str,
    ) -> Result<(String, BitFlags<Reason>, Option<i32>), Error> {
        stored_sequence_with(tx, seq, false)
    }

    /// `stored_sequence` for a loader with the given `disable_filters`.
    fn stored_sequence_with(
        tx: Transcript,
        seq: &str,
        disable_filters: bool,
    ) -> Result<(String, BitFlags<Reason>, Option<i32>), Error> {
        let tx_id = TranscriptId::try_new(TX_ID)?;
        let mut loader = TranscriptLoader::new("GRCh38".to_string(), disable_filters);
        loader.transcript_id_to_transcript.insert(tx_id.clone(), tx);
        loader.fix_cds();

        let mut seq_provider =
            SequenceProvider::FastaMap(HashMap::from([(TX_ID.to_string(), seq.to_string())]));
        let mut seqs =
            filter_transcripts_with_sequence(&mut loader, &mut seq_provider, &NoProgress)?;
        let reason = loader
            .discards
            .get(&Identifier::Transcript(tx_id.clone()))
            .copied()
            .unwrap_or_default();
        let stop_codon = loader.transcript_id_to_transcript[&tx_id].stop_codon;
        Ok((seqs.remove(&tx_id).unwrap_or_default(), reason, stop_codon))
    }

    #[rstest::rstest]
    #[case::complete_cds("NC_000001.11", "ATGAAATAGCC", 9, "ATGAAATAGCC")]
    #[case::stop_codon_at_transcript_end("NC_000001.11", "ATGAAATAG", 9, "ATGAAATAG")]
    fn complete_cds_is_not_padded(
        #[case] contig: &str,
        #[case] seq: &str,
        #[case] stop_codon: i32,
        #[case] expected: &str,
    ) -> Result<(), Error> {
        let tx = one_exon_tx(contig, seq.len(), stop_codon)?;
        assert_eq!(stored_sequence(tx, seq)?.0, expected);
        Ok(())
    }

    /// Polyadenylation completes a stop codon such as `T` or `TA` to `TAA`, on chrMT and
    /// elsewhere. It can do so only at the transcript end.
    #[rstest::rstest]
    #[case::chr_mt("NC_012920.1", "ATGAAAT", 7, "ATGAAATAA", 9)]
    #[case::nuclear("NC_000002.12", "ATGAAATA", 8, "ATGAAATAA", 9)]
    #[case::before_transcript_end("NC_000001.11", "ATGAAATAGCC", 8, "ATGAAATAGCC", 8)]
    fn stop_codon_is_completed_only_at_the_transcript_end(
        #[case] contig: &str,
        #[case] seq: &str,
        #[case] stop_codon: i32,
        #[case] expected_seq: &str,
        #[case] expected_stop_codon: i32,
    ) -> Result<(), Error> {
        let tx = one_exon_tx(contig, seq.len(), stop_codon)?;
        let (stored, reason, stop_codon) = stored_sequence(tx, seq)?;
        assert_eq!(stored, expected_seq);
        assert_eq!(stop_codon, Some(expected_stop_codon));
        // Without a completed stop codon, the CDS has none.
        assert_eq!(
            reason.contains(Reason::MissingStopCodon),
            expected_stop_codon % 3 != 0,
            "{reason:?}"
        );
        Ok(())
    }

    /// A CDS whose end the annotation marks as incomplete keeps its partial last codon. It has
    /// no stop codon.
    #[rstest::rstest]
    #[case::cds_end_nf(Some(Tag::Other("cds_end_NF".into())), None)]
    #[case::partial(None, Some(1))]
    fn incomplete_cds_end_is_not_padded(
        #[case] tag: Option<Tag>,
        #[case] partial: Option<u8>,
    ) -> Result<(), Error> {
        let mut tx = one_exon_tx("NC_000001.11", 7, 7)?;
        for alignment in tx.genome_builds.values_mut() {
            alignment.tag = tag.clone().map(|tag| vec![tag]);
        }
        tx.partial = partial;

        let (seq, reason, stop_codon) = stored_sequence(tx, "ATGAAAT")?;
        assert_eq!(seq, "ATGAAAT");
        assert_eq!(stop_codon, Some(7));
        assert!(reason.contains(Reason::MissingStopCodon), "{reason:?}");
        Ok(())
    }

    /// Annotation reads `MissingStopCodon` as a fact about the CDS, so `db create` sets it also
    /// without filters.
    #[rstest::rstest]
    #[case::with_filters(false)]
    #[case::without_filters(true)]
    fn missing_stop_codon_is_flagged(#[case] disable_filters: bool) -> Result<(), Error> {
        let tx = one_exon_tx("NC_000001.11", 9, 9)?;
        let (_, reason, _) = stored_sequence_with(tx, "ATGAAAAAA", disable_filters)?;
        assert!(reason.contains(Reason::MissingStopCodon), "{reason:?}");
        Ok(())
    }
}
