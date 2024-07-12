//! Transcript database.

use std::cmp::{PartialEq, Reverse};
use std::collections::{HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::{io::Write, path::PathBuf, time::Instant};

use anyhow::{anyhow, Error};
use clap::Parser;
use derive_new::new;
use enumflags2::{bitflags, BitFlag, BitFlags};
use hgvs::data::cdot::json::models;
use hgvs::data::cdot::json::models::{BioType, Gene, Tag, Transcript};
use hgvs::sequences::{translate_cds, TranslationTable};
use itertools::Itertools;
use nom::ToUsize;
use nutype::nutype;
use once_cell::sync::Lazy;
use prost::Message;
use rayon::iter::Either;
use rayon::prelude::*;
use seqrepo::{AliasOrSeqId, Interface, SeqRepo};
use serde::Serialize;
use serde_json::json;
use serde_with::serde_as;
use serde_with::DisplayFromStr;
use strum::Display;
use thousands::Separable;

use crate::common::{trace_rss_now, GenomeRelease};
use crate::pbs::txs::TxSeqDatabase;

/// Mitochondrial accessions.
const MITOCHONDRIAL_ACCESSIONS: &[&str] = &["NC_012920.1"];
const MITOCHONDRIAL_ACCESSION: &str = "NC_012920.1";

/// Command line arguments for `db create txs` sub command.
#[derive(Parser, Debug)]
#[command(about = "Construct mehari transcripts and sequence database", long_about = None)]
pub struct Args {
    /// Genome release to extract transcripts for.
    #[arg(long)]
    pub genome_release: GenomeRelease,

    /// Path to output protobuf file to write to.
    #[arg(long)]
    pub path_out: PathBuf,

    /// Paths to the cdot JSON transcripts to import.
    #[arg(long, required = true)]
    pub path_cdot_json: Vec<PathBuf>,

    /// Path to the seqrepo instance directory to use.
    #[arg(long)]
    pub path_seqrepo_instance: PathBuf,

    /// Path to TSV file for label transfer of transcripts.  Columns are
    /// transcript id (without version), (unused) gene symbol, and label.
    #[arg(long)]
    pub path_mane_txs_tsv: Option<PathBuf>,

    /// Maximal number of transcripts to process. DEPRECATED.
    #[arg(long)]
    pub max_txs: Option<u32>,

    /// Limit transcript database to the following HGNC symbols.  Useful for
    /// building test databases.
    #[arg(long)]
    pub gene_symbols: Option<Vec<String>>,

    /// Number of threads to use for steps supporting parallel processing.
    #[arg(long, default_value = "1")]
    pub threads: usize,
}

/// Helper struct for parsing the label TSV file.
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize)]
struct LabelEntry {
    /// Transcript identifier without version.
    transcript_id: String,
    /// Transcript version.
    transcript_version: usize,
    /// Gene symbol (unused).
    _gene_symbol: String,
    /// Label to transfer.
    label: String,
}
#[nutype(derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, FromStr, Borrow, Deref
))]
struct HgncId(usize);

impl Display for HgncId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "HGNC:{}", self.to_usize())
    }
}

#[nutype(
    sanitize(trim),
    validate(not_empty),
    derive(
        Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, AsRef, Deref, Borrow, Into,
        Display
    )
)]
struct TranscriptId(String);

impl TranscriptId {
    fn without_version(&self) -> Result<Self, Error> {
        Ok(Self::try_new(
            self.as_ref()
                .split('.')
                .next()
                .ok_or(anyhow!("No version to split off"))?,
        )?)
    }

    fn split_version(&self) -> (&str, u32) {
        let (ac, version) = self
            .rsplit_once('.')
            .expect("Invalid accession, expected format 'ac.version'.");
        (ac, version.parse::<u32>().expect("invalid version"))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Hash)]
#[serde(tag = "type", content = "value")]
enum Identifier {
    Hgnc(HgncId),
    TxId(TranscriptId),
}
trait TranscriptExt {
    fn cds_length(&self) -> Option<u32>;

    fn protein_coding(&self) -> bool;

    fn is_on_contig(&self, contig: &str) -> bool;
}
impl TranscriptExt for Transcript {
    fn cds_length(&self) -> Option<u32> {
        self.start_codon
            .and_then(|start| self.stop_codon.map(|stop| stop.abs_diff(start)))
    }

    fn protein_coding(&self) -> bool {
        matches!(self.protein.as_ref().map(|p| !p.is_empty()), Some(true))
            || self
                .biotype
                .as_ref()
                .map_or(false, |bt| bt.iter().any(|b| b.is_protein_coding()))
    }

    fn is_on_contig(&self, contig: &str) -> bool {
        self.genome_builds.values().any(|gb| &gb.contig == contig)
    }
}

trait BioTypeExt {
    fn is_protein_coding(&self) -> bool;
}

impl BioTypeExt for BioType {
    fn is_protein_coding(&self) -> bool {
        // see https://www.ensembl.org/info/genome/genebuild/biotypes.html
        // and https://www.ensembl.org/Help/Faq?id=468
        matches!(
            self,
            BioType::ProteinCoding
                | BioType::MRna
                | BioType::IgCGene
                | BioType::IgDGene
                | BioType::IgVGene
                | BioType::TrCGene
                | BioType::TrDGene
                | BioType::TrJGene
                | BioType::TrVGene
        )
    }
}

static DISCARD_BIOTYPES_TRANSCRIPTS: Lazy<HashSet<BioType>> = Lazy::new(|| {
    HashSet::from([
        // BioType::PseudogenicTranscript,
        BioType::AberrantProcessedTranscript,
        BioType::UnconfirmedTranscript,
        BioType::NmdTranscriptVariant,
    ])
});

// this used to contain `BioType::Pseudogene`,
// but we keep pseudogenes (as there are cases where a protein still is produced in certain individuals)
static DISCARD_BIOTYPES_GENES: Lazy<HashSet<BioType>> = Lazy::new(|| HashSet::from([]));

#[derive(Debug, Clone, Default)]
struct TranscriptLoader {
    genome_release: GenomeRelease,
    cdot_version: String,
    transcript_id_to_transcript: HashMap<TranscriptId, Transcript>,
    hgnc_id_to_gene: HashMap<HgncId, Gene>,
    hgnc_id_to_transcript_ids: HashMap<HgncId, Vec<TranscriptId>>,
    discards: HashMap<Identifier, BitFlags<Reason>>,
    fixes: HashMap<Identifier, BitFlags<Fix>>,
}

impl TranscriptLoader {
    fn new(genome_release: GenomeRelease) -> Self {
        Self {
            genome_release,
            ..Default::default()
        }
    }

    fn merge(&mut self, other: &mut Self) -> &mut Self {
        assert_eq!(self.genome_release, other.genome_release);
        assert_eq!(self.cdot_version, other.cdot_version);
        for (hgnc_id, gene) in other.hgnc_id_to_gene.drain() {
            if let Some(old_gene) = self.hgnc_id_to_gene.insert(hgnc_id, gene.clone()) {
                tracing::warn!(
                    "Overwriting gene: {}\n{:?},\n{:?}",
                    &hgnc_id,
                    &old_gene,
                    &gene
                );
            }
        }

        for (tx_id, tx) in other.transcript_id_to_transcript.drain() {
            if let Some(old_tx) = self.transcript_id_to_transcript.insert(tx_id, tx.clone()) {
                tracing::warn!(
                    "Overwriting transcript: {}\n{:?},\n{:?}",
                    old_tx.id,
                    &old_tx,
                    &tx
                );
            }
        }

        for (hgnc_id, mut tx_ids) in other.hgnc_id_to_transcript_ids.drain() {
            let ids = self.hgnc_id_to_transcript_ids.entry(hgnc_id).or_default();
            ids.append(&mut tx_ids);
            ids.sort_unstable();
            ids.dedup();
        }

        for (id, reason) in other.discards.drain() {
            *self.discards.entry(id).or_default() |= reason;
        }
        for (id, fix) in other.fixes.drain() {
            *self.fixes.entry(id).or_default() |= fix;
        }
        self
    }

    /// Load and extract from cdot JSON.
    fn load_cdot(&mut self, path: impl AsRef<Path>) -> Result<(), Error> {
        let models::Container {
            genes: cdot_genes,
            transcripts: cdot_transcripts,
            cdot_version,
            ..
        } = read_cdot_json(path.as_ref())?;
        let cdot_genes = cdot_genes.into_iter().collect::<HashMap<_, _>>();
        let cdot_transcripts = cdot_transcripts
            .into_iter()
            .map(|(txid, tx)| TranscriptId::try_new(txid).map(|t| (t, tx)))
            .collect::<Result<HashMap<_, _>, _>>()?;
        self.cdot_version = cdot_version;

        for (_, gene) in cdot_genes {
            // We're only considering the entries that have an HGNC ID.
            if let Some(hgnc_id) = gene.hgnc.as_ref() {
                let hgnc_id = hgnc_id.parse()?;
                self.hgnc_id_to_gene.insert(hgnc_id, gene);
                self.hgnc_id_to_transcript_ids.entry(hgnc_id).or_default();
            }
        }

        for (tx_id, tx) in cdot_transcripts {
            if let Some(hgnc_id) = tx.hgnc.as_ref() {
                let hgnc_id = hgnc_id.parse()?;
                self.hgnc_id_to_transcript_ids
                    .entry(hgnc_id)
                    .or_default()
                    .push(tx_id.clone());
            }
            self.transcript_id_to_transcript.insert(tx_id, tx);
        }

        Ok(())
    }

    fn filter_initial_hgnc_entries(&mut self) -> Result<(), Error> {
        tracing::info!("Filtering Hgnc entries …");
        for (hgnc_id, tx_ids) in &self.hgnc_id_to_transcript_ids {
            if tx_ids.is_empty() {
                *self.discards.entry(Identifier::Hgnc(*hgnc_id)).or_default() |=
                    Reason::NoTranscripts;
            }
            if let Some(gene) = self.hgnc_id_to_gene.get(hgnc_id) {
                if gene.biotype.as_ref().map_or(false, |bt| {
                    bt.iter().any(|b| DISCARD_BIOTYPES_GENES.contains(b))
                }) {
                    *self.discards.entry(Identifier::Hgnc(*hgnc_id)).or_default() |=
                        Reason::Biotype;
                }
            } else {
                *self.discards.entry(Identifier::Hgnc(*hgnc_id)).or_default() |=
                    Reason::MissingGene;
            }
        }
        self.discard()?;
        tracing::info!("… done filter Hgnc entries");
        Ok(())
    }

    fn apply_fixes(&mut self, transcript_id_to_tags: &Option<HashMap<TranscriptId, Vec<Tag>>>) {
        self.fix_transcript_genome_builds();
        if let Some(txid_to_label) = transcript_id_to_tags {
            self.update_transcript_tags(txid_to_label);
        }
        self.fix_cds();
    }

    fn filter_genes(&mut self) -> Result<(), Error> {
        tracing::info!("Filtering genes …");
        // Check whether the gene has a gene symbol.
        let missing_symbol = |_hgnc_id: HgncId, gene: &Gene, _txs: &[&Transcript]| -> bool {
            gene.gene_symbol.is_none() || gene.gene_symbol.as_ref().unwrap().is_empty()
        };
        // Check whether the gene has a biotype that should be discarded.
        let biotype = |_hgnc_id: HgncId, gene: &Gene, _txs: &[&Transcript]| -> bool {
            gene.biotype
                .as_ref()
                .map(|bt| bt.iter().any(|bt| DISCARD_BIOTYPES_GENES.contains(bt)))
                .unwrap_or(false)
        };
        // Check whether all transcripts for a gene are predicted.
        let predicted_only = |_hgnc_id: HgncId, _gene: &Gene, txs: &[&Transcript]| -> bool {
            txs.iter().all(|tx| tx.id.starts_with('X'))
        };
        type Filter = fn(HgncId, &Gene, &[&Transcript]) -> bool;
        let filters: [(Filter, Reason); 3] = [
            (missing_symbol, Reason::MissingGeneSymbol),
            (biotype, Reason::Biotype),
            (predicted_only, Reason::PredictedTranscriptsOnly),
        ];

        let discarded_genes = self
            .hgnc_id_to_gene
            .iter()
            .filter_map(|(hgnc_id, gene)| {
                let txs = self.hgnc_id_to_transcript_ids[hgnc_id]
                    .iter()
                    .map(|tx_id| &self.transcript_id_to_transcript[tx_id])
                    .collect_vec();
                let reason = filters
                    .iter()
                    .filter_map(|(f, r)| f(*hgnc_id, gene, &txs).then_some(*r))
                    .fold(BitFlags::<Reason>::default(), |a, b| a | b);
                let empty = reason.is_empty();
                (!empty).then_some((Identifier::Hgnc(*hgnc_id), reason))
            })
            .collect_vec();

        for (id, reason) in discarded_genes {
            self.mark_discarded(&id, reason)?;
        }
        self.discard()?;
        tracing::info!("… done filtering genes");
        Ok(())
    }

    fn filter_selected(&mut self, selected_ids: &Vec<Identifier>) -> Result<(), Error> {
        // Discard everything which is not associated with the ids contained in `selected_ids`.
        for id in selected_ids {
            self.mark_discarded(id, Reason::DeselectedGene.into())?;
            match id {
                Identifier::Hgnc(hgnc_id) => {
                    for tx_id in self.hgnc_id_to_transcript_ids.get(hgnc_id).unwrap() {
                        *self
                            .discards
                            .entry(Identifier::TxId(tx_id.clone()))
                            .or_default() |= Reason::DeselectedGene;
                    }
                }
                Identifier::TxId(_tx_id) => {
                    *self.discards.entry(id.clone()).or_default() |= Reason::DeselectedGene;
                }
            }
        }
        self.discard()?;
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
    fn filter_transcripts(&mut self) -> Result<(), Error> {
        tracing::info!("Filtering transcripts …");

        /// Params used as args for filter functions defined below
        #[derive(new)]
        struct Params<'a> {
            tx: &'a Transcript,
            tx_id: &'a str,
        }

        // Check whether the transcript has an HGNC id.
        let missing_hgnc =
            |p: &Params| -> bool { p.tx.hgnc.is_none() || p.tx.hgnc.as_ref().unwrap().is_empty() };

        // Check whether the transcript has any genome builds.
        let empty_genome_builds = |p: &Params| -> bool { p.tx.genome_builds.is_empty() };

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
                assert_eq!(p.tx.genome_builds.len(), 1);
                let alignment = p.tx.genome_builds.values().next().unwrap();
                let exons = alignment.exons.clone();
                if exons.is_empty() {
                    return true;
                }
                let five_prime_trunc = {
                    let cds_start = alignment.cds_start;
                    let first = exons.iter().min_by_key(|e| e.alt_start_i).unwrap();
                    Some(first.alt_start_i) == cds_start
                };
                let three_prime_trunc = {
                    let cds_end = alignment.cds_end;
                    let last = exons.iter().max_by_key(|e| e.alt_end_i).unwrap();
                    Some(last.alt_end_i) == cds_end
                };
                five_prime_trunc && three_prime_trunc
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
        // `hgnc_grouped_tx_filters`: all transcripts associated to an hgnc id
        // The latter set relies on the grouping of transcripts by hgnc id,
        // e.g. to discard transcripts with an older version within the same hgnc group
        type Filter = fn(&Params) -> bool;
        let tx_filters: [(Filter, Reason); 5] = [
            (missing_hgnc, Reason::MissingHgncId),
            (empty_genome_builds, Reason::EmptyGenomeBuilds),
            (partial, Reason::OnlyPartialAlignmentInRefSeq),
            (predicted, Reason::PredictedTranscript),
            (biotype, Reason::Biotype),
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
            let versioned = Self::by_release_and_version(p.txs);
            let (tx_ac, tx_version) = p.tx_ids[p.idx].split_version();
            for release in p.txs[p.idx].genome_builds.keys() {
                if let Some(other) = versioned.get(&(release.into(), tx_ac.to_string())) {
                    if other.iter().any(|(version, _tx)| *version > tx_version) {
                        return true;
                    }
                }
            }
            false
        };

        // Check whether we have already seen an NM transcript for the gene.
        let nm_instead_of_nr = |p: &GroupParams| -> bool {
            p.tx_ids[p.idx].starts_with("NR_")
                && p.tx_ids.iter().any(|tx_id| tx_id.starts_with("NM_"))
        };

        let hgnc_grouped_tx_filters: [(GroupFilter, Reason); 2] = [
            (old_version, Reason::OldVersion),
            (nm_instead_of_nr, Reason::UseNmTranscriptInsteadOfNr),
        ];

        let start = Instant::now();

        // Apply first set of filters (which do not depend on hgnc grouping)
        let discarded = self
            .transcript_id_to_transcript
            .iter()
            .filter_map(|(tx_id, tx)| {
                let p = Params::new(tx, tx_id);
                let reason = tx_filters
                    .iter()
                    .filter_map(|(f, r)| f(&p).then_some(*r))
                    .fold(BitFlags::<Reason>::default(), |a, b| a | b);
                let empty = reason.is_empty();
                (!empty).then(|| (Identifier::TxId(tx_id.clone()), reason))
            })
            .collect_vec();

        for (id, reason) in discarded {
            self.mark_discarded(&id, reason)?;
        }

        // Apply second set of filters (which depend on hgnc grouping)
        let discarded = self
            .hgnc_id_to_transcript_ids
            .iter()
            .filter_map(|(_hgnc_id, tx_ids)| {
                let txs = tx_ids
                    .iter()
                    .map(|tx_id| self.transcript_id_to_transcript.get(tx_id).unwrap())
                    .collect_vec();
                let tx_reasons = tx_ids
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, tx)| {
                        let p = GroupParams::new(idx, tx_ids, &txs);
                        let reason = hgnc_grouped_tx_filters
                            .iter()
                            .filter_map(|(f, r)| f(&p).then_some(*r))
                            .fold(BitFlags::<Reason>::default(), |a, b| a | b);
                        let empty = reason.is_empty();
                        if !empty {
                            Some((Identifier::TxId(tx.clone()), reason))
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
            self.mark_discarded(&id, reason)?;
        }
        self.discard()?;

        tracing::info!("… done filtering transcripts in {:?}", start.elapsed());
        Ok(())
    }

    fn filter_transcripts_with_sequence(
        &mut self,
        seq_repo: &SeqRepo,
    ) -> Result<HashMap<TranscriptId, String>, Error> {
        tracing::info!("Filtering transcripts with seqrepo …");
        let start = Instant::now();

        let five_prime_truncated = |tx: &Transcript| -> bool {
            let is_mt = tx.is_on_contig(MITOCHONDRIAL_ACCESSION);
            if tx.protein_coding() && !is_mt {
                tx.genome_builds.iter().any(|(_release, alignment)| {
                    let cds_start = alignment.cds_start;
                    let mut exons = alignment.exons.clone();
                    exons.is_empty() || {
                        exons.sort_unstable_by_key(|e| e.alt_start_i);
                        let first = exons.first().unwrap();
                        Some(first.alt_start_i) == cds_start
                    }
                })
            } else {
                false
            }
        };
        let three_prime_truncated = |tx: &Transcript| -> bool {
            let is_mt = tx.is_on_contig(MITOCHONDRIAL_ACCESSION);
            if tx.protein_coding() && !is_mt {
                tx.genome_builds.iter().any(|(_release, alignment)| {
                    let cds_end = alignment.cds_end;
                    let mut exons = alignment.exons.clone();
                    exons.is_empty() || {
                        exons.sort_unstable_by_key(|e| e.alt_end_i);
                        let last = exons.last().unwrap();
                        Some(last.alt_end_i) == cds_end
                    }
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
                tx.cds_length().map_or(true, |l| l % 3 != 0)
            } else {
                false
            }
        };

        let (discards, keeps): (Vec<_>, HashMap<TranscriptId, String>) = self
            .transcript_id_to_transcript
            .par_iter()
            .partition_map(|(tx_id, tx)| {
                let has_invalid_cds_length = invalid_cds_length(tx);

                let namespace: String = if tx_id.starts_with("ENST") {
                    String::from("Ensembl")
                } else {
                    String::from("NCBI")
                };
                let seq = seq_repo.fetch_sequence(&AliasOrSeqId::Alias {
                    value: (*tx_id).to_string(),
                    namespace: Some(namespace.clone()),
                });
                if let Ok(seq) = seq {
                    let is_mt = tx.is_on_contig(MITOCHONDRIAL_ACCESSION);
                    if seq.is_empty() {
                        return Either::Left((
                            Identifier::TxId(tx_id.clone()),
                            Reason::MissingSequence.into(),
                        ));
                    };
                    // Skip transcript if it is coding and the translated CDS does not have a stop codon.
                    let cds = tx
                        .start_codon
                        .and_then(|start| tx.stop_codon.map(|stop| (start as usize, stop as usize)));
                    let cds = if tx.protein_coding() { cds } else { None };
                    if let Some((cds_start, cds_end)) = cds
                    {
                        let cds_length = cds_end - cds_start;
                        let delta = (3 - (cds_length % 3)).max(cds_end.saturating_sub(seq.len()));
                        let seq = append_poly_a(seq, delta);
                        let tx_seq_to_translate = &seq[cds_start..cds_end];
                        let aa_sequence = translate_cds(
                            tx_seq_to_translate,
                            true,
                            "*",
                            TranslationTable::Standard,
                        ).expect("Translation should work, since the length is guaranteed to be a multiple of 3 at this point");

                        let mut reason = Reason::empty();
                        let has_missing_stop_codon = (!is_mt && !aa_sequence.ends_with('*'))
                            || (is_mt && !aa_sequence.contains('*'));
                        if has_missing_stop_codon {
                            reason |= Reason::MissingStopCodon;
                            if five_prime_truncated(tx) {
                                reason |= Reason::FivePrimeEndTruncated;
                            }
                            if three_prime_truncated(tx) {
                                reason |= Reason::ThreePrimeEndTruncated;
                            }
                            if has_invalid_cds_length || self.fixes.get(&Identifier::TxId(tx_id.clone())).map_or(false, |f| f.contains(Fix::Cds)) {
                                reason |= Reason::InvalidCdsLength;
                            }
                        }
                        if reason.is_empty() {
                            Either::Right((tx_id.clone(), seq))
                        } else {
                            Either::Left((Identifier::TxId(tx_id.clone()), reason))
                        }
                    } else {
                        Either::Right((tx_id.clone(), seq))
                    }
                } else {
                    Either::Left((Identifier::TxId(tx_id.clone()), Reason::MissingSequence.into()))
                }
            });

        for (id, reason) in discards {
            self.mark_discarded(&id, reason.into())?;
        }
        self.discard()?;

        tracing::info!(
            "… done filtering transcripts with seqrepo in {:?}",
            start.elapsed()
        );

        Ok(keeps)
    }

    /// Remove all hgnc ids for which there are no transcripts left.
    fn filter_empty_hgnc_mappings(&mut self) -> Result<(), Error> {
        tracing::info!("Removing empty HGNC mappings …");
        for (hgnc_id, txs) in self.hgnc_id_to_transcript_ids.iter() {
            if txs.is_empty() {
                *self.discards.entry(Identifier::Hgnc(*hgnc_id)).or_default() |=
                    Reason::NoTranscriptLeft;
            }
        }
        self.discard()?;
        tracing::info!("… done removing empty HGNC mappings");
        Ok(())
    }

    fn update_pseudogene_status(&mut self) -> Result<(), Error> {
        let empty = Reason::empty();
        for (hgnc_id, _) in &self.hgnc_id_to_transcript_ids {
            if !self
                .discards
                .get(&Identifier::Hgnc(*hgnc_id))
                .unwrap_or(&empty)
                .contains(Reason::NoTranscriptLeft)
            {
                continue;
            }
            if let Some(gene) = self.hgnc_id_to_gene.get(hgnc_id) {
                if gene
                    .biotype
                    .as_ref()
                    .map_or(false, |bt| bt.contains(&BioType::Pseudogene))
                {
                    *self.discards.entry(Identifier::Hgnc(*hgnc_id)).or_default() |=
                        Reason::Pseudogene;
                }
            }
        }
        Ok(())
    }

    /// Find all genes and transcripts on the given contigs.
    fn find_on_contigs(
        &self,
        contig_accessions: &[&str],
    ) -> (HashSet<HgncId>, HashSet<TranscriptId>) {
        self.transcript_id_to_transcript
            .values()
            .filter(|tx| {
                tx.genome_builds
                    .values()
                    .any(|gb| contig_accessions.contains(&gb.contig.as_str()))
            })
            .map(|tx| {
                let hgnc_id: HgncId = tx.hgnc.as_ref().unwrap().parse().unwrap();
                (
                    hgnc_id,
                    TranscriptId::try_new(tx.id.clone()).expect("Invalid TranscriptId"),
                )
            })
            .unzip()
    }

    fn gather_transcript_stats(&self) -> Result<(usize, usize, usize), Error> {
        let (mut n_mane_select, mut n_mane_plus_clinical, mut n_mt) = (0, 0, 0);
        for tx in self.transcript_id_to_transcript.values() {
            let mut is_mane_select = false;
            let mut is_mane_plus_clinical = false;
            for gb in tx.genome_builds.values() {
                if MITOCHONDRIAL_ACCESSIONS.contains(&gb.contig.as_str()) {
                    n_mt += 1;
                }
                if let Some(tag) = &gb.tag {
                    if tag.contains(&Tag::ManeSelect) {
                        is_mane_select = true;
                    }
                    if tag.contains(&Tag::ManePlusClinical) {
                        is_mane_plus_clinical = true;
                    }
                }
            }
            if is_mane_select {
                n_mane_select += 1;
            }
            if is_mane_plus_clinical {
                n_mane_plus_clinical += 1;
            }
        }
        Ok((n_mt, n_mane_select, n_mane_plus_clinical))
    }

    fn update_transcript_tags(&mut self, transcript_id_to_tags: &HashMap<TranscriptId, Vec<Tag>>) {
        self.transcript_id_to_transcript
            .iter_mut()
            .for_each(|(tx_id, tx)| {
                // transfer MANE-related labels from TSV file
                let tx_id_no_version = tx_id.without_version().unwrap();
                if let Some(tags) = transcript_id_to_tags.get(&tx_id_no_version) {
                    tx.genome_builds.iter_mut().for_each(|(_, alignment)| {
                        if let Some(alignment_tag) = &mut alignment.tag {
                            alignment_tag.extend(tags.iter().cloned());
                            alignment_tag.sort();
                            alignment_tag.dedup();
                        }
                    });
                    self.fixes
                        .entry(Identifier::TxId(tx_id.clone()))
                        .or_default()
                        .insert(Fix::Tags);
                }
            });
    }

    fn fix_transcript_genome_builds(&mut self) {
        self.transcript_id_to_transcript
            .values_mut()
            .for_each(|tx| {
                let n = tx.genome_builds.len();
                tx.genome_builds.retain(|key, _| {
                    matches!(
                        (key.as_str(), self.genome_release),
                        ("GRCh37", GenomeRelease::Grch37) | ("GRCh38", GenomeRelease::Grch38)
                    )
                });
                if n != tx.genome_builds.len() {
                    self.fixes
                        .entry(Identifier::TxId(TranscriptId::try_new(&tx.id).unwrap()))
                        .or_default()
                        .insert(Fix::GenomeBuild);
                }
            });
    }

    fn fix_cds(&mut self) {
        self.transcript_id_to_transcript
            .values_mut()
            .filter(|tx| tx.protein_coding())
            .filter_map(|tx| {
                tx.start_codon
                    .and_then(|start| tx.stop_codon.map(|stop| (start, stop)))
                    .and_then(|(cds_start, cds_end)| {
                        let cds_len = cds_end - cds_start;
                        (cds_len % 3 != 0).then_some((tx, cds_start, cds_end, cds_len))
                    })
            })
            .for_each(|(tx, _cds_start, cds_end, cds_len)| {
                assert_eq!(tx.genome_builds.len(), 1);
                for gb in tx.genome_builds.values_mut() {
                    let delta = 3 - (cds_len % 3);
                    if delta == 0 {
                        continue;
                    };
                    tx.stop_codon = Some(cds_end + delta);
                    let exon = gb.exons.iter_mut().max_by_key(|g| g.alt_cds_end_i).unwrap();
                    exon.alt_cds_end_i += delta;
                    exon.cigar.push_str(&format!("{}I", delta));
                    self.fixes
                        .entry(Identifier::TxId(TranscriptId::try_new(&tx.id).unwrap()))
                        .or_default()
                        .insert(Fix::Cds);
                }
            });
    }

    fn symbols_to_id(&self, gene_symbols: &[String]) -> Vec<Identifier> {
        gene_symbols
            .iter()
            .map(|symbol| {
                let (hgnc_id, _) = self
                    .hgnc_id_to_gene
                    .iter()
                    .find(|(_, gene)| gene.gene_symbol == Some(symbol.clone()))
                    .unwrap_or_else(|| panic!("Gene symbol not found: {}", symbol));
                Identifier::Hgnc(*hgnc_id)
            })
            .collect()
    }

    fn gene_name(&self, id: &Identifier) -> Option<String> {
        match id {
            Identifier::Hgnc(hgnc_id) => self
                .hgnc_id_to_gene
                .get(hgnc_id)
                .and_then(|gene| gene.gene_symbol.clone()),
            Identifier::TxId(tx_id) => self
                .transcript_id_to_transcript
                .get(tx_id)
                .and_then(|tx| tx.gene_name.clone()),
        }
    }

    fn tags(&self, id: &Identifier) -> Option<Vec<Tag>> {
        match id {
            Identifier::TxId(tx_id) => self.transcript_id_to_transcript.get(tx_id).map(|tx| {
                tx.genome_builds
                    .values()
                    .flat_map(|gb| gb.tag.clone())
                    .flatten()
                    .sorted_unstable()
                    .dedup()
                    .collect()
            }),
            Identifier::Hgnc(hgnc_id) => {
                self.hgnc_id_to_transcript_ids.get(hgnc_id).map(|tx_ids| {
                    tx_ids
                        .iter()
                        .flat_map(|tx_id| {
                            self.tags(&Identifier::TxId(tx_id.clone()))
                                .unwrap_or_default()
                        })
                        .sorted_unstable()
                        .dedup()
                        .collect()
                })
            }
        }
    }

    /// Group transcripts by release and version (sorted descending).
    fn by_release_and_version<'a>(
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

    fn discard(&mut self) -> Result<(), Error> {
        let (n_transcripts_pre, n_hgnc_ids_pre) = (
            self.transcript_id_to_transcript.len(),
            self.hgnc_id_to_transcript_ids.len(),
        );
        type Counter = HashMap<Reason, usize>;
        #[derive(Debug, Clone, Copy, Display, PartialEq, Eq, Hash)]
        enum Kind {
            Transcript,
            HgncId,
        }
        let mut count_groups = HashMap::<Kind, Counter>::new();
        let kinds = [Kind::Transcript, Kind::HgncId];
        for kind in &kinds {
            let mut counter = Counter::new();
            for r in Reason::all().iter() {
                counter.entry(r).or_default();
            }
            count_groups.insert(*kind, counter);
        }

        for (id, reason) in self.discards.clone() {
            if reason.is_empty() {
                tracing::warn!("Empty discard reason for {:#?}", id);
                continue;
            }
            self._discard_id(&id, reason)?;
        }

        for (id, reason) in &self.discards {
            match id {
                Identifier::TxId(_) => {
                    for r in reason.iter() {
                        *count_groups
                            .entry(Kind::Transcript)
                            .or_default()
                            .entry(r)
                            .or_default() += 1;
                    }
                }
                Identifier::Hgnc(_) => {
                    for r in reason.iter() {
                        *count_groups
                            .entry(Kind::HgncId)
                            .or_default()
                            .entry(r)
                            .or_default() += 1;
                    }
                }
            }
        }
        let (n_transcripts_post, n_hgnc_ids_post) = (
            self.transcript_id_to_transcript.len(),
            self.hgnc_id_to_transcript_ids.len(),
        );
        tracing::info!(
            "Discarded {} transcripts and {} HGNC IDs",
            n_transcripts_pre.abs_diff(n_transcripts_post),
            n_hgnc_ids_pre.abs_diff(n_hgnc_ids_post)
        );
        for kind in kinds {
            tracing::info!(
                "{}: {:?}",
                kind,
                Reason::all()
                    .iter()
                    .filter(|r| count_groups[&kind][r] > 0)
                    .map(|r| format!("{:#?}: {}", r, count_groups[&kind][&r]))
                    .join(", ")
            );
        }
        Ok(())
    }

    fn _discard_id(&mut self, id: &Identifier, reason: BitFlags<Reason>) -> Result<(), Error> {
        *self.discards.entry(id.clone()).or_default() |= reason;
        match id {
            Identifier::Hgnc(hgnc_id) => {
                let _gene = self.hgnc_id_to_gene.remove(hgnc_id);
                let txs = self.hgnc_id_to_transcript_ids.remove(hgnc_id);
                for tx_id in txs.unwrap_or_default() {
                    *self
                        .discards
                        .entry(Identifier::TxId(tx_id.clone()))
                        .or_default() |= reason;
                    let _transcript = self.transcript_id_to_transcript.remove(&tx_id);
                }
            }
            Identifier::TxId(tx_id) => {
                let transcript = self.transcript_id_to_transcript.remove(tx_id);
                if let Some(hgnc) = transcript.and_then(|t| t.hgnc) {
                    let hgnc_id: HgncId = hgnc.parse()?;
                    let txs = self.hgnc_id_to_transcript_ids.get_mut(&hgnc_id);
                    if let Some(txs) = txs {
                        txs.retain(|tx| tx != tx_id);
                        if txs.is_empty() {
                            *self.discards.entry(Identifier::Hgnc(hgnc_id)).or_default() |=
                                Reason::NoTranscriptLeft;
                            self.hgnc_id_to_transcript_ids.remove(&hgnc_id);
                            self.hgnc_id_to_gene.remove(&hgnc_id);
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn mark_discarded(&mut self, id: &Identifier, reason: BitFlags<Reason>) -> Result<(), Error> {
        if reason.is_empty() {
            panic!("Empty discard reason for {:#?}", id);
        }
        *self.discards.entry(id.clone()).or_default() |= reason;
        Ok(())
    }

    /// For each transcript that has been discarded for whatever reason, propagate the reason to its
    /// parent hgnc entry *if* the hgnc entry already exists and has a non-empty reason
    /// (to avoid erroneously discarding hgnc entries for a non-important reason).
    fn propagate_discard_reasons(&mut self, raw: &Self) -> Result<(), Error> {
        let discards = self.discards.clone();
        for (tx_id, reason) in discards.iter().filter_map(|(id, reason)| match id {
            Identifier::TxId(tx_id) => Some((tx_id, reason)),
            _ => None,
        }) {
            let tx = &raw.transcript_id_to_transcript[tx_id];
            if let Some(hgnc_id) = tx.hgnc.as_ref() {
                let hgnc_id = hgnc_id.parse()?;
                if let Some(hgnc_reason) = self.discards.get_mut(&Identifier::Hgnc(hgnc_id)) {
                    if !hgnc_reason.is_empty() {
                        *hgnc_reason |= *reason;
                    }
                }
            }
        }
        Ok(())
    }

    /// Perform protobuf file construction.
    ///
    /// This can be done by simply converting the models from ``hgvs-rs`` to the prost generated data structures.
    fn build_protobuf(
        &mut self,
        sequence_map: &mut HashMap<TranscriptId, String>,
        genome_release: GenomeRelease,
    ) -> Result<TxSeqDatabase, Error> {
        tracing::info!("Constructing protobuf data structures …");
        let start = Instant::now();

        let (hgnc_ids, tx_ids): (Vec<HgncId>, Vec<TranscriptId>) = {
            let mut hgnc_ids = self.hgnc_id_to_transcript_ids.keys().cloned().collect_vec();
            hgnc_ids.sort_unstable();
            let tx_ids = self
                .hgnc_id_to_transcript_ids
                .values()
                .flatten()
                .sorted_unstable()
                .dedup()
                .cloned()
                .collect();
            (hgnc_ids, tx_ids)
        };

        // Construct sequence database.
        tracing::info!("  Constructing sequence database …");
        let seq_db = {
            // Insert into protobuf and keep track of pointers in `Vec`s.
            let mut aliases = Vec::new();
            let mut aliases_idx = Vec::new();
            let mut seqs = Vec::new();

            for tx_id in &tx_ids {
                let seq = sequence_map.remove(tx_id).unwrap();

                // Register sequence into protobuf.
                aliases.push((*tx_id).to_string());
                aliases_idx.push(seqs.len() as u32);
                seqs.push(seq);
            }

            // Finalize by creating `SequenceDb`.
            crate::pbs::txs::SequenceDb {
                aliases,
                aliases_idx,
                seqs,
            }
        };

        tracing::info!("  Creating transcript records for each gene…");
        let data_transcripts = {
            let mut data_transcripts = Vec::new();
            // For each gene (in hgnc id order) ...
            for hgnc_id in &hgnc_ids {
                let tx_ids = self
                    .hgnc_id_to_transcript_ids
                    .get(hgnc_id)
                    .unwrap_or_else(|| panic!("No transcripts for hgnc id {:?}", &hgnc_id));

                // ... for each transcript of the gene ...
                for tx_id in tx_ids.iter().sorted_unstable() {
                    let mut tags: Vec<i32> = Vec::new();
                    let tx = self
                        .transcript_id_to_transcript
                        .get(tx_id)
                        .unwrap_or_else(|| panic!("No transcript for id {:?}", tx_id));
                    // ... build genome alignment for selected:
                    let mut genome_alignments = Vec::new();
                    for (genome_build, alignment) in &tx.genome_builds {
                        // obtain basic properties
                        let genome_build = match genome_build.as_ref() {
                            "GRCh37" => crate::pbs::txs::GenomeBuild::Grch37,
                            "GRCh38" => crate::pbs::txs::GenomeBuild::Grch38,
                            _ => panic!("Unknown genome build {:?}", genome_build),
                        };
                        let models::GenomeAlignment {
                            contig,
                            cds_start,
                            cds_end,
                            ..
                        } = alignment.clone();
                        let strand = match alignment.strand {
                            models::Strand::Plus => crate::pbs::txs::Strand::Plus,
                            models::Strand::Minus => crate::pbs::txs::Strand::Minus,
                        };
                        if let Some(tag) = alignment.tag.as_ref() {
                            for t in tag {
                                let elem = match t {
                                    Tag::Basic => crate::pbs::txs::TranscriptTag::Basic.into(),
                                    Tag::EnsemblCanonical => {
                                        crate::pbs::txs::TranscriptTag::EnsemblCanonical.into()
                                    }
                                    Tag::ManeSelect => {
                                        crate::pbs::txs::TranscriptTag::ManeSelect.into()
                                    }
                                    Tag::ManePlusClinical => {
                                        crate::pbs::txs::TranscriptTag::ManePlusClinical.into()
                                    }
                                    Tag::RefSeqSelect => {
                                        crate::pbs::txs::TranscriptTag::RefSeqSelect.into()
                                    }
                                };
                                if !tags.contains(&elem) {
                                    tags.push(elem);
                                }
                            }
                        }
                        // Look into any "note" string for a selenoprotein marker and
                        // add this as a tag.
                        if let Some(note) = alignment.note.as_ref() {
                            let needle = "UGA stop codon recoded as selenocysteine";
                            if note.contains(needle) {
                                tags.push(crate::pbs::txs::TranscriptTag::Selenoprotein.into());
                            }
                        }
                        // and construct vector of all exons
                        let exons: Vec<_> = alignment
                            .exons
                            .iter()
                            .map(|exon| {
                                let models::Exon {
                                    alt_start_i,
                                    alt_end_i,
                                    ord,
                                    alt_cds_start_i,
                                    alt_cds_end_i,
                                    cigar,
                                } = exon.clone();
                                crate::pbs::txs::ExonAlignment {
                                    alt_start_i,
                                    alt_end_i,
                                    ord,
                                    alt_cds_start_i: if alt_cds_start_i == -1 {
                                        None
                                    } else {
                                        Some(alt_cds_start_i)
                                    },
                                    alt_cds_end_i: if alt_cds_end_i == -1 {
                                        None
                                    } else {
                                        Some(alt_cds_end_i)
                                    },
                                    cigar,
                                }
                            })
                            .collect();
                        // and finally push the genome alignment
                        genome_alignments.push(crate::pbs::txs::GenomeAlignment {
                            genome_build: genome_build.into(),
                            contig,
                            cds_start,
                            cds_end,
                            strand: strand.into(),
                            exons,
                        });
                    }

                    // Now, just obtain the basic properties and create a new `data::Transcript`.
                    let Gene {
                        biotype,
                        gene_symbol,
                        ..
                    } = self.hgnc_id_to_gene.get(hgnc_id).unwrap().clone();
                    let biotype = if biotype.unwrap().contains(&BioType::ProteinCoding) {
                        crate::pbs::txs::TranscriptBiotype::Coding.into()
                    } else {
                        crate::pbs::txs::TranscriptBiotype::NonCoding.into()
                    };
                    let Transcript {
                        protein,
                        start_codon,
                        stop_codon,
                        ..
                    } = tx.clone();

                    tags.sort_unstable();
                    tags.dedup();

                    data_transcripts.push(crate::pbs::txs::Transcript {
                        id: (*tx_id).to_string(),
                        gene_symbol: gene_symbol.expect("missing gene symbol"),
                        gene_id: hgnc_id.to_string(),
                        biotype,
                        tags,
                        protein,
                        start_codon,
                        stop_codon,
                        genome_alignments,
                    });
                }
            }

            data_transcripts
        };

        tracing::info!(" … done creating transcripts in {:#?}", start.elapsed());

        // Build mapping of gene HGNC symbol to transcript IDs.
        tracing::info!("  Build gene symbol to transcript ID mapping …");
        let gene_to_tx = self
            .hgnc_id_to_transcript_ids
            .iter()
            .map(|(gene_id, tx_ids)| crate::pbs::txs::GeneToTxId {
                gene_id: gene_id.to_string(),
                tx_ids: tx_ids.iter().map(|tx_id| tx_id.to_string()).collect(),
            })
            .collect::<Vec<_>>();
        tracing::info!(" … done building gene symbol to transcript ID mapping");

        // Compose transcript database from transcripts and gene to transcript mapping.
        tracing::info!("  Composing transcript seq database …");
        let tx_db = crate::pbs::txs::TranscriptDb {
            transcripts: data_transcripts,
            gene_to_tx,
        };

        // Compose the final transcript and sequence database.
        let tx_seq_db = TxSeqDatabase {
            tx_db: Some(tx_db),
            seq_db: Some(seq_db),
            version: Some(crate::common::version().to_string()),
            genome_release: Some(genome_release.name()),
        };

        tracing::info!(" … done composing transcript seq database");

        Ok(tx_seq_db)
    }
}

fn append_poly_a(seq: String, length: usize) -> String {
    // Append poly-A for chrMT transcripts (which are from ENSEMBL).
    // This also potentially fixes the stop codon.
    let mut seq = seq.into_bytes();
    seq.extend_from_slice(b"A".repeat(length).as_slice());
    String::from_utf8(seq).expect("must be valid UTF-8")
}

#[bitflags]
#[repr(u32)]
#[derive(Debug, Clone, Copy, Serialize, Hash, PartialEq, Eq)]
enum Reason {
    Biotype,
    CdsEndAfterSequenceEnd,
    DeselectedGene,
    EmptyGenomeBuilds,
    FivePrimeEndTruncated,
    InvalidCdsLength,
    MissingGene,
    MissingGeneSymbol,
    MissingHgncId,
    MissingSequence,
    MissingStopCodon,
    NoTranscriptLeft,
    NoTranscripts,
    OldVersion,
    OnlyPartialAlignmentInRefSeq,
    PredictedTranscript,
    PredictedTranscriptsOnly,
    Pseudogene,
    ThreePrimeEndTruncated,
    TranscriptPriority,
    UseNmTranscriptInsteadOfNr,
}

#[bitflags]
#[repr(u8)]
#[derive(Debug, Clone, Copy, Serialize, Hash, PartialEq, Eq)]
enum Fix {
    Cds,
    GenomeBuild,
    Tags,
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "type", content = "value")]
enum ReportEntry {
    Discard(Discard),
    Fix(LogFix),
    Log(serde_json::Value),
}
#[serde_as]
#[derive(Debug, Clone, Serialize)]
struct Discard {
    source: String,
    #[serde_as(as = "DisplayFromStr")]
    reason: BitFlags<Reason>,
    id: Identifier,
    gene_name: Option<String>,
    tags: Option<Vec<Tag>>,
}

#[serde_as]
#[derive(Debug, Clone, Serialize)]
struct LogFix {
    source: String,
    #[serde_as(as = "DisplayFromStr")]
    fix: BitFlags<Fix>,
    id: Identifier,
    gene_name: Option<String>,
    tags: Option<Vec<Tag>>,
}

fn txid_to_label(
    label_tsv_path: impl AsRef<Path>,
) -> Result<HashMap<TranscriptId, Vec<Tag>>, Error> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_path(label_tsv_path.as_ref())?;

    rdr.deserialize()
        .map(|result| {
            result
                .map_err(anyhow::Error::from)
                .and_then(|entry: LabelEntry| {
                    TranscriptId::try_new(entry.transcript_id)
                        .map(|txid| {
                            (
                                txid,
                                entry
                                    .label
                                    .split(',')
                                    .map(models::str_to_tag)
                                    .collect::<Vec<_>>(),
                            )
                        })
                        .map_err(anyhow::Error::from)
                })
        })
        .collect()
}
fn read_cdot_json(path: impl AsRef<Path>) -> Result<models::Container, Error> {
    Ok(if path.as_ref().extension().unwrap_or_default() == "gz" {
        tracing::info!("(from gzip compressed file)");
        serde_json::from_reader(std::io::BufReader::new(flate2::read::GzDecoder::new(
            File::open(path)?,
        )))?
    } else {
        tracing::info!("(from uncompressed file)");
        serde_json::from_reader(std::io::BufReader::new(File::open(path)?))?
    })
}

/// Create file-backed `SeqRepo`.
fn open_seqrepo(path: impl AsRef<Path>) -> Result<SeqRepo, Error> {
    tracing::info!("Opening seqrepo…");
    let start = Instant::now();
    let seqrepo = PathBuf::from(path.as_ref());
    let p = path.as_ref().to_str();
    let path = seqrepo
        .parent()
        .ok_or(anyhow::anyhow!("Could not get parent from {:?}", &p))?
        .to_str()
        .unwrap()
        .to_string();
    let instance = seqrepo
        .file_name()
        .ok_or(anyhow::anyhow!("Could not get basename from {:?}", &p))?
        .to_str()
        .unwrap()
        .to_string();
    let seqrepo = SeqRepo::new(path, &instance)?;
    tracing::info!("… seqrepo opened in {:?}", start.elapsed());
    Ok(seqrepo)
}

/// Load the cdot JSON files.
fn load_cdot_files(args: &Args) -> Result<TranscriptLoader, Error> {
    tracing::info!("Loading cdot JSON files …");
    let start = Instant::now();
    let labels = args
        .path_mane_txs_tsv
        .as_ref()
        .map(txid_to_label)
        .transpose()?;
    let loaders = args
        .path_cdot_json
        .iter()
        .map(|cdot_path| {
            let mut loader = TranscriptLoader::new(args.genome_release);
            loader.load_cdot(cdot_path).map(|_| loader)
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let mut merged = loaders
        .into_iter()
        .reduce(|mut a, mut b| {
            a.merge(&mut b);
            a
        })
        .unwrap();
    merged.apply_fixes(&labels);

    tracing::info!(
        "… done loading cdot JSON files in {:?} -- #transcripts = {}, #hgnc_ids = {}",
        start.elapsed(),
        merged
            .transcript_id_to_transcript
            .len()
            .separate_with_underscores(),
        merged
            .hgnc_id_to_transcript_ids
            .len()
            .separate_with_underscores()
    );

    Ok(merged)
}

/// Main entry point for `db create txs` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), Error> {
    rayon::ThreadPoolBuilder::default()
        .num_threads(args.threads)
        .build_global()?;

    let mut report_file =
        File::create(format!("{}.report.jsonl", args.path_out.display())).map(BufWriter::new)?;
    let mut report = |r: ReportEntry| -> Result<(), Error> {
        writeln!(report_file, "{}", serde_json::to_string(&r)?)?;
        Ok(())
    };
    tracing::info!(
        "Building transcript and sequence database file\ncommon args: {:#?}\nargs: {:#?}",
        common,
        args
    );

    // Load cdot files …
    let mut tx_data = load_cdot_files(args)?;
    for (id, fix) in tx_data.fixes.iter() {
        report(ReportEntry::Fix(LogFix {
            source: "cdot".into(),
            fix: *fix,
            id: id.clone(),
            gene_name: tx_data.gene_name(id),
            tags: tx_data.tags(id),
        }))?;
    }

    let raw_tx_data = tx_data.clone();
    trace_rss_now();
    report(ReportEntry::Log(json!({
        "source": "cdot",
        "total_transcripts": tx_data.transcript_id_to_transcript.len(),
        "total_hgnc_ids": tx_data.hgnc_id_to_transcript_ids.len()
    })))?;

    // … then remove information for certain genes …
    if let Some(ids) = args
        .gene_symbols
        .as_ref()
        .map(|symbols| tx_data.symbols_to_id(symbols))
    {
        tx_data.filter_selected(&ids)?;
    }

    // … then filter hgnc entries with no transcripts to boot …
    tx_data.filter_initial_hgnc_entries()?;
    // … then filter genes (missing hgnc id and/or symbol) …
    tx_data.filter_genes()?;
    // … then filter transcripts …
    tx_data.filter_transcripts()?;
    // … ensure there are no hgnc keys without associated transcripts left …
    tx_data.filter_empty_hgnc_mappings()?;

    // Open seqrepo …
    let seqrepo = open_seqrepo(&args.path_seqrepo_instance)?;
    // … and filter transcripts based on their sequences,
    // e.g. checking whether their translation contains a stop codon …
    let mut sequence_map = tx_data.filter_transcripts_with_sequence(&seqrepo)?;
    // … and, again, ensure there are no hgnc keys without associated transcripts left …
    tx_data.filter_empty_hgnc_mappings()?;
    // … if there are genes with no transcripts left, check whether they are pseudogenes …
    tx_data.update_pseudogene_status()?;

    // … and update all discard annotations …
    tx_data.propagate_discard_reasons(&raw_tx_data)?;
    trace_rss_now();

    report(ReportEntry::Log(json!({
        "source": "cdot_filtered",
        "total_transcripts": tx_data.transcript_id_to_transcript.len(),
        "total_hgnc_ids": tx_data.hgnc_id_to_transcript_ids.len()
    })))?;

    // For some stats, count number of chrMt, MANE Select and MANE Plus Clinical transcripts.
    let (n_mt, n_mane_select, n_mane_plus_clinical) = tx_data.gather_transcript_stats()?;

    report(ReportEntry::Log(json!({
        "source": "cdot_filtered",
        "n_mt": n_mt,
        "n_mane_select": n_mane_select,
        "n_mane_plus_clinical": n_mane_plus_clinical
    })))?;

    // … and finally build protobuf file.
    let tx_db = tx_data.build_protobuf(&mut sequence_map, args.genome_release)?;

    // List all discarded transcripts and genes.
    for (id, reason) in tx_data.discards.into_iter().sorted_unstable() {
        report(ReportEntry::Discard(Discard {
            source: "protobuf".into(),
            reason,
            id: id.clone(),
            gene_name: raw_tx_data.gene_name(&id),
            tags: raw_tx_data.tags(&id),
        }))?;
    }
    trace_rss_now();

    write_tx_db(tx_db, &args.path_out)?;

    tracing::info!("Done building transcript and sequence database file");
    Ok(())
}

fn write_tx_db(tx_db: TxSeqDatabase, path: impl AsRef<Path>) -> Result<(), Error> {
    tracing::info!("Writing out final database …");
    let path = path.as_ref();
    let mut buf = prost::bytes::BytesMut::with_capacity(tx_db.encoded_len());
    tx_db
        .encode(&mut buf)
        .map_err(|e| anyhow!("failed to encode: {}", e))?;
    tracing::info!("  … done constructing final tx and seq database");

    // Write out the final transcript and sequence database.
    tracing::info!("  Writing out final database …");
    // Open file and if necessary, wrap in a decompressor.
    let file = std::fs::File::create(path)
        .map_err(|e| anyhow!("failed to create file {}: {}", path.display(), e))?;
    let ext = path.extension().map(|s| s.to_str());
    let mut writer: Box<dyn Write> = if ext == Some(Some("gz")) {
        Box::new(flate2::write::GzEncoder::new(
            file,
            flate2::Compression::default(),
        ))
    } else if ext == Some(Some("zst")) {
        Box::new(
            zstd::Encoder::new(file, 0)
                .map_err(|e| anyhow!("failed to open zstd encoder for {}: {}", path.display(), e))?
                .auto_finish(),
        )
    } else {
        Box::new(file)
    };
    writer
        .write_all(&buf)
        .map_err(|e| anyhow!("failed to write to {}: {}", path.display(), e))?;
    tracing::info!("  … done writing out final database");

    Ok(())
}

#[cfg(test)]
pub mod test {
    use std::path::{Path, PathBuf};

    use clap_verbosity_flag::Verbosity;
    use itertools::Itertools;
    use temp_testdir::TempDir;

    use crate::common::{Args as CommonArgs, GenomeRelease};
    use crate::db::create::TranscriptLoader;
    use crate::db::dump;

    use super::{run, Args};

    #[test]
    fn filter_transcripts_brca1() -> Result<(), anyhow::Error> {
        let path_tsv = Path::new("tests/data/db/create/txs/txs_main.tsv");
        let mut tx_data = TranscriptLoader::new(GenomeRelease::Grch37);
        let labels = super::txid_to_label(path_tsv)?;
        tx_data.load_cdot(Path::new(
            "tests/data/db/create/txs/cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json",
        ))?;
        tx_data.apply_fixes(&Some(labels));

        eprintln!("{:#?}", &tx_data.hgnc_id_to_transcript_ids);
        insta::assert_yaml_snapshot!(tx_data
            .hgnc_id_to_transcript_ids
            .get(&1100)
            .unwrap()
            .iter()
            .map(|s| s.as_str())
            .sorted_unstable()
            .collect::<Vec<_>>());

        tx_data.filter_transcripts()?;
        insta::assert_yaml_snapshot!(tx_data
            .hgnc_id_to_transcript_ids
            .get(&1100)
            .unwrap()
            .iter()
            .map(|s| s.as_str())
            .sorted_unstable()
            .collect::<Vec<_>>());

        insta::assert_snapshot!(&tx_data.cdot_version);

        Ok(())
    }

    #[test]
    fn run_smoke_brca1_opa1() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            path_out: tmp_dir.join("out.bin.zst"),
            path_cdot_json: vec![PathBuf::from(
                "tests/data/db/create/txs/cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json",
            )],
            path_mane_txs_tsv: Some(PathBuf::from("tests/data/db/create/txs/txs_main.tsv")),
            path_seqrepo_instance: PathBuf::from("tests/data/db/create/txs/latest"),
            genome_release: GenomeRelease::Grch38,
            max_txs: None,
            gene_symbols: None,
            threads: 1,
        };

        run(&common_args, &args)?;

        let mut buf: Vec<u8> = Vec::new();
        dump::run_with_write(
            &Default::default(),
            &dump::Args {
                path_db: tmp_dir.join("out.bin.zst"),
            },
            &mut buf,
        )?;
        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }

    #[test]
    fn run_smoke_selenoproteins() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            path_out: tmp_dir.join("out.bin.zst"),
            path_cdot_json: vec![PathBuf::from(
                "tests/data/db/create/seleonoproteins/cdot-0.2.22.refseq.grch38.selenon.json",
            )],
            path_mane_txs_tsv: Some(PathBuf::from("tests/data/db/create/txs/txs_main.tsv")),
            path_seqrepo_instance: PathBuf::from("tests/data/db/create/seleonoproteins/latest"),
            genome_release: GenomeRelease::Grch38,
            max_txs: None,
            gene_symbols: None,
            threads: 1,
        };

        run(&common_args, &args)?;

        let mut buf: Vec<u8> = Vec::new();
        dump::run_with_write(
            &Default::default(),
            &dump::Args {
                path_db: tmp_dir.join("out.bin.zst"),
            },
            &mut buf,
        )?;
        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[test]
    fn run_smoke_mitochondrial() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(5, 0),
        };
        let args = Args {
            path_out: tmp_dir.join("out.bin.zst"),
            path_cdot_json: vec![PathBuf::from(
                "tests/data/db/create/mitochondrial/cdot-0.2.23.ensembl.chrMT.grch37.gff3.json",
            )],
            path_mane_txs_tsv: None,
            path_seqrepo_instance: PathBuf::from("tests/data/db/create/mitochondrial/latest"),
            genome_release: GenomeRelease::Grch37,
            max_txs: None,
            gene_symbols: None,
            threads: 1,
        };

        run(&common_args, &args)?;

        let mut buf: Vec<u8> = Vec::new();
        dump::run_with_write(
            &Default::default(),
            &dump::Args {
                path_db: tmp_dir.join("out.bin.zst"),
            },
            &mut buf,
        )?;
        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }
}
