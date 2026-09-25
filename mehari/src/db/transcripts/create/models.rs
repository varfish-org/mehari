use crate::annotate::seqvars::consequence::terms::FeatureTag;
use crate::db::transcripts::create::DisplayFromStr;
use crate::db::transcripts::create::filter::MITOCHONDRIAL_ACCESSIONS;
use anyhow::{Error, anyhow};
use enumflags2::bitflags;
use enumflags2::{BitFlag, BitFlags};
use hgvs::data::cdot::json::models::{BioType, Exon, Gene, Tag, Transcript};
use itertools::Itertools;
use nutype::nutype;
use serde::Serialize;
use serde_with::serde_as;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::str::FromStr;
use strum::Display;

/// Helper struct for parsing the label TSV file.
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize)]
pub struct LabelEntry {
    /// Transcript identifier without version.
    pub(crate) transcript_id: String,
    /// Transcript version.
    transcript_version: usize,
    /// Gene symbol (unused).
    _gene_symbol: String,
    /// Label to transfer.
    pub(crate) label: String,
}

#[nutype(
    sanitize(trim),
    validate(not_empty),
    derive(
        Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, AsRef, Deref, Borrow, Into,
        Display
    )
)]
pub struct TranscriptId(String);

impl TranscriptId {
    pub(crate) fn without_version(&self) -> Result<Self, Error> {
        Ok(Self::try_new(
            self.as_ref()
                .split('.')
                .next()
                .ok_or(anyhow!("No version to split off"))?,
        )?)
    }

    /// Split into accession and version. An ID without a numeric version has version 0.
    pub(crate) fn split_version(&self) -> (&str, u32) {
        self.rsplit_once('.')
            .and_then(|(ac, version)| Some((ac, version.parse().ok()?)))
            .unwrap_or((self.as_ref(), 0))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Hash)]
#[serde(tag = "type", content = "value")]
pub enum Identifier {
    Gene(GeneId),
    Transcript(TranscriptId),
}

#[bitflags]
#[repr(u32)]
#[derive(Debug, Clone, Copy, Serialize, Hash, PartialEq, Eq)]
pub enum Reason {
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
    CdsStartOrEndNotConfirmed,
}

impl Reason {
    /// Reasons that make a transcript completely unusable (Hard Filter).
    pub fn hard() -> BitFlags<Reason> {
        Reason::MissingSequence
            | Reason::EmptyGenomeBuilds
            | Reason::DeselectedGene
            | Reason::NoTranscripts
            | Reason::NoTranscriptLeft
    }

    #[allow(dead_code)]
    /// Reasons that mean the transcript is flawed but salvageable (Soft Filter).
    pub fn soft() -> BitFlags<Reason> {
        BitFlags::all() ^ Self::hard()
    }
}

#[bitflags]
#[repr(u8)]
#[derive(Debug, Clone, Copy, Serialize, Hash, PartialEq, Eq)]
pub enum Fix {
    Cds,
    GenomeBuild,
    Tags,
    UnalignedBases,
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "type", content = "value")]
pub enum ReportEntry {
    Discard(Discard),
    SoftFilter(SoftFilter),
    Fix(LogFix),
    Log(serde_json::Value),
}

#[serde_as]
#[derive(Debug, Clone, Serialize)]
pub struct Discard {
    pub source: String,
    #[serde_as(as = "DisplayFromStr")]
    pub reason: BitFlags<Reason>,
    pub id: Identifier,
    pub gene_name: Option<String>,
    pub tags: Option<Vec<Tag>>,
}

#[serde_as]
#[derive(Debug, Clone, Serialize)]
pub struct SoftFilter {
    pub source: String,
    #[serde_as(as = "DisplayFromStr")]
    pub reason: BitFlags<Reason>,
    pub id: Identifier,
    pub gene_name: Option<String>,
    pub tags: Option<Vec<Tag>>,
}

#[serde_as]
#[derive(Debug, Clone, Serialize)]
pub struct LogFix {
    pub source: String,
    #[serde_as(as = "DisplayFromStr")]
    pub fix: BitFlags<Fix>,
    pub id: Identifier,
    pub gene_name: Option<String>,
    pub tags: Option<Vec<Tag>>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub enum GeneId {
    Hgnc(usize),
    Gene(String),
}

impl Display for GeneId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            GeneId::Hgnc(id) => write!(f, "HGNC:{}", id),
            GeneId::Gene(id) => write!(f, "GENE:{}", id),
        }
    }
}

impl GeneId {
    /// The ID as written to the database: `HGNC:1100` for an HGNC ID, else the ID without the
    /// `GENE:` prefix, e.g. `ENSG00000182378.15`.
    pub(crate) fn to_db_string(&self) -> String {
        match self {
            GeneId::Hgnc(_) => self.to_string(),
            GeneId::Gene(id) => id.clone(),
        }
    }
}

impl FromStr for GeneId {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(num_str) = s.strip_prefix("HGNC:") {
            num_str.parse::<usize>().map(GeneId::Hgnc)
        } else if let Some(gene_str) = s.strip_prefix("GENE:") {
            Ok(GeneId::Gene(gene_str.to_string()))
        } else {
            s.parse::<usize>().map(GeneId::Hgnc)
        }
    }
}

pub trait TranscriptExt {
    fn cds_length(&self) -> Option<u32>;

    fn protein_coding(&self) -> bool;

    fn is_on_contig(&self, contig: &str) -> bool;

    /// Whether the annotation marks the CDS end as incomplete: GENCODE tag `cds_end_NF` or
    /// RefSeq `partial`.
    fn cds_end_incomplete(&self) -> bool;
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
                .is_some_and(|bt| bt.iter().any(|b| b.is_protein_coding()))
    }

    fn is_on_contig(&self, contig: &str) -> bool {
        self.genome_builds.values().any(|gb| gb.contig == contig)
    }

    fn cds_end_incomplete(&self) -> bool {
        self.partial == Some(1)
            || self.genome_builds.values().any(|gb| {
                gb.tag
                    .as_ref()
                    .is_some_and(|tags| tags.contains(&Tag::Other("cds_end_NF".into())))
            })
    }
}

trait BioTypeExt {
    fn is_protein_coding(&self) -> bool;
}

impl BioTypeExt for BioType {
    fn is_protein_coding(&self) -> bool {
        // see https://www.ensembl.org/info/genome/genebuild/biotypes.html
        // and https://www.ensembl.org/Help/Faq?id=468
        matches!(self, BioType::ProteinCoding | BioType::MRna)
    }
}

#[derive(Debug, Clone, Default)]
pub struct TranscriptLoader {
    pub(crate) genome_release: String,
    pub(crate) annotation_version: String,
    pub(crate) transcript_id_to_transcript: HashMap<TranscriptId, Transcript>,
    pub(crate) gene_id_to_gene: HashMap<GeneId, Gene>,
    pub(crate) gene_id_to_transcript_ids: HashMap<GeneId, Vec<TranscriptId>>,
    pub(crate) discards: HashMap<Identifier, BitFlags<Reason>>,
    pub(crate) fixes: HashMap<Identifier, BitFlags<Fix>>,
    pub(crate) disable_filters: bool,
}

impl TranscriptLoader {
    pub(crate) fn new(genome_release: String, disable_filters: bool) -> Self {
        Self {
            genome_release,
            disable_filters,
            ..Default::default()
        }
    }

    pub(crate) fn merge(&mut self, other: &mut Self) -> &mut Self {
        assert_eq!(self.genome_release, other.genome_release);
        self.disable_filters |= other.disable_filters;

        // Use annotation_version from whichever provides it, preferring non-empty
        if self.annotation_version.is_empty() {
            self.annotation_version
                .clone_from(&other.annotation_version);
        }

        for (gene_id, gene) in other.gene_id_to_gene.drain() {
            if let Some(old_gene) = self.gene_id_to_gene.insert(gene_id.clone(), gene.clone()) {
                tracing::warn!(
                    "Overwriting gene: {}\n{:?},\n{:?}",
                    &gene_id,
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

        for (gene_id, mut tx_ids) in other.gene_id_to_transcript_ids.drain() {
            let ids = self.gene_id_to_transcript_ids.entry(gene_id).or_default();
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

    pub(crate) fn apply_fixes(
        &mut self,
        transcript_id_to_tags: &Option<HashMap<TranscriptId, Vec<FeatureTag>>>,
    ) {
        self.fix_transcript_genome_builds();
        if let Some(txid_to_label) = transcript_id_to_tags {
            self.update_transcript_tags(txid_to_label);
        }
        self.fix_unaligned_bases();
        self.fix_cds();
    }

    pub(crate) fn update_pseudogene_status(&mut self) -> Result<(), Error> {
        let empty = Reason::empty();
        for (gene_id, gene) in &self.gene_id_to_gene {
            if !self
                .discards
                .get(&Identifier::Gene(gene_id.clone()))
                .unwrap_or(&empty)
                .contains(Reason::NoTranscriptLeft)
            {
                continue;
            }
            if gene
                .biotype
                .as_ref()
                .is_some_and(|bt| bt.contains(&BioType::Pseudogene))
            {
                *self
                    .discards
                    .entry(Identifier::Gene(gene_id.clone()))
                    .or_default() |= Reason::Pseudogene;
            }
        }
        Ok(())
    }

    pub(crate) fn gather_transcript_stats(&self) -> Result<(usize, usize, usize), Error> {
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

    pub(crate) fn update_transcript_tags(
        &mut self,
        transcript_id_to_tags: &HashMap<TranscriptId, Vec<FeatureTag>>,
    ) {
        self.transcript_id_to_transcript
            .iter_mut()
            .for_each(|(tx_id, tx)| {
                // transfer MANE-related labels from TSV file
                let tx_id_no_version = tx_id.without_version().unwrap();
                if let Some(tags) = transcript_id_to_tags.get(&tx_id_no_version) {
                    let mut any_change = false;
                    tx.genome_builds.iter_mut().for_each(|(_, alignment)| {
                        // Initialize alignment.tag to Some(Vec::new()) when None
                        if alignment.tag.is_none() {
                            alignment.tag = Some(Vec::new());
                        }

                        if let Some(alignment_tag) = &mut alignment.tag {
                            let tags_to_add: Vec<Tag> = tags
                                .iter()
                                .filter_map(|t| {
                                    let standard = Tag::from(t.clone());
                                    if alignment_tag.contains(&standard) {
                                        None
                                    } else {
                                        let backported = Tag::from(t.to_backported());
                                        if alignment_tag.contains(&backported) {
                                            None
                                        } else {
                                            Some(backported)
                                        }
                                    }
                                })
                                .collect();

                            if !tags_to_add.is_empty() {
                                alignment_tag.extend(tags_to_add);
                                alignment_tag.sort();
                                alignment_tag.dedup();
                                any_change = true;
                            }
                        }
                    });

                    if any_change {
                        self.fixes
                            .entry(Identifier::Transcript(tx_id.clone()))
                            .or_default()
                            .insert(Fix::Tags);
                    }
                }
            });
    }

    pub(crate) fn fix_transcript_genome_builds(&mut self) {
        let genome_release = self.genome_release.clone();
        self.transcript_id_to_transcript
            .values_mut()
            .for_each(|tx| {
                let n = tx.genome_builds.len();
                tx.genome_builds
                    .retain(|key, _| key.eq_ignore_ascii_case(&genome_release));
                if n != tx.genome_builds.len() {
                    self.fixes
                        .entry(Identifier::Transcript(
                            TranscriptId::try_new(&tx.id).unwrap(),
                        ))
                        .or_default()
                        .insert(Fix::GenomeBuild);
                }
            });
    }

    /// Complete the stop codon of a CDS whose length is not a multiple of 3 with `A` bases, as
    /// the poly-A tail does. This applies only if the CDS runs to the transcript end and the
    /// annotation does not mark the CDS end as incomplete.
    pub(crate) fn fix_cds(&mut self) {
        self.transcript_id_to_transcript
            .values_mut()
            .filter(|tx| tx.protein_coding() && !tx.cds_end_incomplete())
            .filter_map(|tx| {
                tx.start_codon
                    .and_then(|start| tx.stop_codon.map(|stop| (start, stop)))
                    .and_then(|(cds_start, cds_end)| {
                        let cds_len = cds_end - cds_start;
                        (cds_len % 3 != 0).then_some((tx, cds_start, cds_end, cds_len))
                    })
            })
            .for_each(|(tx, _cds_start, cds_end, cds_len)| {
                // Handle multiple genome builds safely
                if tx.genome_builds.is_empty() {
                    tracing::warn!(
                        "Transcript {} has no genome builds, skipping CDS fix",
                        tx.id
                    );
                    return;
                }

                for gb in tx.genome_builds.values_mut() {
                    let delta = 3 - (cds_len % 3);
                    // `alt_cds_end_i` is the last transcript position of an exon.
                    let tx_end = gb.exons.iter().map(|exon| exon.alt_cds_end_i).max();
                    if delta == 0 || tx_end != Some(cds_end) {
                        continue;
                    };
                    tx.stop_codon = Some(cds_end + delta);
                    let exon = gb
                        .exons
                        .iter_mut()
                        .max_by_key(|g| g.alt_cds_end_i)
                        .expect("No exons found during fix_cds");
                    exon.alt_cds_end_i += delta;
                    // The padding bases exist in the transcript only. In the CIGAR
                    // convention of the hgvs mapper that is `D`; `I` advances the genome.
                    exon.cigar.push_str(&format!("{}D", delta));
                    self.fixes
                        .entry(Identifier::Transcript(
                            TranscriptId::try_new(&tx.id).unwrap(),
                        ))
                        .or_default()
                        .insert(Fix::Cds);
                }
            });
    }

    pub(crate) fn fix_unaligned_bases(&mut self) {
        for (tx_id, tx) in self.transcript_id_to_transcript.iter_mut() {
            let mut changed = false;
            for alignment in tx.genome_builds.values_mut() {
                changed |= fill_unaligned_bases(&mut alignment.exons);
            }
            if changed {
                self.fixes
                    .entry(Identifier::Transcript(tx_id.clone()))
                    .or_default()
                    .insert(Fix::UnalignedBases);
            }
        }
    }

    pub(crate) fn gene_name(&self, id: &Identifier) -> Option<String> {
        match id {
            Identifier::Gene(gene_id) => self
                .gene_id_to_gene
                .get(gene_id)
                .and_then(|gene| gene.gene_symbol.clone()),
            Identifier::Transcript(tx_id) => self
                .transcript_id_to_transcript
                .get(tx_id)
                .and_then(|tx| tx.gene_name.clone()),
        }
    }

    pub(crate) fn tags(&self, id: &Identifier) -> Option<Vec<Tag>> {
        match id {
            Identifier::Transcript(tx_id) => {
                self.transcript_id_to_transcript.get(tx_id).map(|tx| {
                    tx.genome_builds
                        .values()
                        .flat_map(|gb| gb.tag.clone())
                        .flatten()
                        .sorted_unstable()
                        .dedup()
                        .collect()
                })
            }
            Identifier::Gene(gene_id) => {
                self.gene_id_to_transcript_ids.get(gene_id).map(|tx_ids| {
                    tx_ids
                        .iter()
                        .flat_map(|tx_id| {
                            self.tags(&Identifier::Transcript(tx_id.clone()))
                                .unwrap_or_default()
                        })
                        .sorted_unstable()
                        .dedup()
                        .collect()
                })
            }
        }
    }

    pub(crate) fn discard(&mut self, remove: bool) -> Result<(), Error> {
        let (_n_transcripts_pre, _n_gene_ids_pre) = (
            self.transcript_id_to_transcript.len(),
            self.gene_id_to_transcript_ids.len(),
        );
        type Counter = HashMap<Reason, usize>;
        #[derive(Debug, Clone, Copy, Display, PartialEq, Eq, Hash)]
        enum Kind {
            Transcript,
            GeneId,
        }
        let mut count_groups = HashMap::<Kind, Counter>::new();
        let kinds = [Kind::Transcript, Kind::GeneId];
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
            self._discard_id(&id, reason, remove)?;
        }

        for (id, reason) in &self.discards {
            match id {
                Identifier::Transcript(_) => {
                    for r in reason.iter() {
                        *count_groups
                            .entry(Kind::Transcript)
                            .or_default()
                            .entry(r)
                            .or_default() += 1;
                    }
                }
                Identifier::Gene(_) => {
                    for r in reason.iter() {
                        *count_groups
                            .entry(Kind::GeneId)
                            .or_default()
                            .entry(r)
                            .or_default() += 1;
                    }
                }
            }
        }

        let mut hard_txs = 0;
        let mut soft_txs = 0;
        let mut hard_genes = 0;
        let mut soft_genes = 0;

        for (id, reason) in &self.discards {
            let is_hard = reason.intersects(Reason::hard());
            let is_soft = reason.intersects(Reason::soft());

            match id {
                Identifier::Transcript(_) => {
                    if is_hard {
                        hard_txs += 1;
                    } else if is_soft {
                        soft_txs += 1;
                    }
                }
                Identifier::Gene(_) => {
                    if is_hard {
                        hard_genes += 1;
                    } else if is_soft {
                        soft_genes += 1;
                    }
                }
            }
        }

        tracing::info!(
            "Filtered {} transcripts ({} hard discarded, {} soft flagged)",
            hard_txs + soft_txs,
            hard_txs,
            soft_txs
        );
        tracing::info!(
            "Filtered {} Gene IDs ({} hard discarded, {} soft flagged)",
            hard_genes + soft_genes,
            hard_genes,
            soft_genes
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

    fn _discard_id(
        &mut self,
        id: &Identifier,
        reason: BitFlags<Reason>,
        remove: bool,
    ) -> Result<(), Error> {
        *self.discards.entry(id.clone()).or_default() |= reason;
        match id {
            Identifier::Gene(gene_id) => {
                let txs = if remove {
                    self.gene_id_to_transcript_ids.remove(gene_id)
                } else {
                    self.gene_id_to_transcript_ids.get(gene_id).cloned()
                };
                if remove {
                    self.gene_id_to_gene.remove(gene_id);
                }
                for tx_id in txs.unwrap_or_default() {
                    *self
                        .discards
                        .entry(Identifier::Transcript(tx_id.clone()))
                        .or_default() |= reason;
                    if remove {
                        let _transcript = self.transcript_id_to_transcript.remove(&tx_id);
                    }
                }
            }
            Identifier::Transcript(tx_id) => {
                let transcript = if remove {
                    self.transcript_id_to_transcript.remove(tx_id)
                } else {
                    self.transcript_id_to_transcript.get(tx_id).cloned()
                };

                if let Some(hgnc) = transcript.and_then(|t| t.hgnc) {
                    let gene_id: GeneId = hgnc.parse().expect("Invalid GeneId");

                    if remove {
                        let txs = self.gene_id_to_transcript_ids.get_mut(&gene_id);
                        if let Some(txs) = txs {
                            txs.retain(|tx| tx != tx_id);
                            if txs.is_empty() {
                                *self
                                    .discards
                                    .entry(Identifier::Gene(gene_id.clone()))
                                    .or_default() |= Reason::NoTranscriptLeft;
                                self.gene_id_to_transcript_ids.remove(&gene_id);
                                self.gene_id_to_gene.remove(&gene_id);
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    pub(crate) fn mark_discarded(
        &mut self,
        id: &Identifier,
        reason: BitFlags<Reason>,
    ) -> Result<(), Error> {
        if reason.is_empty() {
            panic!("Empty discard reason for {:#?}", id);
        }
        *self.discards.entry(id.clone()).or_default() |= reason;
        Ok(())
    }
}

/// Adds transcript bases that do not align to the genome to the exon that follows them, as
/// transcript-only bases (`D`). RefSeq aligns some transcripts only in parts, which leaves
/// holes in the transcript coordinates of their exons. The hgvs mapper needs exons that cover
/// the transcript from position 1 without holes: it rejects holes between exons, and it
/// counts transcript positions from the start of the first exon.
///
/// Returns whether any exon changed.
fn fill_unaligned_bases(exons: &mut [Exon]) -> bool {
    let mut in_tx_order: Vec<_> = exons.iter_mut().collect();
    in_tx_order.sort_by_key(|exon| exon.ord);
    let mut changed = false;
    let mut next_tx_start = 1;
    for exon in in_tx_order {
        let unaligned = exon.alt_cds_start_i - next_tx_start;
        if unaligned > 0 {
            // The CIGAR runs in transcript direction, also on `-`.
            exon.cigar = format!("{unaligned}D{}", exon.cigar);
            exon.alt_cds_start_i = next_tx_start;
            changed = true;
        }
        next_tx_start = exon.alt_cds_end_i + 1;
    }
    changed
}

#[cfg(test)]
mod tests {
    use super::*;
    use hgvs::data::interface::TxExonsRecord;
    use hgvs::mapper::alignment::build_tx_cigar;
    use hgvs::mapper::cigar::CigarMapper;

    /// Each exon must still map to the transcript positions of its alignment, and the
    /// exons must cover the transcript from position 1 without holes.
    #[rstest::rstest]
    #[case::hole_plus(1, [(1, 100), (140, 239)], true)]
    #[case::hole_minus(-1, [(1, 100), (140, 239)], true)]
    #[case::leading_plus(1, [(8, 107), (108, 207)], true)]
    #[case::leading_minus(-1, [(8, 107), (108, 207)], true)]
    #[case::aligned(1, [(1, 100), (101, 200)], false)]
    fn unaligned_bases_become_transcript_only_bases(
        #[case] strand: i16,
        #[case] tx_ranges: [(i32, i32); 2],
        #[case] changed: bool,
    ) -> Result<(), anyhow::Error> {
        // Two exons of 100 bases around an intron of 200 bases. `ord` counts in transcript
        // direction, so on `-` the first exon is the one at the genomic end.
        let genomic = if strand == 1 {
            [(1000, 1100), (1300, 1400)]
        } else {
            [(1300, 1400), (1000, 1100)]
        };
        let aligned: Vec<_> = (0..)
            .zip(genomic.into_iter().zip(tx_ranges))
            .map(
                |(ord, ((alt_start_i, alt_end_i), (alt_cds_start_i, alt_cds_end_i)))| Exon {
                    alt_start_i,
                    alt_end_i,
                    ord,
                    alt_cds_start_i,
                    alt_cds_end_i,
                    cigar: "100M".into(),
                },
            )
            .collect();
        let mut exons = aligned.clone();
        let is_changed = fill_unaligned_bases(&mut exons);

        // The alignment as the hgvs mapper sees it, see `Provider::get_tx_exons`.
        exons.sort_by_key(|exon| exon.alt_start_i);
        let records: Vec<_> = exons
            .iter()
            .map(|exon| TxExonsRecord {
                alt_start_i: exon.alt_start_i,
                alt_end_i: exon.alt_end_i,
                cigar: exon.cigar.clone(),
                ..Default::default()
            })
            .collect();
        let mapper = CigarMapper::new(&build_tx_cigar(&records, strand)?);
        let tx_len = tx_ranges[1].1;
        assert_eq!(mapper.tgt_len, tx_len);
        // The first genomic base of each exon maps to its transcript position. The mapper
        // counts transcript positions from the genomic start, also on `-`.
        for exon in &aligned {
            let tgt_pos = if strand == 1 {
                exon.alt_cds_start_i - 1
            } else {
                tx_len - exon.alt_cds_end_i
            };
            let mapped = mapper.map_ref_to_tgt(exon.alt_start_i - 1000, "start", true)?;
            assert_eq!(mapped.pos, tgt_pos, "exon {}", exon.ord);
        }

        // The hgvs mapper rejects holes between exons, see `NonAdjacentExons`.
        exons.sort_by_key(|exon| exon.ord);
        assert_eq!(exons[0].alt_cds_start_i, 1);
        assert_eq!(exons[1].alt_cds_start_i, exons[0].alt_cds_end_i + 1);
        assert_eq!(is_changed, changed);

        Ok(())
    }

    /// GFF3 transcripts without `transcript_id` (e.g. RefSeq tRNAs) keep their `ID`,
    /// which has no version.
    #[rstest::rstest]
    #[case("NM_000001.12", ("NM_000001", 12))]
    #[case("TRNAV-CAC-1-1", ("TRNAV-CAC-1-1", 0))]
    #[case("gene.t1", ("gene.t1", 0))]
    fn split_version(#[case] tx_id: &str, #[case] expected: (&str, u32)) -> Result<(), Error> {
        assert_eq!(TranscriptId::try_new(tx_id)?.split_version(), expected);
        Ok(())
    }

    #[rstest::rstest]
    #[case(GeneId::Hgnc(1100), "HGNC:1100")]
    #[case(GeneId::Gene("ENSG00000182378.15".into()), "ENSG00000182378.15")]
    #[case(GeneId::Gene("85358".into()), "85358")]
    fn gene_id_to_db_string(#[case] gene_id: GeneId, #[case] expected: &str) {
        assert_eq!(gene_id.to_db_string(), expected);
    }
}
