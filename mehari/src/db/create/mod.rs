//! Transcript database.

use crate::annotate::seqvars::ann::FeatureTag;
use crate::common::trace_rss_now;
use crate::pbs::txs::{SourceVersion, TxSeqDatabase};
use anyhow::{Error, anyhow};
use cli::Args;
use enumflags2::{BitFlag, BitFlags};
use hgvs::data::cdot::json::models as cdot_models;
use hgvs::data::cdot::json::models::{BioType, Gene, Tag, Transcript};
use itertools::Itertools;
use models::{GeneId, Identifier, LabelEntry, ReportEntry};
use once_cell::sync::Lazy;
use prost::Message;
use serde_json::json;
use serde_with::DisplayFromStr;
use std::cmp::{PartialEq, Reverse};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::{io::Write, time::Instant};
use strum::Display;
use thousands::Separable;

mod build;
mod cdot;
pub mod cli;
mod filter;
mod gff3;
pub mod models;
pub mod reference;

use crate::db::create::build::build_protobuf;
use crate::db::create::cdot::load_cdot;
use crate::db::create::filter::{
    filter_empty_gene_id_mappings, filter_genes, filter_initial_gene_id_entries,
    filter_transcripts, filter_transcripts_with_sequence,
};
use crate::db::create::gff3::load_gff3;
use models::*;

/// Mitochondrial accessions.
const MITOCHONDRIAL_ACCESSIONS: &[&str] = &["NC_012920.1", "NC_001807.4"];

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
                .is_some_and(|bt| bt.iter().any(|b| b.is_protein_coding()))
    }

    fn is_on_contig(&self, contig: &str) -> bool {
        self.genome_builds.values().any(|gb| gb.contig == contig)
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

#[derive(Debug, Clone, Default)]
struct TranscriptLoader {
    genome_release: String,
    annotation_version: String,
    transcript_id_to_transcript: HashMap<TranscriptId, Transcript>,
    gene_id_to_gene: HashMap<GeneId, Gene>,
    gene_id_to_transcript_ids: HashMap<GeneId, Vec<TranscriptId>>,
    discards: HashMap<Identifier, BitFlags<Reason>>,
    fixes: HashMap<Identifier, BitFlags<Fix>>,
    disable_filters: bool,
}

impl TranscriptLoader {
    fn new(genome_release: String, disable_filters: bool) -> Self {
        Self {
            genome_release,
            disable_filters,
            ..Default::default()
        }
    }

    fn merge(&mut self, other: &mut Self) -> &mut Self {
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

    fn apply_fixes(
        &mut self,
        transcript_id_to_tags: &Option<HashMap<TranscriptId, Vec<FeatureTag>>>,
    ) {
        self.fix_transcript_genome_builds();
        if let Some(txid_to_label) = transcript_id_to_tags {
            self.update_transcript_tags(txid_to_label);
        }
        self.fix_cds();
    }

    fn update_pseudogene_status(&mut self) -> Result<(), Error> {
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

    fn update_transcript_tags(
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

    fn fix_transcript_genome_builds(&mut self) {
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
                    let exon = gb
                        .exons
                        .iter_mut()
                        .max_by_key(|g| g.alt_cds_end_i)
                        .expect("No exons found during fix_cds");
                    exon.alt_cds_end_i += delta;
                    exon.cigar.push_str(&format!("{}I", delta));
                    self.fixes
                        .entry(Identifier::Transcript(
                            TranscriptId::try_new(&tx.id).unwrap(),
                        ))
                        .or_default()
                        .insert(Fix::Cds);
                }
            });
    }

    fn gene_name(&self, id: &Identifier) -> Option<String> {
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

    fn tags(&self, id: &Identifier) -> Option<Vec<Tag>> {
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

    fn discard(&mut self, remove: bool) -> Result<(), Error> {
        let (n_transcripts_pre, n_gene_ids_pre) = (
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
        let (n_transcripts_post, n_gene_ids_post) = (
            self.transcript_id_to_transcript.len(),
            self.gene_id_to_transcript_ids.len(),
        );
        tracing::info!(
            "Discarded {} transcripts and {} Gene IDs",
            n_transcripts_pre.abs_diff(n_transcripts_post),
            n_gene_ids_pre.abs_diff(n_gene_ids_post)
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

    fn mark_discarded(&mut self, id: &Identifier, reason: BitFlags<Reason>) -> Result<(), Error> {
        if reason.is_empty() {
            panic!("Empty discard reason for {:#?}", id);
        }
        *self.discards.entry(id.clone()).or_default() |= reason;
        Ok(())
    }

    /// For each transcript that has been discarded for whatever reason, propagate the reason to its
    /// parent gene id entry *if* the gene id entry already exists and has a non-empty reason
    /// (to avoid erroneously discarding gene id entries for a non-important reason).
    fn propagate_discard_reasons(&mut self, _raw: &Self) -> Result<(), Error> {
        // First check whether all transcripts of a gene have been marked as discarded.
        for (gene_id, _) in self.gene_id_to_gene.iter() {
            let tx_ids = self
                .gene_id_to_transcript_ids
                .get(gene_id)
                .map(|v| v.as_slice())
                .unwrap_or_default();
            if !tx_ids.is_empty()
                && tx_ids.iter().all(|tx_id| {
                    self.discards
                        .get(&Identifier::Transcript(tx_id.clone()))
                        .is_some_and(|d| d.intersects(Reason::hard()))
                })
            {
                *self
                    .discards
                    .entry(Identifier::Gene(gene_id.clone()))
                    .or_default() |= Reason::NoTranscriptLeft;
            }
        }

        Ok(())
    }
}

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

fn append_poly_a(seq: String, length: usize) -> String {
    // Append poly-A for chrMT transcripts (which are from ENSEMBL).
    // This also potentially fixes the stop codon.
    let mut seq = seq.into_bytes();
    seq.extend_from_slice(b"A".repeat(length).as_slice());
    String::from_utf8(seq).expect("must be valid UTF-8")
}

fn txid_to_label(
    label_tsv_path: impl AsRef<Path>,
) -> Result<HashMap<TranscriptId, Vec<FeatureTag>>, Error> {
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
                                    .map(|s| {
                                        let cdot_tag = cdot_models::str_to_tag(s);
                                        FeatureTag::from(cdot_tag)
                                    })
                                    .collect::<Vec<_>>(),
                            )
                        })
                        .map_err(anyhow::Error::from)
                })
        })
        .collect()
}
pub(crate) fn read_cdot_json(path: impl AsRef<Path>) -> Result<cdot_models::Container, Error> {
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

/// Load the annotations (JSON or GFF3).
fn load_annotations(args: &Args) -> Result<TranscriptLoader, Error> {
    tracing::info!("Loading annotations …");
    let start = Instant::now();
    let labels = args
        .mane_transcripts
        .as_ref()
        .map(txid_to_label)
        .transpose()?;
    let loaders = args
        .annotation
        .iter()
        .map(|path| {
            let mut loader = TranscriptLoader::new(args.assembly.clone(), args.disable_filters);

            let ext = path.extension().unwrap_or_default().to_string_lossy();
            let file_stem = path.file_stem().unwrap_or_default().to_string_lossy();
            let is_gff3 = ext.ends_with("gff")
                || ext.ends_with("gff3")
                || (ext == "gz" && (file_stem.ends_with("gff3") || file_stem.ends_with("gff")));

            if is_gff3 {
                load_gff3(&mut loader, path).map(|_| loader)
            } else {
                load_cdot(&mut loader, path).map(|_| loader)
            }
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
        "… done loading annotations in {:?} -- #transcripts = {}, #gene_ids = {}",
        start.elapsed(),
        merged
            .transcript_id_to_transcript
            .len()
            .separate_with_underscores(),
        merged
            .gene_id_to_transcript_ids
            .len()
            .separate_with_underscores()
    );

    Ok(merged)
}

/// Main entry point for `db create txs` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), Error> {
    fn _run(common: &crate::common::Args, args: &Args) -> Result<(), Error> {
        let mut report_file =
            File::create(format!("{}.report.jsonl", args.output.display())).map(BufWriter::new)?;
        let mut report = |r: ReportEntry| -> Result<(), Error> {
            writeln!(report_file, "{}", serde_json::to_string(&r)?)?;
            Ok(())
        };
        tracing::info!(
            "Building transcript and sequence database file\ncommon args: {:#?}\nargs: {:#?}",
            common,
            args
        );

        let mut tx_data = load_annotations(args)?;
        for (id, fix) in tx_data.fixes.iter() {
            report(ReportEntry::Fix(LogFix {
                source: "annotations".into(),
                fix: *fix,
                id: id.clone(),
                gene_name: tx_data.gene_name(id),
                tags: tx_data.tags(id),
            }))?;
        }

        let raw_tx_data = tx_data.clone();
        trace_rss_now();
        report(ReportEntry::Log(json!({
            "source": "annotations",
            "total_transcripts": tx_data.transcript_id_to_transcript.len(),
            "total_gene_ids": tx_data.gene_id_to_transcript_ids.len()
        })))?;

        let tx_data_ = &mut tx_data;

        // … then filter gene id entries with no transcripts to boot …
        filter_initial_gene_id_entries(tx_data_)?;
        // … then filter genes (missing gene id and/or symbol) …
        filter_genes(tx_data_)?;
        // … then filter transcripts …
        filter_transcripts(tx_data_)?;
        // … ensure there are no gene keys without associated transcripts left …
        filter_empty_gene_id_mappings(tx_data_)?;

        // Open seqrepo / FASTA …
        let mut seq_provider = reference::open_sequence_provider(args)?;
        // … and filter transcripts based on their sequences,
        // e.g. checking whether their translation contains a stop codon …
        let mut sequence_map = filter_transcripts_with_sequence(tx_data_, &mut seq_provider)?;
        filter_empty_gene_id_mappings(tx_data_)?;
        // … if there are genes with no transcripts left, check whether they are pseudogenes …
        tx_data_.update_pseudogene_status()?;

        // … trigger the discard routine, but do not remove anything, just make sure they are consistent,
        // and report the stats.
        let remove = false;
        tx_data_.discard(remove)?;

        // … and update all discard annotations …
        tx_data_.propagate_discard_reasons(&raw_tx_data)?;

        trace_rss_now();

        report(ReportEntry::Log(json!({
            "source": "annotations_filtered",
            "total_transcripts": tx_data_.transcript_id_to_transcript.len(),
            "total_gene_ids": tx_data_.gene_id_to_transcript_ids.len()
        })))?;

        // For some stats, count number of chrMt, MANE Select and MANE Plus Clinical transcripts.
        let (n_mt, n_mane_select, n_mane_plus_clinical) = tx_data_.gather_transcript_stats()?;

        report(ReportEntry::Log(json!({
            "source": "annotations_filtered",
            "n_mt": n_mt,
            "n_mane_select": n_mane_select,
            "n_mane_plus_clinical": n_mane_plus_clinical
        })))?;

        let source_version = bundle_source_version_information(args);

        // … and finally construct protobuf txdb data structures.
        let tx_db = build_protobuf(tx_data_, &mut sequence_map, source_version)?;

        // List all discarded transcripts and genes.
        for (id, reason) in tx_data.discards.into_iter().sorted_unstable() {
            if reason.intersects(Reason::hard()) {
                // if it has at least one hard reason → actually discarded.
                report(ReportEntry::Discard(Discard {
                    source: "protobuf".into(),
                    reason,
                    id: id.clone(),
                    gene_name: raw_tx_data.gene_name(&id),
                    tags: raw_tx_data.tags(&id),
                }))?;
            } else {
                // only soft reasons → kept but flagged
                report(ReportEntry::SoftFilter(SoftFilter {
                    source: "protobuf".into(),
                    reason,
                    id: id.clone(),
                    gene_name: raw_tx_data.gene_name(&id),
                    tags: raw_tx_data.tags(&id),
                }))?;
            }
        }
        trace_rss_now();

        write_tx_db(tx_db, &args.output, args.compression_level)?;

        tracing::info!("Done building transcript and sequence database file");
        Ok(())
    }

    fn bundle_source_version_information(args: &Args) -> SourceVersion {
        let assembly = args.assembly.clone();

        let assembly_version = args.assembly_version.clone();

        let source_name = args.transcript_source.clone();

        let source_version = args.transcript_source_version.clone().unwrap_or("".into());
        let annotation_version = args.annotation_version.clone().unwrap_or("".into());
        let annotation_name = args
            .annotation
            .iter()
            .map(|p| p.file_stem().unwrap().to_string_lossy())
            .collect::<Vec<_>>()
            .join(",");

        SourceVersion {
            mehari_version: crate::common::version().to_string(),
            assembly: assembly.to_string(),
            assembly_version,
            source_name: source_name.to_string(),
            source_version,
            annotation_version,
            annotation_name,
        }
    }

    let threadpool = rayon::ThreadPoolBuilder::default()
        .num_threads(args.threads)
        .build()?;

    threadpool.install(|| _run(common, args))
}

pub(crate) fn write_tx_db(
    tx_db: TxSeqDatabase,
    path: impl AsRef<Path>,
    compression_level: i32,
) -> Result<(), Error> {
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
    let mut writer: Box<dyn Write> = if ext == Some(Some("zst")) {
        Box::new(
            zstd::Encoder::new(file, compression_level)
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
    use rstest::rstest;
    use temp_testdir::TempDir;

    use crate::common::Args as CommonArgs;
    use crate::db::create::TranscriptLoader;
    use crate::db::create::cdot::load_cdot;
    use crate::db::create::cli::Args;
    use crate::db::create::filter::filter_transcripts;
    use crate::db::create::models::GeneId;
    use crate::db::dump;

    use super::run;

    #[test]
    fn filter_transcripts_brca1() -> Result<(), anyhow::Error> {
        let path_tsv = Path::new("tests/data/db/create/txs/txs_main.tsv");
        let mut tx_data = TranscriptLoader::new("GRCh37".to_string(), false);
        let labels = super::txid_to_label(path_tsv)?;
        load_cdot(
            &mut tx_data,
            Path::new("tests/data/db/create/txs/cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json"),
        )?;
        tx_data.apply_fixes(&Some(labels));

        eprintln!("{:#?}", &tx_data.gene_id_to_transcript_ids);
        insta::assert_yaml_snapshot!(
            tx_data
                .gene_id_to_transcript_ids
                .get(&GeneId::Hgnc(1100))
                .unwrap()
                .iter()
                .map(|s| s.as_str())
                .sorted_unstable()
                .collect::<Vec<_>>()
        );

        filter_transcripts(&mut tx_data)?;
        insta::assert_yaml_snapshot!(
            tx_data
                .gene_id_to_transcript_ids
                .get(&GeneId::Hgnc(1100))
                .unwrap()
                .iter()
                .map(|s| s.as_str())
                .sorted_unstable()
                .collect::<Vec<_>>()
        );

        insta::assert_snapshot!(&tx_data.annotation_version);

        Ok(())
    }

    #[rstest]
    #[case("grch37")]
    #[case("grch38")]
    fn run_smoke_brca1_opa1(#[case] assembly: String) -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            output: tmp_dir.join("out.bin.zst"),
            annotation: vec![PathBuf::from(
                "tests/data/db/create/txs/cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json",
            )],
            mane_transcripts: Some(PathBuf::from("tests/data/db/create/txs/txs_main.tsv")),
            seqrepo: Some(PathBuf::from("tests/data/db/create/txs/latest")),
            transcript_sequences: None,
            assembly: assembly.clone(),
            assembly_version: None,
            transcript_source: "refseq".to_string(),
            transcript_source_version: None,
            disable_filters: false,
            threads: 1,
            annotation_version: Some("0.2.22".to_string()),
            compression_level: 19,
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
        crate::common::set_snapshot_suffix!("{}", assembly.to_lowercase());
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
            output: tmp_dir.join("out.bin.zst"),
            annotation: vec![PathBuf::from(
                "tests/data/db/create/seleonoproteins/cdot-0.2.22.refseq.grch38.selenon.json",
            )],
            mane_transcripts: Some(PathBuf::from("tests/data/db/create/txs/txs_main.tsv")),
            seqrepo: Some(PathBuf::from("tests/data/db/create/seleonoproteins/latest")),
            transcript_sequences: None,
            assembly: "grch38".to_string(),
            assembly_version: None,
            transcript_source: "refseq".to_string(),
            transcript_source_version: None,
            disable_filters: false,
            threads: 1,
            annotation_version: Some("0.2.22".to_string()),
            compression_level: 19,
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
            output: tmp_dir.join("out.bin.zst"),
            annotation: vec![PathBuf::from(
                "tests/data/db/create/mitochondrial/cdot-0.2.23.ensembl.chrMT.grch37.gff3.json",
            )],
            mane_transcripts: None,
            seqrepo: Some(PathBuf::from("tests/data/db/create/mitochondrial/latest")),
            transcript_sequences: None,
            assembly: "grch37".to_string(),
            assembly_version: None,
            transcript_source: "ensembl".to_string(),
            transcript_source_version: Some("98".into()),
            disable_filters: false,
            threads: 1,
            annotation_version: Some("0.2.23".to_string()),
            compression_level: 19,
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
