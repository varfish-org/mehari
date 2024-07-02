//! Transcript database.

use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::{io::Write, path::PathBuf, time::Instant};

use anyhow::{anyhow, Error};
use clap::Parser;
use derive_new::new;
use enumflags2::{bitflags, BitFlags};
use hgvs::data::cdot::json::models;
use hgvs::data::cdot::json::models::{Gene, Tag, Transcript};
use hgvs::sequences::{translate_cds, TranslationTable};
use indexmap::{IndexMap, IndexSet};
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use nom::ToUsize;
use nutype::nutype;
use once_cell::sync::Lazy;
use prost::Message;
use seqrepo::{AliasOrSeqId, Interface, SeqRepo};
use serde::Serialize;
use serde_json::json;
use serde_with::serde_as;
use serde_with::DisplayFromStr;
use thousands::Separable;

use crate::common::{trace_rss_now, GenomeRelease};
use crate::pbs::txs::TxSeqDatabase;

/// Progress bar style to use.
pub static PROGRESS_STYLE: Lazy<ProgressStyle> = Lazy::new(|| {
    ProgressStyle::with_template(
        "[{elapsed_precise}] [{wide_bar:.cyan/blue}] {human_pos}/{human_len} ({eta})",
    )
    .unwrap()
});

/// Mitochondrial accessions.
const MITOCHONDRIAL_ACCESSIONS: &[&str] = &["NC_012920.1"];

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
    /// Maximal number of transcripts to process.
    #[arg(long)]
    pub max_txs: Option<u32>,
    /// Limit transcript database to the following HGNC symbols.  Useful for
    /// building test databases.
    #[arg(long)]
    pub gene_symbols: Option<Vec<String>>,
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
struct GeneId(String);

#[nutype(
    sanitize(trim),
    validate(not_empty),
    derive(
        Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, AsRef, Deref, Borrow, Into,
        Display
    )
)]
struct GeneSymbol(String);

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
        Ok(Self::new(
            self.as_ref()
                .split('.')
                .next()
                .ok_or(anyhow!("No version to split off"))?,
        )?)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Serialize, Hash)]
#[serde(tag = "type", content = "value")]
enum Identifier {
    Hgnc(HgncId),
    TxId(TranscriptId),
    GeneId(GeneId),
}
trait TranscriptExt {
    fn cds_length(&self) -> Option<u32>;
}
impl TranscriptExt for Transcript {
    fn cds_length(&self) -> Option<u32> {
        self.start_codon
            .and_then(|start| self.stop_codon.map(|stop| stop.abs_diff(start)))
    }
}

#[derive(Debug, Clone, Default)]
struct TranscriptLoader {
    genome_release: GenomeRelease,
    cdot_version: String,
    transcripts: IndexMap<TranscriptId, Transcript>,
    hgnc_id_to_gene_id: IndexMap<HgncId, GeneId>,
    hgnc_id_to_transcript_ids: IndexMap<HgncId, Vec<TranscriptId>>,
    gene_id_to_gene: IndexMap<GeneId, Gene>,
    discards: IndexMap<Identifier, BitFlags<Reason>>,
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
        self.hgnc_id_to_gene_id
            .extend(other.hgnc_id_to_gene_id.drain(..));
        self.gene_id_to_gene.extend(other.gene_id_to_gene.drain(..));
        self.transcripts.extend(other.transcripts.drain(..));
        self.hgnc_id_to_transcript_ids
            .extend(other.hgnc_id_to_transcript_ids.drain(..));
        for (id, reason) in other.discards.drain(..) {
            *self.discards.entry(id).or_default() |= reason;
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
        } = load_cdot_transcripts(path.as_ref())?;
        let cdot_genes = cdot_genes
            .into_iter()
            .map(|(gene_id, gene)| GeneId::new(gene_id).map(|g| (g, gene)))
            .collect::<Result<IndexMap<_, _>, _>>()?;
        let cdot_transcripts = cdot_transcripts
            .into_iter()
            .map(|(txid, tx)| TranscriptId::new(txid).map(|t| (t, tx)))
            .collect::<Result<IndexMap<_, _>, _>>()?;
        self.cdot_version = cdot_version;

        for (gene_id, gene) in cdot_genes {
            self.gene_id_to_gene.insert(gene_id.clone(), gene.clone());
            if let Some(hgnc_id) = gene.hgnc.as_ref().map(|hgnc| hgnc.parse()).transpose()? {
                self.hgnc_id_to_transcript_ids.entry(hgnc_id).or_default();
                self.hgnc_id_to_gene_id.insert(hgnc_id, gene_id.clone());
            }
        }

        for (tx_id, tx) in cdot_transcripts {
            self.transcripts.insert(tx_id.clone(), tx.clone());
            if let Some(hgnc_id) = tx.hgnc.as_ref().map(|hgnc| hgnc.parse()).transpose()? {
                self.hgnc_id_to_transcript_ids
                    .entry(hgnc_id)
                    .or_default()
                    .push(tx_id.clone());
            }
        }
        Ok(())
    }

    fn apply_fixes(&mut self, transcript_id_to_tags: &Option<IndexMap<TranscriptId, Vec<Tag>>>) {
        self.fix_transcript_genome_builds();
        if let Some(txid_to_label) = transcript_id_to_tags {
            self.update_transcript_tags(txid_to_label);
        }
        self.fix_mitochondrial_cds();
    }

    fn filter_genes(&mut self) -> Result<(), Error> {
        tracing::info!("Filtering genes …");
        let missing_hgnc = |_gene_id: &str, gene: &Gene| -> bool {
            gene.hgnc.is_none() || gene.hgnc.as_ref().unwrap().is_empty()
        };
        let filters: [(&dyn Fn(&str, &Gene) -> bool, Reason); 1] =
            [(&missing_hgnc, Reason::MissingHgncId)];

        let discarded_genes = self
            .gene_id_to_gene
            .iter()
            .filter_map(|(gene_id, gene)| {
                let reason = filters
                    .iter()
                    .filter_map(|(f, r)| f(gene_id, gene).then(|| *r))
                    .fold(BitFlags::<Reason>::default(), |a, b| a | b);
                (!reason)
                    .is_empty()
                    .then(|| (Identifier::GeneId(gene_id.clone()), reason))
            })
            .collect_vec();

        for (id, reason) in discarded_genes {
            self.discard_id(&id, reason)?;
        }
        self.discard()?;
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
    fn filter_transcripts(
        &mut self,
        max_genes: Option<u32>,
        selected_hgnc_ids: &Option<Vec<Identifier>>,
    ) -> Result<(), Error> {
        tracing::info!("Filtering transcripts …");

        /// Params used as args for filter functions defined below
        #[derive(new)]
        struct Params<'a> {
            ac: &'a str,
            release: &'a str,
            tx: &'a Transcript,
            seen_ac: &'a IndexSet<(String, String)>,
            seen_nm: bool,
            hgnc_id_to_gene_id: &'a IndexMap<HgncId, GeneId>,
        }

        let missing_hgnc =
            |p: &Params| -> bool { p.tx.hgnc.is_none() || p.tx.hgnc.as_ref().unwrap().is_empty() };
        let deselected_gene = |p: &Params| -> bool {
            !p.hgnc_id_to_gene_id
                .contains_key(&p.tx.hgnc.as_ref().unwrap().parse::<HgncId>().unwrap())
        };
        let empty_genome_builds = |p: &Params| -> bool { p.tx.genome_builds.is_empty() };
        let partial = |p: &Params| -> bool { matches!(p.tx.partial, Some(1)) };

        // Check whether we have already seen the same version of the transcript.
        let old_version = |p: &Params| -> bool {
            p.seen_ac
                .contains(&(p.ac.to_string(), p.release.to_string()))
        };

        // Check whether we have already seen an NM transcript for the gene.
        let nm_instead_of_nr = |p: &Params| -> bool { p.ac.starts_with("NR_") && p.seen_nm };

        // Check whether the transcript is predicted.
        let predicted = |p: &Params| -> bool { p.ac.starts_with('X') };

        // Check transcript's CDS length for being a multiple of 3.
        //
        // Note that the chrMT transcripts have been fixed earlier already to
        // accomodate for how they are fixed by poly-A tailing.
        let invalid_cds_length =
            |p: &Params| -> bool { p.tx.cds_length().map_or(true, |l| l % 3 != 0) };

        type Filter = fn(&Params) -> bool;
        let filters: [(Filter, Reason); 8] = [
            (missing_hgnc, Reason::MissingHgncId),
            (deselected_gene, Reason::DeselectedGene),
            (empty_genome_builds, Reason::EmptyGenomeBuilds),
            (partial, Reason::OnlyPartialAlignmentInRefSeq),
            (old_version, Reason::OldVersion),
            (nm_instead_of_nr, Reason::UseNmTranscriptInsteadOfNr),
            (predicted, Reason::PredictedTranscript),
            (invalid_cds_length, Reason::InvalidCdsLength),
        ];

        let start = Instant::now();

        // Potentially limit number of genes.
        let hgnc_id_to_transcript_ids: Box<dyn Iterator<Item = (&HgncId, &Vec<TranscriptId>)>> =
            if let Some(max_genes) = max_genes {
                tracing::warn!("Limiting to {} genes!", max_genes);
                Box::new(
                    self.hgnc_id_to_transcript_ids
                        .iter()
                        .take(max_genes as usize),
                )
            } else {
                Box::new(self.hgnc_id_to_transcript_ids.iter())
            };

        // Filter map from gene symbol to Vec of chosen transcript identifiers.
        for (hgnc_id, tx_ids) in hgnc_id_to_transcript_ids {
            // Skip transcripts where the gene symbol is not contained in `selected_hgnc_ids`.
            if !selected_hgnc_ids
                .as_ref()
                .map(|ids| ids.contains(&Identifier::Hgnc(*hgnc_id)))
                .unwrap_or(true)
            {
                tracing::trace!("skipping {} / {:?}, because not selected", hgnc_id, tx_ids);
                for tx_id in tx_ids {
                    *self
                        .discards
                        .entry(Identifier::TxId(tx_id.clone()))
                        .or_default() |= Reason::DeselectedGene;
                }
                continue;
            }

            // Only select the highest version of each transcript.
            //
            // First, look for NM transcript,
            // and split off transcript versions from accessions and sort them.
            let seen_nm = tx_ids.iter().any(|tx_id| tx_id.starts_with("NM_"));
            let versioned = Self::split_off_versions(tx_ids);

            // Build `next_tx_ids`.
            let mut seen_ac = IndexSet::new();
            let mut next_tx_ids = Vec::new();

            for (ac, version) in versioned {
                let full_ac = TranscriptId::new(format!("{}.{}", &ac, version))?;
                let ac = ac.to_string();

                let releases = self
                    .transcripts
                    .get(&full_ac)
                    .map(|tx| tx.genome_builds.keys().cloned().collect::<Vec<_>>())
                    .unwrap_or_default();

                for release in releases {
                    let tx = self
                        .transcripts
                        .get(&full_ac)
                        .expect("must exist; accession taken from map earlier");
                    let p = Params::new(
                        &ac,
                        &release,
                        tx,
                        &seen_ac,
                        seen_nm,
                        &self.hgnc_id_to_gene_id,
                    );
                    let reason = filters
                        .iter()
                        .filter_map(|(f, r)| f(&p).then(|| *r))
                        .fold(BitFlags::<Reason>::default(), |a, b| a | b);
                    if !reason.is_empty() {
                        *self
                            .discards
                            .entry(Identifier::TxId(full_ac.clone()))
                            .or_default() |= reason;
                        continue;
                    }
                    // Otherwise, mark transcript as included by storing its accession.
                    next_tx_ids.push(full_ac.clone());
                    seen_ac.insert((ac.clone(), release));
                }
            }

            if next_tx_ids.is_empty() {
                *self.discards.entry(Identifier::Hgnc(*hgnc_id)).or_default() |=
                    Reason::NoTranscript;
            }
        }
        self.discard()?;

        tracing::info!("… done filtering transcripts in {:?}", start.elapsed());
        Ok(())
    }

    fn filter_transcripts_with_sequence(&mut self, seq_repo: &SeqRepo) -> Result<(), Error> {
        tracing::info!("Filtering transcripts with seqrepo …");
        let start = Instant::now();

        let (_, mt_txs) = self.find_on_contigs(MITOCHONDRIAL_ACCESSIONS);
        for (tx_id, tx) in &self.transcripts {
            let cds_length = tx.cds_length().expect("must be some");
            if cds_length % 3 != 0 {
                *self
                    .discards
                    .entry(Identifier::TxId(tx_id.clone()))
                    .or_default() |= Reason::InvalidCdsLength;
                continue;
            }

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
                let is_mt = mt_txs.contains(tx_id);
                let seq = if seq.is_empty() {
                    *self
                        .discards
                        .entry(Identifier::TxId(tx_id.clone()))
                        .or_default() |= Reason::MissingSequence;
                    continue;
                } else {
                    if is_mt {
                        append_poly_a(seq)
                    } else {
                        seq
                    }
                };
                // Skip transcript if it is coding and the translated CDS does not have a stop codon.
                if let Some((cds_start, cds_end)) = tx
                    .start_codon
                    .and_then(|start| tx.stop_codon.map(|stop| (start, stop)))
                {
                    let (cds_start, cds_end) = (cds_start as usize, cds_end as usize);
                    if cds_end > seq.len() {
                        *self
                            .discards
                            .entry(Identifier::TxId(tx_id.clone()))
                            .or_default() |= Reason::CdsEndAfterSequenceEnd;
                        continue;
                    }
                    let tx_seq_to_translate = &seq[cds_start..cds_end];
                    let aa_sequence =
                        translate_cds(tx_seq_to_translate, true, "*", TranslationTable::Standard)?;

                    if (!is_mt && !aa_sequence.ends_with('*'))
                        || (is_mt && !aa_sequence.contains('*'))
                    {
                        *self
                            .discards
                            .entry(Identifier::TxId(tx_id.clone()))
                            .or_default() |= Reason::MissingStopCodon;
                        continue;
                    }
                }
            } else {
                *self
                    .discards
                    .entry(Identifier::TxId(tx_id.clone()))
                    .or_default() |= Reason::MissingSequence;
                continue;
            }
        }
        self.discard()?;

        tracing::info!(
            "… done filtering transcripts with seqrepo in {:?}",
            start.elapsed()
        );

        Ok(())
    }

    fn find_on_contigs(
        &self,
        contig_accessions: &[&str],
    ) -> (IndexSet<GeneId>, IndexSet<TranscriptId>) {
        self.transcripts
            .values()
            .filter(|tx| {
                tx.genome_builds
                    .values()
                    .any(|gb| contig_accessions.contains(&gb.contig.as_str()))
            })
            .map(|tx| {
                let hgnc_id: HgncId = tx.hgnc.as_ref().unwrap().parse().unwrap();
                let gene_id = self.hgnc_id_to_gene_id.get(&hgnc_id).expect(&format!(
                    "No gene ID found for hgnc id {:#?}",
                    tx.hgnc.as_ref()
                ));
                (
                    GeneId::new(gene_id.clone()).expect("Invalid GeneId"),
                    TranscriptId::new(tx.id.clone()).expect("Invalid TranscriptId"),
                )
            })
            .unzip()
    }

    fn gather_transcript_stats(&self) -> Result<(usize, usize, usize), Error> {
        let (mut n_mane_select, mut n_mane_plus_clinical, mut n_mt) = (0, 0, 0);
        for tx in self.transcripts.values() {
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

    fn update_transcript_tags(&mut self, transcript_id_to_tags: &IndexMap<TranscriptId, Vec<Tag>>) {
        self.transcripts.iter_mut().for_each(|(tx_id, tx)| {
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
            }
        });
    }

    fn fix_transcript_genome_builds(&mut self) {
        self.transcripts.values_mut().for_each(|tx| {
            tx.genome_builds.retain(|key, _| {
                matches!(
                    (key.as_str(), self.genome_release),
                    ("GRCh37", GenomeRelease::Grch37) | ("GRCh38", GenomeRelease::Grch38)
                )
            });
        });
    }

    fn fix_mitochondrial_cds(&mut self) {
        self.transcripts.values_mut().for_each(|tx| {
            // fix coding mitochondrial transcripts that have a CDS that is not a multiple of 3
            if let Some(cds_start) = tx.start_codon {
                let cds_end = tx.stop_codon.expect("must be some if start_codon is some");
                let cds_len = cds_end - cds_start;
                if cds_len % 3 != 0 {
                    assert_eq!(
                        tx.genome_builds.len(),
                        1,
                        "only one genome build expected at this point"
                    );
                    let gb = tx.genome_builds.iter_mut().next().unwrap().1;
                    if MITOCHONDRIAL_ACCESSIONS.contains(&gb.contig.as_ref()) {
                        assert_eq!(gb.exons.len(), 1, "only single-exon genes assumed on chrMT");
                        let delta = 3 - cds_len % 3;
                        tx.stop_codon = Some(cds_end + delta);
                        let exon = gb.exons.iter_mut().next().unwrap();
                        exon.alt_cds_end_i += delta;
                        exon.cigar.push_str(&format!("{}I", delta));
                    }
                }
            }
        });
    }

    fn symbols_to_id(&self, gene_symbols: &[String]) -> Vec<Identifier> {
        gene_symbols
            .iter()
            .map(|symbol| {
                let (hgnc_id, _) = self
                    .hgnc_id_to_gene_id
                    .iter()
                    .find(|(_, gene_id)| {
                        self.gene_id_to_gene
                            .get(*gene_id)
                            .map(|g| g.gene_symbol == Some(symbol.clone()))
                            .unwrap_or(false)
                    })
                    .expect(&format!("Gene symbol not found: {}", symbol));
                Identifier::Hgnc(*hgnc_id)
            })
            .collect()
    }

    /// Split off versions from transcript identifiers, then sort descending by version.
    fn split_off_versions(tx_ids: &Vec<TranscriptId>) -> Vec<(&str, u32)> {
        let mut versioned: Vec<_> = tx_ids
            .iter()
            .map(|tx_id| {
                let s: Vec<_> = tx_id.split('.').collect();
                (s[0], s[1].parse::<u32>().expect("invalid version"))
            })
            .collect();
        versioned.sort_unstable_by_key(|(k, v)| (k.to_string(), std::cmp::Reverse(*v)));
        versioned
    }

    fn discard(&mut self) -> Result<(), Error> {
        let (n_transcripts_pre, n_genes_pre, n_hgnc_ids_pre) = (
            self.transcripts.len(),
            self.gene_id_to_gene.len(),
            self.hgnc_id_to_transcript_ids.len(),
        );
        for (id, reason) in self.discards.clone() {
            self.discard_id(&id, reason)?;
        }
        let (n_transcripts_post, n_genes_post, n_hgnc_ids_post) = (
            self.transcripts.len(),
            self.gene_id_to_gene.len(),
            self.hgnc_id_to_transcript_ids.len(),
        );
        tracing::info!(
            "Discarded {} transcripts, {} genes, and {} HGNC IDs",
            n_transcripts_pre - n_transcripts_post,
            n_genes_pre - n_genes_post,
            n_hgnc_ids_pre - n_hgnc_ids_post
        );
        Ok(())
    }

    fn discard_id(&mut self, id: &Identifier, reason: BitFlags<Reason>) -> Result<(), Error> {
        *self.discards.entry(id.clone()).or_default() |= reason;
        match id {
            Identifier::Hgnc(hgnc_id) => {
                let gene_id = self
                    .hgnc_id_to_gene_id
                    .swap_remove(hgnc_id)
                    .ok_or_else(|| anyhow!("Gene not found for {}", hgnc_id))?;
                let _gene = self.gene_id_to_gene.swap_remove(&gene_id);
                let tx_ids = self
                    .hgnc_id_to_transcript_ids
                    .swap_remove(hgnc_id)
                    .ok_or_else(|| anyhow!("No transcripts for {}", hgnc_id))?;
                for tx_id in tx_ids {
                    self.transcripts.swap_remove(&tx_id);
                    *self.discards.entry(Identifier::TxId(tx_id)).or_default() |= reason;
                }
            }
            Identifier::TxId(tx_id) => {
                let _tx = self
                    .transcripts
                    .swap_remove(tx_id)
                    .ok_or_else(|| anyhow!("Transcript not found: {}", tx_id))?;
            }
            Identifier::GeneId(gene_id) => {
                let _gene = self
                    .gene_id_to_gene
                    .swap_remove(gene_id)
                    .ok_or_else(|| anyhow!("Gene not found: {}", gene_id))?;
                if let Some(hgnc_id) = self.hgnc_id_to_gene_id.iter().find_map(|(hgnc_id, g)| {
                    if g == gene_id {
                        Some(hgnc_id)
                    } else {
                        None
                    }
                }) {
                    *self.discards.entry(Identifier::Hgnc(*hgnc_id)).or_default() |= reason;
                    let tx_ids = self
                        .hgnc_id_to_transcript_ids
                        .swap_remove(hgnc_id)
                        .ok_or_else(|| anyhow!("Gene not found: {}", hgnc_id))?;
                    for tx_id in tx_ids {
                        self.transcripts.swap_remove(&tx_id);
                        *self.discards.entry(Identifier::TxId(tx_id)).or_default() |= reason;
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
        seqrepo: &SeqRepo,
        is_silent: bool,
        genome_release: GenomeRelease,
    ) -> Result<TxSeqDatabase, Error> {
        let (_, mt_transcript_ids) = self.find_on_contigs(MITOCHONDRIAL_ACCESSIONS);

        tracing::info!("Constructing protobuf data structures …");
        trace_rss_now();
        let start = Instant::now();

        // Construct sequence database.
        tracing::info!("  Constructing sequence database …");
        let seq_db = {
            // Insert into protobuf and keep track of pointers in `Vec`s.
            let mut aliases = Vec::new();
            let mut aliases_idx = Vec::new();
            let mut seqs = Vec::new();
            let pb = if is_silent {
                ProgressBar::hidden()
            } else {
                ProgressBar::new(self.transcripts.len() as u64)
            };
            pb.set_style(PROGRESS_STYLE.clone());
            for (tx_id, _tx) in &self.transcripts {
                pb.inc(1);
                let is_mt = mt_transcript_ids.contains(tx_id);
                let namespace: String = if tx_id.starts_with("ENST") {
                    String::from("Ensembl")
                } else {
                    String::from("NCBI")
                };
                let res_seq = seqrepo.fetch_sequence(&AliasOrSeqId::Alias {
                    value: (*tx_id).to_string(),
                    namespace: Some(namespace.clone()),
                });
                let seq = if let Ok(seq) = res_seq {
                    if is_mt {
                        append_poly_a(seq)
                    } else {
                        seq
                    }
                } else {
                    continue;
                };

                // Register sequence into protobuf.
                aliases.push((*tx_id).to_string());
                aliases_idx.push(seqs.len() as u32);
                seqs.push(seq);
            }
            pb.finish_and_clear();
            // Finalize by creating `SequenceDb`.
            crate::pbs::txs::SequenceDb {
                aliases,
                aliases_idx,
                seqs,
            }
        };

        trace_rss_now();

        tracing::info!("  Creating transcript records for each gene…");
        let data_transcripts = {
            let hgnc_ids = {
                let mut hgnc_ids: Vec<_> = self.hgnc_id_to_gene_id.keys().cloned().collect();
                hgnc_ids.sort();
                hgnc_ids
            };
            let mut data_transcripts = Vec::new();
            // For each gene (in hgnc id order) ...
            for hgnc_id in &hgnc_ids {
                let gene_id = self.hgnc_id_to_gene_id.get(hgnc_id).unwrap();
                let tx_ids = self
                    .hgnc_id_to_transcript_ids
                    .get(hgnc_id)
                    .unwrap_or_else(|| panic!("No transcripts for hgnc id {:?}", &hgnc_id));
                let tx_ids = tx_ids
                    .iter()
                    .filter(|&tx_id| {
                        self.discards
                            .entry(Identifier::TxId(tx_id.clone()))
                            .or_default()
                            .is_empty()
                    })
                    .collect::<Vec<_>>();
                if tx_ids.is_empty() {
                    *self
                        .discards
                        .entry(Identifier::Hgnc(hgnc_id.clone()))
                        .or_default() |= Reason::NoTranscript;
                    continue;
                }

                // ... for each transcript of the gene ...
                for tx_id in tx_ids {
                    let mut tags: Vec<i32> = Vec::new();
                    let tx = self
                        .transcripts
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
                        hgnc,
                        gene_symbol,
                        ..
                    } = self.gene_id_to_gene.get(gene_id).unwrap().clone();
                    let biotype = if biotype.unwrap().contains(&models::BioType::ProteinCoding) {
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

                    if gene_symbol.is_none() {
                        *self
                            .discards
                            .entry(Identifier::TxId(tx_id.clone()))
                            .or_default() |= Reason::MissingGeneSymbol;
                        continue;
                    }

                    data_transcripts.push(crate::pbs::txs::Transcript {
                        id: (*tx_id).to_string(),
                        gene_symbol: gene_symbol.expect("missing gene symbol"),
                        gene_id: hgnc
                            .map(|s| s.parse::<HgncId>().expect("Invalid HGNC Id"))
                            .expect("missing HGNC ID")
                            .to_string(),
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

        trace_rss_now();

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

        trace_rss_now();

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

        trace_rss_now();

        Ok(tx_seq_db)
    }
}

fn append_poly_a(seq: String) -> String {
    // Append poly-A for chrMT transcripts (which are from ENSEMBL).
    // This also potentially fixes the stop codon.
    let mut seq = seq.into_bytes();
    seq.extend_from_slice(b"A".repeat(300).as_slice());
    String::from_utf8(seq).expect("must be valid UTF-8")
}

#[bitflags]
#[repr(u16)]
#[derive(Debug, Clone, Copy, Serialize)]
enum Reason {
    MissingHgncId,
    // NotMTandMissingMapLocation,
    DeselectedGene,
    EmptyGenomeBuilds,
    OldVersion,
    UseNmTranscriptInsteadOfNr,
    PredictedTranscript,
    InvalidCdsLength,
    NoTranscript,
    MissingSequence,
    CdsEndAfterSequenceEnd,
    MissingStopCodon,
    TranscriptPriority,
    OnlyPartialAlignmentInRefSeq,
    MissingGeneSymbol,
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "type", content = "value")]
enum ReportEntry {
    Discard(Discard),
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

fn txid_to_label(
    label_tsv_path: impl AsRef<Path>,
) -> Result<IndexMap<TranscriptId, Vec<Tag>>, Error> {
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
                    TranscriptId::new(entry.transcript_id)
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
fn load_cdot_transcripts(path: impl AsRef<Path>) -> Result<models::Container, Error> {
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
            loader.load_cdot(cdot_path).unwrap();
            loader
        })
        .collect_vec();
    let mut merged = loaders
        .into_iter()
        .reduce(|mut a, mut b| {
            // because we're loading multiple cdot files,
            // it may happen that there are multiple entries from different sources
            // for the same hgnc id
            // in that case, we report which transcripts are "overwritten" by the merge
            for hgnc_id in b.hgnc_id_to_transcript_ids.keys() {
                if let Some(tx_ids) = a.hgnc_id_to_transcript_ids.get(hgnc_id) {
                    for tx_id in tx_ids {
                        *a.discards
                            .entry(Identifier::TxId(tx_id.clone()))
                            .or_default() |= Reason::TranscriptPriority;
                    }
                }
            }
            a.merge(&mut b);
            a
        })
        .unwrap();
    merged.apply_fixes(&labels);

    tracing::info!(
        "… done loading cdot JSON files in {:?} -- #genes = {}, #transcripts = {}, #transcript_ids_for_gene = {}",
        start.elapsed(),
        merged.gene_id_to_gene.len().separate_with_underscores(),
        merged.transcripts.len().separate_with_underscores(),
        merged.hgnc_id_to_transcript_ids.len().separate_with_underscores()
    );

    Ok(merged)
}

/// Main entry point for `db create txs` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
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

    // Open seqrepo …
    let seqrepo = open_seqrepo(&args.path_seqrepo_instance)?;

    // … then load cdot files …
    let mut tx_data = load_cdot_files(args)?;
    report(ReportEntry::Log(json!({
        "source": "cdot",
        "total_genes": tx_data.gene_id_to_gene.len(),
        "total_transcripts": tx_data.transcripts.len(),
        "total_hgnc_ids": tx_data.hgnc_id_to_transcript_ids.len()
    })))?;

    // … then remove information for certain genes …
    let selected_hgnc_ids = args.gene_symbols.as_ref().map(|symbols| {
        let r = tx_data.symbols_to_id(symbols);
        tracing::info!("Will limit to {:?}", r);
        r
    });

    tx_data.filter_genes()?;
    tx_data.filter_transcripts(args.max_txs, &selected_hgnc_ids)?;
    tx_data.filter_transcripts_with_sequence(&seqrepo)?;

    report(ReportEntry::Log(json!({
        "source": "cdot_filtered",
        "total_genes": tx_data.gene_id_to_gene.len(),
        "total_transcripts": tx_data.transcripts.len(),
        "total_hgnc_ids": tx_data.hgnc_id_to_transcript_ids.len()
    })))?;

    // Count number of chrMt, MANE Select and MANE Plus Clinical transcripts.
    let (n_mt, n_mane_select, n_mane_plus_clinical) = tx_data.gather_transcript_stats()?;

    report(ReportEntry::Log(json!({
        "source": "cdot_filtered",
        "n_mt": n_mt,
        "n_mane_select": n_mane_select,
        "n_mane_plus_clinical": n_mane_plus_clinical
    })))?;

    // … and finally build protobuf file.
    let tx_db =
        tx_data.build_protobuf(&seqrepo, common.verbose.is_silent(), args.genome_release)?;

    for (id, reason) in tx_data.discards {
        report(ReportEntry::Discard(Discard {
            source: "protobuf".into(),
            reason,
            id,
            gene_name: None,
            tags: None,
        }))?;
    }

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

    trace_rss_now();

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
    trace_rss_now();

    Ok(())
}

#[cfg(test)]
pub mod test {
    use std::path::{Path, PathBuf};

    use clap_verbosity_flag::Verbosity;
    use hgvs::sequences::{translate_cds, TranslationTable};
    use seqrepo::{AliasOrSeqId, Interface};
    use temp_testdir::TempDir;

    use crate::common::{Args as CommonArgs, GenomeRelease};
    use crate::db::create::{open_seqrepo, TranscriptLoader};
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
            .collect::<Vec<_>>());

        tx_data.filter_transcripts(None, &None)?;
        insta::assert_yaml_snapshot!(tx_data
            .hgnc_id_to_transcript_ids
            .get(&1100)
            .unwrap()
            .iter()
            .map(|s| s.as_str())
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
