//! Consequence prediction for structural variants.

use std::{collections::HashMap, sync::Arc};

use hgvs::{data::interface::Provider, static_data::Assembly};

use crate::{
    annotate::seqvars::provider::{MehariProvider, TxIntervalTrees},
    db::create::txs::data::{Strand, Transcript, TxSeqDatabase},
};

/// Enumeration for effect on transcript.
#[derive(
    serde::Serialize, serde::Deserialize, PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Copy,
)]
#[serde(rename_all = "snake_case")]
pub enum TranscriptEffect {
    /// Affects the full transcript.
    TranscriptVariant,
    /// An exon is affected by the SV.
    ExonVariant,
    /// The splice region is affected by the SV.
    SpliceRegionVariant,
    /// The intron is affected by the SV.
    IntronVariant,
    /// The upstream region of the transcript is affected.
    UpstreamVariant,
    /// The downstream region of the transcript is affected.
    DownstreamVariant,
    /// Only intergenic regions is affected,
    IntergenicVariant,
}

impl TranscriptEffect {
    /// Return vector with all transcript effects.
    pub fn vec_all() -> Vec<TranscriptEffect> {
        vec![
            TranscriptEffect::TranscriptVariant,
            TranscriptEffect::ExonVariant,
            TranscriptEffect::SpliceRegionVariant,
            TranscriptEffect::IntronVariant,
            TranscriptEffect::UpstreamVariant,
            TranscriptEffect::DownstreamVariant,
            TranscriptEffect::IntergenicVariant,
        ]
    }
}

/// Helpful types for using structural variants.
pub mod interface {
    /// Structural Variant type.
    #[derive(
        serde::Serialize, serde::Deserialize, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash,
    )]
    pub enum StrucVarType {
        #[serde(rename = "DEL")]
        Del,
        #[serde(rename = "DUP")]
        Dup,
        #[serde(rename = "INS")]
        Ins,
        #[serde(rename = "INV")]
        Inv,
        #[serde(rename = "BND")]
        Bnd,
    }

    /// Strand orientation of a structural variant.
    #[derive(
        serde::Serialize,
        serde::Deserialize,
        Debug,
        Clone,
        Default,
        PartialEq,
        Eq,
        PartialOrd,
        Ord,
        Hash,
    )]
    pub enum StrandOrientation {
        #[serde(rename = "3to3")]
        ThreeToThree,
        #[serde(rename = "5to5")]
        FiveToFive,
        #[serde(rename = "3to5")]
        ThreeToFive,
        #[serde(rename = "5to3")]
        FiveToThree,
        #[serde(rename = "NtoN")]
        #[default]
        NotApplicable,
    }

    /// Trait for generic description of structural variant.
    ///
    /// The values of ``stop()`` and ``chrom2()`` are unused then they are set to the
    /// corresponding other value: ``start()`` and ``chrom()``.
    pub trait StrucVar {
        /// Chromosome name.
        fn chrom(&self) -> String;
        /// The second involved chromosome, or the same value as `chrom()`.
        fn chrom2(&self) -> String;

        /// 1-based start position of the variant (or position on first chromosome for break-ends)
        fn start(&self) -> i32;
        /// 1-based end position of the variant.
        fn stop(&self) -> i32;

        /// Type of the structural variant
        fn sv_type(&self) -> StrucVarType;
        /// The strand orientation of the structural variant, if applicable.
        fn strand_orientation(&self) -> StrandOrientation;
    }
}

/// Ad-hoc data structure for `tx_regions`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct TxRegion {
    // 0-based begin position
    begin: i32,
    // 0-based end position
    end: i32,
    // "arbitrary" number
    no: usize,
    // effect of the transcript (encodes region type)
    effect: TranscriptEffect,
}

/// Explanation of transcript effect per individual gene.
#[derive(serde::Serialize, serde::Deserialize, Debug, Clone)]
pub struct GeneTranscriptEffects {
    /// HGNC identifier
    hgnc_id: String,
    /// Transcript effects for the gene.
    transcript_effects: Vec<TranscriptEffect>,
}

/// Length of the upstream/downstream region.
static X_STREAM: i32 = 5000;

/// Return list of half-open intervals for a given transcript.
fn tx_regions(tx: &Transcript) -> Vec<TxRegion> {
    assert_eq!(
        tx.genome_alignments.len(),
        1,
        "only one alignment supported"
    );
    if tx.genome_alignments[0].exons.is_empty() {
        // no exons? skip!
        return Vec::new();
    }

    let mut result = Vec::new();
    let mut tx_start = None;
    let mut tx_end = None;
    let genome_alignment = &tx.genome_alignments[0];

    // Loop over all exons to determine leftmost and rightmost genome position.
    for exon_alignment in &genome_alignment.exons {
        if let Some(value) = tx_start {
            if exon_alignment.alt_start_i < value {
                tx_start = Some(exon_alignment.alt_start_i);
            }
        } else {
            tx_start = Some(exon_alignment.alt_start_i);
        }
        if let Some(value) = tx_end {
            if exon_alignment.alt_end_i > value {
                tx_end = Some(exon_alignment.alt_end_i);
            }
        } else {
            tx_end = Some(exon_alignment.alt_end_i);
        }
    }

    let tx_start = tx_start.expect("must have been set");
    let tx_end = tx_end.expect("must have been set");

    // Perform the actual region extraction.
    let mut prev_alt_end_i = 0;
    for (no, exon_alignment) in genome_alignment.exons.iter().enumerate() {
        if exon_alignment.alt_start_i == tx_start {
            // is first, register upstream/downstream
            result.push(TxRegion {
                begin: exon_alignment.alt_start_i - X_STREAM,
                end: exon_alignment.alt_start_i - 1,
                no,
                effect: if genome_alignment.strand == Strand::Plus as i32 {
                    TranscriptEffect::UpstreamVariant
                } else {
                    TranscriptEffect::DownstreamVariant
                },
            });
        } else {
            // is not first, register splice region on the left boundary
            result.push(TxRegion {
                begin: (exon_alignment.alt_start_i - 1) - 8,
                end: (exon_alignment.alt_start_i - 1) + 3,
                no,
                effect: TranscriptEffect::SpliceRegionVariant,
            })
        }

        if exon_alignment.alt_end_i == tx_end {
            // is last, register upstream/downstream
            result.push(TxRegion {
                begin: exon_alignment.alt_end_i,
                end: exon_alignment.alt_end_i + X_STREAM,
                no,
                effect: if genome_alignment.strand == Strand::Plus as i32 {
                    TranscriptEffect::DownstreamVariant
                } else {
                    TranscriptEffect::UpstreamVariant
                },
            });
        } else {
            // is not last, register splice region on the right boundary
            result.push(TxRegion {
                begin: exon_alignment.alt_end_i - 3,
                end: exon_alignment.alt_end_i + 8,
                no,
                effect: TranscriptEffect::SpliceRegionVariant,
            })
        }

        // register the exon
        result.push(TxRegion {
            begin: exon_alignment.alt_start_i - 1,
            end: exon_alignment.alt_end_i,
            no,
            effect: TranscriptEffect::ExonVariant,
        });

        if exon_alignment.alt_start_i != tx_start {
            // is not first exon, register intron "right" of it
            result.push(TxRegion {
                begin: prev_alt_end_i,
                end: exon_alignment.alt_start_i - 1,
                no,
                effect: TranscriptEffect::IntronVariant,
            });
        }

        // store end of prev exon for next intron's start
        prev_alt_end_i = exon_alignment.alt_end_i;
    }

    result
}

/// Return the transcript region / effect for the given breakpoint.
fn gene_tx_effects_for_bp(tx: &Transcript, pos: i32) -> Vec<TranscriptEffect> {
    // Obtain list of regions for transcript.
    let regions = tx_regions(tx);

    // Determine how this relates to the breakpoint.
    let pos = pos - 1; // 1-based to 0-based
    let mut result = regions
        .iter()
        .filter(|r| r.begin <= pos && pos < r.end)
        .map(|r| r.effect)
        .collect::<Vec<_>>();
    if result.is_empty() {
        result.push(TranscriptEffect::IntergenicVariant);
    } else {
        result.sort();
        result.dedup();
    }
    result
}

/// Return the transcript region / effect for the given range.
fn gene_tx_effect_for_range(tx: &Transcript, start: i32, stop: i32) -> Vec<TranscriptEffect> {
    // Obtain list of regions for transcript.
    let regions = tx_regions(tx);

    // Determine how this relates to the left and right breakpoints.
    let pos = start - 1; // 1-based to 0-based
    let mut result = regions
        .iter()
        .filter(|region| pos <= region.end && region.begin <= stop)
        .map(|region| region.effect)
        .collect::<Vec<_>>();

    // Remove any duplicates.
    result.sort();
    result.dedup();

    // If we have both upstream and downstream then the full transcript is affected.
    if result.contains(&TranscriptEffect::UpstreamVariant)
        && result.contains(&TranscriptEffect::DownstreamVariant)
    {
        result.push(TranscriptEffect::TranscriptVariant);
    }

    result
}

/// Helper that computes effects on transcripts for a single breakend, e.g., one side of BND or INS.
fn compute_tx_effects_for_breakpoint(
    sv: &impl interface::StrucVar,
    mehari_tx_db: &TxSeqDatabase,
    mehari_tx_idx: &TxIntervalTrees,
    chrom_to_acc: &HashMap<String, String>,
) -> Vec<GeneTranscriptEffects> {
    // Shortcut to the `TranscriptDb`.
    let tx_db = mehari_tx_db
        .tx_db
        .as_ref()
        .expect("transcripts must be present");
    // Compute canonical chromosome name and map to accession.
    let chrom = chrom_to_acc.get(&annonars::common::cli::canonicalize(&sv.chrom()));
    if chrom.is_none() {
        return Default::default();
    }
    let chrom = chrom.expect("chromosome must be known at this point");
    // Create range to query the interval trees for.
    let query = (sv.start() - X_STREAM)..(sv.start() + X_STREAM);

    if let Some(idx) = mehari_tx_idx.contig_to_idx.get(chrom) {
        let mut effects_by_gene: HashMap<_, Vec<_>> = HashMap::new();

        // Collect all transcripts that overlap the INS and compute the effect of the INS on
        // the transcript.
        let tree = &mehari_tx_idx.trees[*idx];
        for it in tree.find(query) {
            let tx = &tx_db.transcripts[*it.data() as usize];
            let hgnc_id = format!("HGNC:{}", &tx.gene_id);
            effects_by_gene
                .entry(hgnc_id)
                .or_default()
                .extend(gene_tx_effects_for_bp(tx, sv.start()));
        }

        // Deduplicate results.
        effects_by_gene.iter_mut().for_each(|(_, v)| {
            v.sort();
            v.dedup()
        });

        // Convert the results into the final format.
        effects_by_gene
            .into_iter()
            .map(|(hgnc_id, transcript_effects)| GeneTranscriptEffects {
                hgnc_id,
                transcript_effects,
            })
            .collect()
    } else {
        // We do not have any transcripts for this chromosome.
        Default::default()
    }
}

/// Compute effect for linear SVs.
fn compute_tx_effects_for_linear(
    sv: &impl interface::StrucVar,
    mehari_tx_db: &TxSeqDatabase,
    mehari_tx_idx: &TxIntervalTrees,
    chrom_to_acc: &HashMap<String, String>,
) -> Vec<GeneTranscriptEffects> {
    // Shortcut to the `TranscriptDb`.
    let tx_db = mehari_tx_db
        .tx_db
        .as_ref()
        .expect("transcripts must be present");
    // Compute canonical chromosome name and map to accession.
    let chrom = chrom_to_acc.get(&annonars::common::cli::canonicalize(&sv.chrom()));
    if chrom.is_none() {
        return Default::default();
    }
    let chrom = chrom.expect("chromosome must be known at this point");
    // Create range to query the interval trees for.
    let query = (sv.start() - X_STREAM)..(sv.stop() + X_STREAM);

    if let Some(idx) = mehari_tx_idx.contig_to_idx.get(chrom) {
        let mut effects_by_gene: HashMap<_, Vec<_>> = HashMap::new();

        // Collect all transcripts that overlap the SV and compute the effect of the SV on
        // the transcript.
        let tree = &mehari_tx_idx.trees[*idx];
        for it in tree.find(query) {
            let tx = &tx_db.transcripts[*it.data() as usize];
            let hgnc_id = format!("HGNC:{}", &tx.gene_id);
            effects_by_gene
                .entry(hgnc_id)
                .or_default()
                .extend(gene_tx_effect_for_range(tx, sv.start(), sv.stop()));
        }

        // Deduplicate results.
        effects_by_gene.iter_mut().for_each(|(_, v)| {
            v.sort();
            v.dedup()
        });

        // Convert the results into the final format.
        effects_by_gene
            .into_iter()
            .map(|(hgnc_id, transcript_effects)| GeneTranscriptEffects {
                hgnc_id,
                transcript_effects,
            })
            .collect()
    } else {
        // We do not have any transcripts for this chromosome.
        Default::default()
    }
}

/// Wrap mapper, provider, and map for consequence prediction.
#[derive(derivative::Derivative)]
#[derivative(Debug)]
pub struct ConsequencePredictor {
    /// The internal transcript provider for locating transcripts.
    #[derivative(Debug = "ignore")]
    provider: Arc<MehariProvider>,
    /// Mapping from chromosome name to accession.
    #[derivative(Debug = "ignore")]
    chrom_to_acc: HashMap<String, String>,
}

impl ConsequencePredictor {
    pub fn new(provider: Arc<MehariProvider>, assembly: Assembly) -> Self {
        let acc_to_chrom = provider.get_assembly_map(assembly);
        let mut chrom_to_acc = HashMap::new();
        for (acc, chrom) in &acc_to_chrom {
            let chrom = if chrom.starts_with("chr") {
                chrom.strip_prefix("chr").unwrap()
            } else {
                chrom
            };
            chrom_to_acc.insert(chrom.to_string(), acc.clone());
            chrom_to_acc.insert(format!("chr{}", chrom), acc.clone());
        }

        ConsequencePredictor {
            provider,
            chrom_to_acc,
        }
    }

    /// Compute effect(s) of `sv` on transcript of genes.
    pub fn compute_tx_effects(
        &self,
        sv: &impl interface::StrucVar,
        // mehari_tx_db: &TxSeqDatabase,
        // mehari_tx_idx: &TxIntervalTrees,
        // chrom_to_acc: &HashMap<String, String>,
    ) -> Vec<GeneTranscriptEffects> {
        match sv.sv_type() {
            interface::StrucVarType::Ins | interface::StrucVarType::Bnd => {
                compute_tx_effects_for_breakpoint(
                    sv,
                    &self.provider.tx_seq_db,
                    &self.provider.tx_trees,
                    &self.chrom_to_acc,
                )
            }
            interface::StrucVarType::Del
            | interface::StrucVarType::Dup
            | interface::StrucVarType::Inv => compute_tx_effects_for_linear(
                sv,
                &self.provider.tx_seq_db,
                &self.provider.tx_trees,
                &self.chrom_to_acc,
            ),
        }
    }

    /// Return data version string (if set).
    pub fn data_version(&self) -> Option<String> {
        self.provider.as_ref().tx_seq_db.version.clone()
    }
}
