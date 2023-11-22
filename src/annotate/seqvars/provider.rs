//! Implementation of `hgvs` Provider interface based on protobuf.

use std::collections::HashMap;

use annonars::common::cli::CANONICAL;
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use biocommons_bioutils::assemblies::{Assembly, ASSEMBLY_INFOS};
use hgvs::{
    data::error::Error,
    data::{
        cdot::json::NCBI_ALN_METHOD,
        interface::{
            Provider as ProviderInterface, TxExonsRecord, TxForRegionRecord, TxIdentityInfo,
            TxInfoRecord, TxMappingOptionsRecord,
        },
    },
    sequences::{seq_md5, TranslationTable},
};

use crate::{
    annotate::seqvars::csq::ALT_ALN_METHOD,
    pbs::txs::{GeneToTxId, Strand, Transcript, TranscriptTag, TxSeqDatabase},
};

type IntervalTree = ArrayBackedIntervalTree<i32, u32>;

pub struct TxIntervalTrees {
    /// Mapping from contig accession to index in `trees`.
    pub contig_to_idx: HashMap<String, usize>,
    /// Interval tree to index in `TxSeqDatabase::tx_db::transcripts`, for each contig.
    pub trees: Vec<IntervalTree>,
}

impl TxIntervalTrees {
    pub fn new(db: &TxSeqDatabase, assembly: Assembly) -> Self {
        let (contig_to_idx, trees) = Self::build_indices(db, assembly);
        Self {
            contig_to_idx,
            trees,
        }
    }

    fn build_indices(
        db: &TxSeqDatabase,
        assembly: Assembly,
    ) -> (HashMap<String, usize>, Vec<IntervalTree>) {
        let mut contig_to_idx = HashMap::new();
        let mut trees: Vec<IntervalTree> = Vec::new();

        let mut txs = 0;

        // Pre-create interval trees for canonical contigs.
        ASSEMBLY_INFOS[assembly].sequences.iter().for_each(|seq| {
            if CANONICAL.contains(&seq.name.as_str()) {
                let contig_idx = *contig_to_idx
                    .entry(seq.refseq_ac.clone())
                    .or_insert(trees.len());
                if contig_idx >= trees.len() {
                    trees.push(IntervalTree::new());
                }
            }
        });

        for (tx_id, tx) in db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts
            .iter()
            .enumerate()
        {
            for genome_alignment in &tx.genome_alignments {
                let contig = &genome_alignment.contig;
                if let Some(contig_idx) = contig_to_idx.get(contig) {
                    let mut start = std::i32::MAX;
                    let mut stop = std::i32::MIN;
                    for exon in &genome_alignment.exons {
                        start = std::cmp::min(start, exon.alt_start_i);
                        stop = std::cmp::max(stop, exon.alt_end_i);
                    }
                    trees[*contig_idx].insert(start..stop, tx_id as u32);
                }
            }

            txs += 1;
        }

        tracing::debug!("Loaded {} transcript", txs);
        trees.iter_mut().for_each(|t| t.index());

        (contig_to_idx, trees)
    }
}

/// Configuration for constructing the `Provider`.
#[derive(Debug, Clone, Default, derive_builder::Builder)]
#[builder(pattern = "immutable")]
pub struct Config {
    /// Whether to use transcript picking.  When enabled, only use (a)
    /// ManeSelect+ManePlusClinical, (b) ManeSelect, (c) longest transcript
    /// (the first available).
    #[builder(default = "false")]
    pub transcript_picking: bool,
}

/// Provider based on the protobuf `TxSeqDatabase`.
pub struct Provider {
    /// Database of transcripts and sequences as deserialized from protobuf.
    pub tx_seq_db: TxSeqDatabase,
    /// Interval trees for the tanscripts.
    pub tx_trees: TxIntervalTrees,
    /// Mapping from gene identifier to index in `TxSeqDatabase::tx_db::gene_to_tx`.
    gene_map: HashMap<String, u32>,
    /// Mapping from transcript accession to index in `TxSeqDatabase::tx_db::transcripts`.
    tx_map: HashMap<String, u32>,
    /// Mapping from sequence accession to index in `TxSeqDatabase::seq_db::seqs`.
    seq_map: HashMap<String, u32>,
    /// When transcript picking is enabled, contains the `GeneToTxIdx` entries
    /// for each gene; the order matches the one of `tx_seq_db.gene_to_tx`.
    picked_gene_to_tx_id: Option<Vec<GeneToTxId>>,
}

impl Provider {
    /// Create a new `MehariProvider` from a `TxSeqDatabase`.
    ///
    /// # Arguments
    ///
    /// * `tx_seq_db` - The `TxSeqDatabase` to use.
    /// * `assembly` - The assembly to use.
    pub fn new(mut tx_seq_db: TxSeqDatabase, assembly: Assembly, config: Config) -> Self {
        let tx_trees = TxIntervalTrees::new(&tx_seq_db, assembly);
        let gene_map = HashMap::from_iter(
            tx_seq_db
                .tx_db
                .as_ref()
                .expect("no tx_db?")
                .gene_to_tx
                .iter()
                .enumerate()
                .map(|(idx, entry)| (entry.gene_id.clone(), idx as u32)),
        );
        let tx_map = HashMap::from_iter(
            tx_seq_db
                .tx_db
                .as_ref()
                .expect("no tx_db?")
                .transcripts
                .iter()
                .enumerate()
                .map(|(idx, tx)| (tx.id.clone(), idx as u32)),
        );
        let seq_map = HashMap::from_iter(
            tx_seq_db
                .seq_db
                .as_ref()
                .expect("no seq_db?")
                .aliases
                .iter()
                .zip(
                    tx_seq_db
                        .seq_db
                        .as_ref()
                        .expect("no seq_db?")
                        .aliases_idx
                        .iter(),
                )
                .map(|(alias, idx)| (alias.clone(), *idx)),
        );

        // When transcript picking is enabled, restrict to ManeSelect and ManePlusClinical if
        // we have any such transcript.  Otherwise, fall back to the longest transcript.
        let picked_gene_to_tx_id = if config.transcript_picking {
            if let Some(tx_db) = tx_seq_db.tx_db.as_mut() {
                // The new gene-to-txid mapping we will build.
                let mut new_gene_to_tx = Vec::new();

                // Process each gene.
                for entry in tx_db.gene_to_tx.iter() {
                    // First, determine whether we have any MANE transcripts.
                    let mane_tx_ids = entry
                        .tx_ids
                        .iter()
                        .filter(|tx_id| {
                            tx_map
                                .get(*tx_id)
                                .map(|tx_idx| {
                                    let tx = &tx_db.transcripts[*tx_idx as usize];
                                    tx.tags.contains(&TranscriptTag::ManePlusClinical.into())
                                        || tx.tags.contains(&TranscriptTag::ManeSelect.into())
                                })
                                .unwrap_or_default()
                        })
                        .cloned()
                        .collect::<Vec<_>>();

                    // Now, construct gene-to-txid mapping entry.
                    let new_entry = if !mane_tx_ids.is_empty() {
                        // For the case that we have MANE transcripts.
                        GeneToTxId {
                            gene_id: entry.gene_id.clone(),
                            tx_ids: mane_tx_ids,
                        }
                    } else {
                        // Otherwise, determine the longest transcript's length.
                        let (_, tx_id) = entry
                            .tx_ids
                            .iter()
                            .map(|tx_id| {
                                tx_map
                                    .get(tx_id)
                                    .map(|tx_idx| {
                                        // A slight complication, we need to look at all genome alignments...
                                        let tx = &tx_db.transcripts[*tx_idx as usize];
                                        let mut max_tx_length = 0;
                                        for genome_alignment in tx.genome_alignments.iter() {
                                            // We just count length in reference so we don't have to look
                                            // into the CIGAR string.
                                            let mut tx_length = 0;
                                            for exon_alignment in genome_alignment.exons.iter() {
                                                tx_length += exon_alignment.alt_cds_end_i()
                                                    - exon_alignment.alt_cds_start_i();
                                            }
                                            if tx_length > max_tx_length {
                                                max_tx_length = tx_length;
                                            }
                                        }
                                        (max_tx_length, tx_id.clone())
                                    })
                                    .unwrap_or_default()
                            })
                            .max()
                            .unwrap_or_else(|| panic!("no length for gene {}", &entry.gene_id));

                        GeneToTxId {
                            gene_id: entry.gene_id.clone(),
                            tx_ids: vec![tx_id],
                        }
                    };

                    tracing::trace!(
                        "picked transcripts {:?} for gene {}",
                        new_entry.tx_ids,
                        new_entry.gene_id
                    );
                    new_gene_to_tx.push(new_entry);
                }

                Some(new_gene_to_tx)
            } else {
                None
            }
        } else {
            None
        };

        Self {
            tx_seq_db,
            tx_trees,
            gene_map,
            tx_map,
            seq_map,
            picked_gene_to_tx_id,
        }
    }

    /// Return whether transcript picking is enabled.
    pub fn transcript_picking(&self) -> bool {
        self.picked_gene_to_tx_id.is_some()
    }

    /// Return the picked transcript IDs for a gene.
    ///
    /// # Args
    ///
    /// * `gene_name` - The gene HGNC ID.
    ///
    /// # Returns
    ///
    /// The picked transcript IDs, or None if the gene is not found.
    pub fn get_picked_transcripts(&self, hgnc_id: &str) -> Option<Vec<String>> {
        self.gene_map.get(hgnc_id).map(|gene_idx| {
            let gene_to_tx = if let Some(picked_gene_to_tx_id) = self.picked_gene_to_tx_id.as_ref()
            {
                picked_gene_to_tx_id
            } else {
                &self.tx_seq_db.tx_db.as_ref().expect("no tx_db?").gene_to_tx
            };

            // tracing::trace!(
            //     "get_picked_transcripts({}) = {:?}",
            //     hgnc_id,
            //     &gene_to_tx[*gene_idx as usize].tx_ids
            // );
            gene_to_tx[*gene_idx as usize].tx_ids.clone()
        })
    }

    /// Return `Transcript` for the given transcript accession.
    ///
    /// # Args
    ///
    /// * `tx_id` - The transcript accession.
    ///
    /// # Returns
    ///
    /// The `Transcript` for the given accession, or None if the accession was not found.
    pub fn get_tx(&self, tx_id: &str) -> Option<Transcript> {
        self.tx_map.get(tx_id).map(|idx| {
            let result = self
                .tx_seq_db
                .tx_db
                .as_ref()
                .expect("no tx_db?")
                .transcripts[*idx as usize]
                .clone();

            // tracing::trace!("get_tx({}) = {:?}", tx_id, &result);
            result
        })
    }
}

impl ProviderInterface for Provider {
    fn data_version(&self) -> &str {
        panic!("not implemented");
    }

    fn schema_version(&self) -> &str {
        panic!("not implemented");
    }

    fn get_assembly_map(
        &self,
        assembly: biocommons_bioutils::assemblies::Assembly,
    ) -> indexmap::IndexMap<String, String> {
        indexmap::IndexMap::from_iter(
            ASSEMBLY_INFOS[assembly]
                .sequences
                .iter()
                .map(|record| (record.refseq_ac.clone(), record.name.clone())),
        )
    }

    fn get_gene_info(&self, _hgnc: &str) -> Result<hgvs::data::interface::GeneInfoRecord, Error> {
        panic!("not implemented");
    }

    fn get_pro_ac_for_tx_ac(&self, tx_ac: &str) -> Result<Option<String>, Error> {
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx_idx = tx_idx as usize;
        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx];
        Ok(tx.protein.clone())
    }

    fn get_seq_part(
        &self,
        ac: &str,
        begin: Option<usize>,
        end: Option<usize>,
    ) -> Result<String, Error> {
        let seq_idx = *self
            .seq_map
            .get(ac)
            .ok_or(Error::NoSequenceRecord(ac.to_string()))?;
        let seq_idx = seq_idx as usize;

        let seq = &self.tx_seq_db.seq_db.as_ref().expect("no seq_db?").seqs[seq_idx];
        match (begin, end) {
            (Some(begin), Some(end)) => {
                let begin = std::cmp::min(begin, seq.len());
                let end = std::cmp::min(end, seq.len());
                Ok(seq[begin..end].to_string())
            }
            (Some(begin), None) => {
                let begin = std::cmp::min(begin, seq.len());
                Ok(seq[begin..].to_string())
            }
            (None, Some(end)) => Ok(seq[..end].to_string()),
            (None, None) => Ok(seq.clone()),
        }
    }

    fn get_acs_for_protein_seq(&self, seq: &str) -> Result<Vec<String>, Error> {
        Ok(vec![format!("MD5_{}", seq_md5(seq, true)?)])
    }

    fn get_similar_transcripts(
        &self,
        _tx_ac: &str,
    ) -> Result<Vec<hgvs::data::interface::TxSimilarityRecord>, Error> {
        panic!("not implemented");
    }

    fn get_tx_exons(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        _alt_aln_method: &str,
    ) -> Result<Vec<TxExonsRecord>, Error> {
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx_idx = tx_idx as usize;

        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx];
        for genome_alignment in &tx.genome_alignments {
            if genome_alignment.contig == alt_ac {
                return Ok(genome_alignment
                    .exons
                    .iter()
                    .map(|exon| TxExonsRecord {
                        hgnc: tx.gene_id.clone(),
                        tx_ac: tx_ac.to_string(),
                        alt_ac: alt_ac.to_string(),
                        alt_aln_method: ALT_ALN_METHOD.to_string(),
                        alt_strand: match Strand::try_from(genome_alignment.strand)
                            .expect("invalid strand")
                        {
                            Strand::Plus => 1,
                            Strand::Minus => -1,
                            _ => unreachable!("invalid strand {}", &genome_alignment.strand),
                        },
                        ord: exon.ord,
                        tx_start_i: exon.alt_cds_start_i.map(|val| val - 1).unwrap_or(-1),
                        tx_end_i: exon.alt_cds_end_i.unwrap_or(-1),
                        alt_start_i: exon.alt_start_i,
                        alt_end_i: exon.alt_end_i,
                        cigar: exon.cigar.clone(),
                        tx_aseq: None,
                        alt_aseq: None,
                        tx_exon_set_id: std::i32::MAX,
                        alt_exon_set_id: std::i32::MAX,
                        tx_exon_id: std::i32::MAX,
                        alt_exon_id: std::i32::MAX,
                        exon_aln_id: std::i32::MAX,
                    })
                    .collect());
            }
        }

        Err(Error::NoAlignmentFound(
            tx_ac.to_string(),
            alt_ac.to_string(),
        ))
    }

    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<TxInfoRecord>, Error> {
        let tx_acs = self.get_picked_transcripts(gene).ok_or_else(|| {
            tracing::warn!("gene ID not found {}", gene);
            Error::NoGeneFound(gene.to_string())
        })?;

        tx_acs
            .iter()
            .map(|tx_ac| -> Result<TxInfoRecord, Error> {
                let tx_idx = *self
                    .tx_map
                    .get(tx_ac)
                    .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
                let tx_idx = tx_idx as usize;
                let tx = &self
                    .tx_seq_db
                    .tx_db
                    .as_ref()
                    .expect("no tx_db?")
                    .transcripts[tx_idx];

                if let Some(genome_alignment) = tx.genome_alignments.first() {
                    Ok(TxInfoRecord {
                        hgnc: tx.gene_id.clone(),
                        cds_start_i: genome_alignment.cds_start,
                        cds_end_i: genome_alignment.cds_end,
                        tx_ac: tx.id.clone(),
                        alt_ac: genome_alignment.contig.to_string(),
                        alt_aln_method: "splign".into(),
                    })
                } else {
                    Err(Error::NoTranscriptFound(gene.to_string()))
                }
            })
            .collect::<Result<Vec<_>, _>>()
    }

    fn get_tx_for_region(
        &self,
        alt_ac: &str,
        _alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, Error> {
        let contig_idx = *self
            .tx_trees
            .contig_to_idx
            .get(alt_ac)
            .ok_or(Error::NoTranscriptFound(alt_ac.to_string()))?;
        let query = start_i..end_i;
        let tx_idxs = self.tx_trees.trees[contig_idx].find(query);

        Ok(tx_idxs
            .iter()
            .map(|entry| {
                let tx = &self
                    .tx_seq_db
                    .tx_db
                    .as_ref()
                    .expect("no tx_db?")
                    .transcripts[*entry.data() as usize];
                assert_eq!(
                    tx.genome_alignments.len(),
                    1,
                    "Can only have one alignment in Mehari"
                );
                let alt_strand = tx.genome_alignments.first().unwrap().strand;
                TxForRegionRecord {
                    tx_ac: tx.id.clone(),
                    alt_ac: alt_ac.to_string(),
                    alt_strand: match Strand::try_from(alt_strand).expect("invalid strand") {
                        Strand::Plus => 1,
                        Strand::Minus => -1,
                        _ => unreachable!("invalid strand {}", alt_strand),
                    },
                    alt_aln_method: ALT_ALN_METHOD.to_string(),
                    start_i,
                    end_i,
                }
            })
            .collect())
    }

    fn get_tx_identity_info(&self, tx_ac: &str) -> Result<TxIdentityInfo, Error> {
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx_idx = tx_idx as usize;
        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx];
        let is_selenoprotein = tx.tags.contains(&(TranscriptTag::Selenoprotein as i32));

        let hgnc = tx.gene_id.clone();

        let mut tmp = tx
            .genome_alignments
            .first()
            .unwrap()
            .exons
            .iter()
            .map(|exon| {
                (
                    exon.ord,
                    exon.alt_cds_end_i
                        .map(|alt_cds_end_i| alt_cds_end_i + 1 - exon.alt_cds_start_i.unwrap())
                        .unwrap_or_default(),
                )
            })
            .collect::<Vec<(i32, i32)>>();
        tmp.sort();

        let lengths = tmp.into_iter().map(|(_, length)| length).collect();
        Ok(TxIdentityInfo {
            tx_ac: tx_ac.to_string(),
            alt_ac: tx_ac.to_string(), // sic(!)
            alt_aln_method: String::from("transcript"),
            cds_start_i: tx.start_codon.unwrap_or_default(),
            cds_end_i: tx.stop_codon.unwrap_or_default(),
            lengths,
            hgnc,
            translation_table: if is_selenoprotein {
                TranslationTable::Selenocysteine
            } else {
                TranslationTable::Standard
            },
        })
    }

    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        _alt_aln_method: &str,
    ) -> Result<TxInfoRecord, Error> {
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx_idx = tx_idx as usize;
        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx];

        for genome_alignment in &tx.genome_alignments {
            if genome_alignment.contig == alt_ac {
                return Ok(TxInfoRecord {
                    hgnc: tx.gene_id.clone(),
                    cds_start_i: genome_alignment.cds_start,
                    cds_end_i: genome_alignment.cds_end,
                    tx_ac: tx.id.clone(),
                    alt_ac: alt_ac.to_string(),
                    alt_aln_method: ALT_ALN_METHOD.to_string(),
                });
            }
        }

        Err(Error::NoAlignmentFound(
            tx_ac.to_string(),
            alt_ac.to_string(),
        ))
    }

    fn get_tx_mapping_options(&self, tx_ac: &str) -> Result<Vec<TxMappingOptionsRecord>, Error> {
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx_idx = tx_idx as usize;

        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx];

        let genome_alignment = tx.genome_alignments.first().unwrap();
        Ok(vec![TxMappingOptionsRecord {
            tx_ac: tx_ac.to_string(),
            alt_ac: genome_alignment.contig.clone(),
            alt_aln_method: NCBI_ALN_METHOD.to_string(),
        }])
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn test_sync() {
        fn is_sync<T: Sync>() {}
        is_sync::<super::Provider>();
    }
}
