//! Implementation of `hgvs` Provider interface based on protobuf.

use std::collections::HashMap;

use annonars::common::cli::CANONICAL;
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use hgvs::{
    data::error::Error,
    data::{
        cdot::json::NCBI_ALN_METHOD,
        interface::{
            Provider as ProviderInterface, TxExonsRecord, TxForRegionRecord, TxIdentityInfo,
            TxInfoRecord, TxMappingOptionsRecord,
        },
    },
    sequences::seq_md5,
    static_data::{Assembly, ASSEMBLY_INFOS},
};

use crate::{
    annotate::seqvars::csq::ALT_ALN_METHOD,
    db::create::txs::data::{Strand, Transcript, TxSeqDatabase},
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

pub struct MehariProvider {
    pub tx_seq_db: TxSeqDatabase,
    pub tx_trees: TxIntervalTrees,
    tx_map: HashMap<String, u32>,
    seq_map: HashMap<String, u32>,
}

impl MehariProvider {
    pub fn new(tx_seq_db: TxSeqDatabase, assembly: Assembly) -> Self {
        let tx_trees = TxIntervalTrees::new(&tx_seq_db, assembly);
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

        Self {
            tx_seq_db,
            tx_trees,
            tx_map,
            seq_map,
        }
    }

    pub fn get_tx(&self, tx_id: &str) -> Option<Transcript> {
        self.tx_map.get(tx_id).map(|idx| {
            self.tx_seq_db
                .tx_db
                .as_ref()
                .expect("no tx_db?")
                .transcripts[*idx as usize]
                .clone()
        })
    }
}

impl ProviderInterface for MehariProvider {
    fn data_version(&self) -> &str {
        panic!("not implemented");
    }

    fn schema_version(&self) -> &str {
        panic!("not implemented");
    }

    fn get_assembly_map(
        &self,
        assembly: hgvs::static_data::Assembly,
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
                        hgnc: tx.gene_name.clone(),
                        tx_ac: tx_ac.to_string(),
                        alt_ac: alt_ac.to_string(),
                        alt_aln_method: ALT_ALN_METHOD.to_string(),
                        alt_strand: match Strand::try_from(genome_alignment.strand)
                            .expect("invalid strand")
                        {
                            Strand::Plus => 1,
                            Strand::Minus => -1,
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

    fn get_tx_for_gene(&self, _gene: &str) -> Result<Vec<TxInfoRecord>, Error> {
        panic!("not implemented");
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

        let hgnc = tx.gene_name.clone();

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
                    hgnc: tx.gene_name.clone(),
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
        is_sync::<super::MehariProvider>();
    }
}
