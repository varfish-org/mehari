//! Implementation of `hgvs` Provider interface based on flatbuffers.

use std::collections::HashMap;

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use hgvs::{data::interface::Provider as ProviderInterface, static_data::ASSEMBLY_INFOS};
use linked_hash_map::LinkedHashMap;

use crate::db::create::txs::data::{TxSeqDatabase, Transcript};

type IntervalTree = ArrayBackedIntervalTree<i32, u32>;

pub struct TxIntervalTrees {
    /// Mapping from contig accession to index in `trees`.
    pub contig_to_idx: HashMap<String, usize>,
    /// Interval tree to index in `TxSeqDatabase::tx_db::transcripts`, for each contig.
    pub trees: Vec<IntervalTree>,
}

impl TxIntervalTrees {
    pub fn new(db: &TxSeqDatabase) -> Self {
        let (contig_to_idx, trees) = Self::build_indices(db);
        Self {
            contig_to_idx,
            trees,
        }
    }

    fn build_indices(db: &TxSeqDatabase) -> (HashMap<String, usize>, Vec<IntervalTree>) {
        let mut contig_to_idx = HashMap::new();
        let mut trees: Vec<IntervalTree> = Vec::new();

        let mut txs = 0;

        for (tx_id, tx) in db
            .tx_db
            .transcripts
            .iter()
            .enumerate()
        {
            for genome_alignment in &tx.genome_alignments {
                let contig = &genome_alignment.contig;
                let contig_idx = *contig_to_idx
                    .entry(contig.clone())
                    .or_insert(trees.len());
                if contig_idx >= trees.len() {
                    trees.push(IntervalTree::new());
                }
                let mut start = std::i32::MAX;
                let mut stop = std::i32::MIN;
                for exon in &genome_alignment.exons {
                    start = std::cmp::min(start, exon.alt_start_i - 1);
                    stop = std::cmp::max(stop, exon.alt_end_i);
                }
                trees[contig_idx].insert(start..stop, tx_id as u32);
            }

            txs += 1;
        }

        tracing::debug!("Loaded {} transcript", txs);
        trees.iter_mut().for_each(|t| t.index());

        (contig_to_idx, trees)
    }
}

pub struct MehariProvider {
    pub tx_db: TxSeqDatabase,
    pub tx_trees: TxIntervalTrees,
    tx_map: HashMap<String, u32>,
}

impl MehariProvider {
    pub fn new(tx_db: TxSeqDatabase) -> Self {
        let tx_trees = TxIntervalTrees::new(&tx_db);
        let tx_map = HashMap::from_iter(
            tx_db
                .tx_db
                .transcripts
                .iter()
                .enumerate()
                .map(|(idx, tx)| (tx.id.clone(), idx as u32)),
        );

        Self {
            tx_db,
            tx_trees,
            tx_map,
        }
    }

    pub fn get_tx(&self, tx_id: &str) -> Option<Transcript> {
        if let Some(idx) = self.tx_map.get(tx_id) {
            Some(
                self.tx_db
                    .tx_db
                    .transcripts[*idx as usize]
                    .clone()
            )
        } else {
            None
        }
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
    ) -> LinkedHashMap<String, String> {
        LinkedHashMap::from_iter(
            ASSEMBLY_INFOS[assembly]
                .sequences
                .iter()
                .map(|record| (record.refseq_ac.clone(), record.name.clone())),
        )
    }

    fn get_gene_info(
        &self,
        _hgnc: &str,
    ) -> Result<hgvs::data::interface::GeneInfoRecord, anyhow::Error> {
        panic!("not implemented");
    }

    fn get_pro_ac_for_tx_ac(&self, _tx_ac: &str) -> Result<Option<String>, anyhow::Error> {
        panic!("not implemented");
    }

    fn get_seq_part(
        &self,
        _ac: &str,
        _begin: Option<usize>,
        _end: Option<usize>,
    ) -> Result<String, anyhow::Error> {
        panic!("not implemented");
    }

    fn get_acs_for_protein_seq(&self, _seq: &str) -> Result<Vec<String>, anyhow::Error> {
        panic!("not implemented");
    }

    fn get_similar_transcripts(
        &self,
        _tx_ac: &str,
    ) -> Result<Vec<hgvs::data::interface::TxSimilarityRecord>, anyhow::Error> {
        panic!("not implemented");
    }

    fn get_tx_exons(
        &self,
        _tx_ac: &str,
        _alt_ac: &str,
        _alt_aln_method: &str,
    ) -> Result<Vec<hgvs::data::interface::TxExonsRecord>, anyhow::Error> {
        panic!("not implemented");
    }

    fn get_tx_for_gene(
        &self,
        _gene: &str,
    ) -> Result<Vec<hgvs::data::interface::TxInfoRecord>, anyhow::Error> {
        panic!("not implemented");
    }

    fn get_tx_for_region(
        &self,
        _alt_ac: &str,
        _alt_aln_method: &str,
        _start_i: i32,
        _end_i: i32,
    ) -> Result<Vec<hgvs::data::interface::TxForRegionRecord>, anyhow::Error> {
        panic!("not implemented");
    }

    fn get_tx_identity_info(
        &self,
        _tx_ac: &str,
    ) -> Result<hgvs::data::interface::TxIdentityInfo, anyhow::Error> {
        panic!("not implemented");
    }

    fn get_tx_info(
        &self,
        _tx_ac: &str,
        _alt_ac: &str,
        _alt_aln_method: &str,
    ) -> Result<hgvs::data::interface::TxInfoRecord, anyhow::Error> {
        panic!("not implemented");
    }

    fn get_tx_mapping_options(
        &self,
        _tx_ac: &str,
    ) -> Result<Vec<hgvs::data::interface::TxMappingOptionsRecord>, anyhow::Error> {
        panic!("not implemented");
    }
}
