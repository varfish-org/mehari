//! Implementation of `hgvs` Provider interface based on protobuf.

use crate::annotate::seqvars::reference::{
    InMemoryFastaAccess, ReferenceReader, UnbufferedIndexedFastaAccess,
};
use crate::common::contig::ContigManager;
use crate::db::TranscriptDatabase;
use crate::db::transcripts::create::models::Reason;
use crate::{
    annotate::seqvars::consequence::logic::ALT_ALN_METHOD,
    pbs::txs::{GeneToTxId, Strand, Transcript, TranscriptTag, TxSeqDatabase},
};
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use enumflags2::BitFlags;
use hgvs::{
    data::error::Error,
    data::{
        cdot::json::NCBI_ALN_METHOD,
        interface::{
            Provider as ProviderInterface, TranslationException, TxExonsRecord, TxForRegionRecord,
            TxIdentityInfo, TxInfoRecord, TxMappingOptionsRecord,
        },
    },
    sequences::{TranslationTable, seq_md5},
};
use indexmap::IndexMap;
use itertools::Itertools;
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::Path;
use std::sync::Arc;

type IntervalTree = ArrayBackedIntervalTree<i32, u32>;

pub struct TxIntervalTrees {
    /// Mapping from contig accession to index in `trees`.
    pub contig_to_idx: FxHashMap<String, usize>,
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

    fn build_indices(db: &TxSeqDatabase) -> (FxHashMap<String, usize>, Vec<IntervalTree>) {
        let assembly_name = db.assembly();
        let mut contig_to_idx = FxHashMap::default();
        let mut trees: Vec<IntervalTree> = Vec::new();

        let tx_db = db.tx_db.as_ref().expect("no tx_db?");
        let mut txs_counted = 0;

        for (tx_id, tx) in tx_db.transcripts.iter().enumerate() {
            for genome_alignment in &tx.genome_alignments {
                // Only index alignments matching the current assembly/build
                if genome_alignment
                    .genome_build
                    .eq_ignore_ascii_case(&assembly_name)
                {
                    let contig = &genome_alignment.contig;

                    let contig_idx = *contig_to_idx.entry(contig.clone()).or_insert_with(|| {
                        let idx = trees.len();
                        trees.push(IntervalTree::new());
                        idx
                    });

                    let mut start = i32::MAX;
                    let mut stop = i32::MIN;
                    for exon in &genome_alignment.exons {
                        start = std::cmp::min(start, exon.alt_start_i);
                        stop = std::cmp::max(stop, exon.alt_end_i);
                    }

                    if start <= stop {
                        trees[contig_idx].insert(start..stop, tx_id as u32);
                    }
                }
            }
            txs_counted += 1;
        }

        tracing::debug!(
            "Indexed {} transcripts across {} contigs for assembly {}",
            txs_counted,
            trees.len(),
            assembly_name
        );

        trees.iter_mut().for_each(|t| t.index());

        (contig_to_idx, trees)
    }

    pub fn get_tx_for_region(
        &self,
        tx_seq_db: &TxSeqDatabase,
        alt_ac: &str,
        _alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, Error> {
        let contig_idx = match self.contig_to_idx.get(alt_ac) {
            Some(idx) => *idx,
            None => return Ok(Vec::new()),
        };
        let query = start_i..end_i;
        let tx_idxs = self.trees[contig_idx].find(query);

        Ok(tx_idxs
            .iter()
            .map(|entry| {
                let tx = &tx_seq_db.tx_db.as_ref().expect("no tx_db?").transcripts
                    [*entry.data() as usize];
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
}

/// Extension trait to easily check the filter status of a transcript.
pub trait PbsTranscriptExt {
    /// Returns true if the transcript has absolutely no filter flags (hard or soft).
    fn is_clean(&self) -> bool;

    /// Returns true if the transcript is incomplete (e.g. missing stop codon).
    fn is_incomplete_3p(&self) -> bool;

    /// Returns true if the CDS is incomplete at its 5' end, i.e. lacks its start codon.
    fn is_incomplete_5p(&self) -> bool;

    /// Returns the transcript length: the last transcript position that an exon covers.
    ///
    /// The sum of the genomic exon lengths would count genome-only bases (CIGAR `I`) and miss
    /// transcript-only bases (CIGAR `D`). The stored transcript sequence can end with `A` bases
    /// that `db create` appends. Bases after the last exon, such as an unaligned poly-A tail,
    /// do not count.
    fn tx_len(&self) -> i32;

    /// Computes the available CDS length. If the stop codon is missing,
    /// it uses the provided transcript length to bound the remaining sequence.
    fn available_cds_len(&self, tx_len: i32) -> Option<i32>;

    /// The number of bases after the last exon that complete the stop codon. `db create`
    /// adds them as `A` bases, as the poly-A tail does. No total counts them.
    fn stop_codon_padding(&self, tx_len: i32) -> i32;

    /// `seq` without the `A` bases that `db create` appended to the stored sequence to
    /// complete the stop codon. `ref_len` is the length of the stored sequence. If it runs past
    /// the stop codon, e.g. into a poly-A tail, `db create` appended nothing.
    ///
    /// A RefSeq sequence can have an unaligned 3' tail that ends inside the completed stop codon.
    /// `db create` pads it to the stop codon end as well, so this drops the tail bases too. No
    /// such transcript is known.
    fn without_stop_codon_padding(&self, seq: &str, ref_len: usize, tx_len: i32) -> String;
}

impl PbsTranscriptExt for Transcript {
    fn is_clean(&self) -> bool {
        self.filter_reason.unwrap_or(0) == 0
    }

    fn is_incomplete_3p(&self) -> bool {
        let flags = BitFlags::<Reason>::from_bits_truncate(self.filter_reason.unwrap_or(0));
        flags.intersects(Reason::MissingStopCodon | Reason::ThreePrimeEndTruncated)
    }

    fn is_incomplete_5p(&self) -> bool {
        let flags = BitFlags::<Reason>::from_bits_truncate(self.filter_reason.unwrap_or(0));
        flags.contains(Reason::FivePrimeEndTruncated)
    }

    fn tx_len(&self) -> i32 {
        // `alt_cds_end_i` is the 1-based transcript position of the exon's last base.
        self.genome_alignments
            .iter()
            .flat_map(|alignment| &alignment.exons)
            .filter_map(|exon| exon.alt_cds_end_i)
            .max()
            .unwrap_or_default()
    }

    fn available_cds_len(&self, tx_len: i32) -> Option<i32> {
        self.start_codon
            .map(|start| self.stop_codon.unwrap_or(tx_len) - start)
    }

    fn stop_codon_padding(&self, tx_len: i32) -> i32 {
        self.stop_codon.map_or(0, |stop| (stop - tx_len).max(0))
    }

    fn without_stop_codon_padding(&self, seq: &str, ref_len: usize, tx_len: i32) -> String {
        let appended =
            if self.stop_codon.and_then(|stop| usize::try_from(stop).ok()) == Some(ref_len) {
                usize::try_from(self.stop_codon_padding(tx_len)).unwrap_or_default()
            } else {
                0
            };
        seq.get(..seq.len().saturating_sub(appended))
            .unwrap_or(seq)
            .to_string()
    }
}

/// Provider based on the protobuf `TxSeqDatabase`.
pub struct Provider {
    /// Database of transcripts and sequences as deserialized from protobuf.
    pub tx_seq_db: TxSeqDatabase,

    /// Interval trees for the tanscripts.
    pub tx_trees: TxIntervalTrees,

    /// The contig name manager.
    pub contig_manager: Arc<ContigManager>,

    /// Mapping from HGNC ID to the indices in `TxSeqDatabase::tx_db::gene_to_tx`. Two genes can
    /// have one HGNC ID, e.g. a RefSeq and an Ensembl gene.
    gene_map: FxHashMap<String, Vec<u32>>,

    /// Mapping from transcript accession to index in `TxSeqDatabase::tx_db::transcripts`.
    tx_map: FxHashMap<String, u32>,

    /// Index in `TxSeqDatabase::tx_db::gene_to_tx` of the gene of each transcript, by index in
    /// `TxSeqDatabase::tx_db::transcripts`.
    tx_genes: Vec<Option<u32>>,

    /// Mapping from sequence accession to index in `TxSeqDatabase::seq_db::seqs`.
    seq_map: FxHashMap<String, u32>,

    /// Map from contig accession to common name (e.g., "NC_000001.11" -> "1")
    assembly_map: IndexMap<String, String>,

    /// Maps any known contig alias (e.g. "chr1", "NC_000001.11") to the exact
    /// string used in the database (e.g. "1")
    contig_alias_map: FxHashMap<String, String>,

    /// The contigs of the database that are mitochondrial, see `translation_table`.
    mitochondrial_contigs: FxHashSet<String>,

    /// for reading parts of reference sequences
    reference_reader: Option<ReferenceReaderImpl>,

    /// The data version.
    data_version: String,

    /// The schema version.
    schema_version: String,
}

enum ReferenceReaderImpl {
    InMemory(InMemoryFastaAccess),
    Unbuffered(UnbufferedIndexedFastaAccess),
}

impl ReferenceReader for ReferenceReaderImpl {
    fn get(
        &self,
        ac: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> anyhow::Result<Option<Vec<u8>>> {
        match self {
            ReferenceReaderImpl::InMemory(v) => v.get(ac, start, end),
            ReferenceReaderImpl::Unbuffered(v) => v.get(ac, start, end),
        }
    }
}

impl Provider {
    /// Create a new `MehariProvider` from a `TxSeqDatabase`.
    ///
    /// # Arguments
    ///
    /// * `tx_seq_db` - The `TxSeqDatabase` to use.
    pub fn new(
        tx_seq_db: TxSeqDatabase,
        reference: Option<impl AsRef<Path>>,
        in_memory_reference: bool,
    ) -> Self {
        let assembly_name = tx_seq_db.assembly().clone();

        // ContigManager is still useful for name resolution (e.g. NC -> chr name)
        // If it's a custom assembly, it might just return the accession.
        let contig_manager = Arc::new(ContigManager::new(&assembly_name));

        // Build the assembly map from alignments in the DB
        let mut assembly_map = IndexMap::new();
        if let Some(tx_db) = &tx_seq_db.tx_db {
            for tx in &tx_db.transcripts {
                for aln in &tx.genome_alignments {
                    if aln.genome_build.eq_ignore_ascii_case(&assembly_name)
                        && !assembly_map.contains_key(&aln.contig)
                    {
                        let common_name = contig_manager
                            .get_primary_name(&aln.contig)
                            .cloned()
                            .unwrap_or_else(|| aln.contig.clone());

                        assembly_map.insert(aln.contig.clone(), common_name);
                    }
                }
            }
        }

        let tx_trees = TxIntervalTrees::new(&tx_seq_db);

        let mut contig_alias_map = FxHashMap::default();
        if let Some(tx_db) = &tx_seq_db.tx_db {
            for tx in &tx_db.transcripts {
                for aln in &tx.genome_alignments {
                    let db_contig = aln.contig.clone();
                    contig_alias_map.insert(db_contig.clone(), db_contig.clone());

                    if let Some(primary) = contig_manager.get_primary_name(&db_contig) {
                        contig_alias_map.insert(primary.to_string(), db_contig.clone());
                    }
                    if let Some(acc) = contig_manager.get_accession(&db_contig) {
                        contig_alias_map.insert(acc.to_string(), db_contig.clone());
                    }
                    if let Some(info) = contig_manager.get_contig_info(&db_contig) {
                        contig_alias_map.insert(info.name_with_chr.clone(), db_contig.clone());
                        contig_alias_map.insert(info.name_without_chr.clone(), db_contig.clone());
                    }
                }
            }
        }
        // The keys include the contigs of all alignments.
        let mitochondrial_contigs = (contig_alias_map.keys())
            .filter(|contig| contig_manager.is_mitochondrial_alias(contig))
            .cloned()
            .collect();

        let mut gene_map: FxHashMap<String, Vec<u32>> = FxHashMap::default();
        for (idx, entry) in tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .gene_to_tx
            .iter()
            // .filter(|gene| !gene.filtered.unwrap_or(false))
            .enumerate()
        {
            if let Some(hgnc_id) = &entry.hgnc_id {
                gene_map
                    .entry(hgnc_id.clone())
                    .or_default()
                    .push(idx as u32);
            }
        }
        let tx_map = FxHashMap::from_iter(
            tx_seq_db
                .tx_db
                .as_ref()
                .expect("no tx_db?")
                .transcripts
                .iter()
                // .filter(|tx| !tx.filtered.unwrap_or(false))
                .enumerate()
                .map(|(idx, tx)| (tx.id.clone(), idx as u32)),
        );
        let mut tx_genes = Vec::new();
        if let Some(tx_db) = &tx_seq_db.tx_db {
            tx_genes.resize(tx_db.transcripts.len(), None);
            for (gene_idx, gene) in tx_db.gene_to_tx.iter().enumerate() {
                for tx_id in &gene.tx_ids {
                    if let Some(tx_gene) = tx_map
                        .get(tx_id)
                        .and_then(|tx_idx| tx_genes.get_mut(*tx_idx as usize))
                    {
                        tx_gene.get_or_insert(gene_idx as u32);
                    }
                }
            }
        }
        let seq_map = FxHashMap::from_iter(
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

        // TODO obtain or construct data_version and schema_version somehow
        // for now, these are just set to a combination of things that make this provider instance
        // unique, such that `data_version` and `schema_version` can be used as keys for
        // caching.
        let data_version = format!(
            "{}{}",
            tx_seq_db.version.as_ref().unwrap_or(&"".to_string()),
            tx_seq_db
                .source_version
                .iter()
                .map(|v| format!("{:#?}", v))
                .join(","),
        );
        let schema_version = data_version.clone();

        let reference_reader = match (reference, in_memory_reference) {
            (Some(path), false) => Some(ReferenceReaderImpl::Unbuffered(
                UnbufferedIndexedFastaAccess::from_path(path, contig_manager.clone())
                    .expect("Failed to open reference FASTA file"),
            )),
            (Some(path), true) => Some(ReferenceReaderImpl::InMemory(
                InMemoryFastaAccess::from_path(path, contig_manager.clone())
                    .expect("Failed to open reference FASTA file"),
            )),
            // Once per assembly, as predictors can share a provider.
            (None, _) => {
                tracing::warn!(
                    "No reference genome given for {}: indels in introns or at exon edges may \
                     not be 3'-shifted in their HGVS descriptions, and their consequences can be \
                     off",
                    assembly_name
                );
                None
            }
        };

        Self {
            tx_seq_db,
            tx_trees,
            contig_manager,
            contig_alias_map,
            mitochondrial_contigs,
            gene_map,
            tx_map,
            tx_genes,
            seq_map,
            assembly_map,
            reference_reader,
            data_version,
            schema_version,
        }
    }

    /// Return the assembly of the provider.
    ///
    /// # Returns
    ///
    /// The assembly of the provider.
    pub fn assembly(&self) -> String {
        self.tx_seq_db.assembly().clone()
    }

    /// Return the picked transcript IDs of all genes with the HGNC ID.
    ///
    /// # Args
    ///
    /// * `hgnc_id` - The gene HGNC ID.
    ///
    /// # Returns
    ///
    /// The picked transcript IDs, or None if no unfiltered gene has the HGNC ID.
    pub fn get_picked_transcripts(&self, hgnc_id: &str) -> Option<Vec<String>> {
        let gene_to_tx = &self.tx_seq_db.tx_db.as_ref()?.gene_to_tx;
        let genes: Vec<_> = self
            .gene_map
            .get(hgnc_id)?
            .iter()
            .filter_map(|gene_idx| gene_to_tx.get(*gene_idx as usize))
            .filter(|gene| gene.filtered != Some(true))
            .collect();
        (!genes.is_empty()).then(|| {
            genes
                .iter()
                .flat_map(|gene| gene.tx_ids.iter().cloned())
                .collect()
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
    pub fn get_tx(&self, tx_id: &str) -> Option<&Transcript> {
        self.tx_map.get(tx_id).and_then(|idx| {
            let tx = &self
                .tx_seq_db
                .tx_db
                .as_ref()
                .expect("no tx_db?")
                .transcripts[*idx as usize];
            if let Some(true) = tx.filtered {
                None
            } else {
                Some(tx)
            }
        })
    }

    /// Return the gene entry of the transcript `tx_ac`, if any.
    pub fn gene_of(&self, tx_ac: &str) -> Option<&GeneToTxId> {
        let tx_idx = *self.tx_map.get(tx_ac)?;
        let gene_idx = (*self.tx_genes.get(tx_idx as usize)?)?;
        self.tx_seq_db
            .tx_db
            .as_ref()?
            .gene_to_tx
            .get(gene_idx as usize)
    }

    /// Return the key that groups `tx` by gene: the gene ID of its gene entry, else the gene ID of
    /// `tx`. The versions of a gene share the entry. A RefSeq and an Ensembl gene stay apart, also
    /// if they have the same HGNC ID.
    pub fn gene_key<'a>(&'a self, tx: &'a Transcript) -> &'a str {
        self.gene_of(&tx.id)
            .map(|gene| gene.gene_id.as_str())
            .unwrap_or(&tx.gene_id)
    }

    pub fn reference_available(&self) -> bool {
        self.reference_reader.is_some()
    }

    /// The genetic code of `tx`. A selenoprotein without the positions of its selenocysteine
    /// codons reads every UGA as selenocysteine, except at the CDS end, see
    /// `get_tx_translation_exceptions`.
    fn translation_table(&self, tx: &Transcript) -> TranslationTable {
        let is_selenoprotein = tx.tags.contains(&(TranscriptTag::Selenoprotein as i32));
        let has_selenocysteine = tx
            .translation_exceptions
            .iter()
            .any(|exception| exception.amino_acid == "U");
        let is_mitochondrial = (tx.genome_alignments.first())
            .is_some_and(|alignment| self.mitochondrial_contigs.contains(&alignment.contig));
        if is_mitochondrial {
            TranslationTable::VertebrateMitochondrial
        } else if is_selenoprotein && !has_selenocysteine {
            // With selenocysteine in the translation exceptions, UGA reads as
            // selenocysteine only there.
            TranslationTable::Selenocysteine
        } else {
            TranslationTable::Standard
        }
    }
}

impl ProviderInterface for Provider {
    fn data_version(&self) -> &str {
        // TODO replace with proper data_version (see comment in `Provider::new`)
        &self.data_version
    }

    fn schema_version(&self) -> &str {
        // TODO replace with proper schema_version (see comment in `Provider::new`)
        &self.schema_version
    }

    fn get_assembly_map(
        &self,
        _assembly: &str,
    ) -> Result<IndexMap<String, String>, hgvs::data::error::Error> {
        // We ignore the `assembly` argument because this provider
        // is bound to a specific build during construction anyway.
        Ok(self.assembly_map.clone())
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
        let is_contig = self.contig_alias_map.contains_key(ac);
        let seq = if is_contig && self.reference_available() {
            let reader = self.reference_reader.as_ref().unwrap();
            let seq = reader
                .get(ac, begin.map(|x| x as u64), end.map(|x| x as u64))
                .map_err(|e| {
                    tracing::error!("Failed to fetch sequence for {}: {}", ac, e);
                    Error::NoSequenceRecord(format!("{} - {}", ac, e))
                })?
                .ok_or_else(|| Error::NoSequenceRecord(ac.to_string()))?;

            return String::from_utf8(seq).map_err(|_| {
                Error::NoSequenceRecord("Failed converting seq to UTF-8.".to_string())
            });
        } else {
            let seq_idx = *self
                .seq_map
                .get(ac)
                .ok_or_else(|| Error::NoSequenceRecord(ac.to_string()))?;
            let seq_idx = seq_idx as usize;
            self.tx_seq_db.seq_db.as_ref().expect("no seq_db?").seqs[seq_idx].as_bytes()
        };

        let slice = match (begin, end) {
            (Some(begin), Some(end)) => {
                let begin = std::cmp::min(begin, seq.len());
                let end = std::cmp::min(end, seq.len());
                &seq[begin..end]
            }
            (Some(begin), None) => {
                let begin = std::cmp::min(begin, seq.len());
                &seq[begin..]
            }
            (None, Some(end)) => &seq[..end],
            (None, None) => seq,
        };
        Ok(String::from_utf8_lossy(slice).to_string())
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
        let db_contig = self
            .contig_alias_map
            .get(alt_ac)
            .map(String::as_str)
            .unwrap_or(alt_ac);
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

        let hgnc = tx.hgnc_id.clone().unwrap_or_default();
        let tx_ac_str = tx_ac.to_string();
        let alt_ac_str = db_contig.to_string();
        let alt_aln_method_str = ALT_ALN_METHOD.to_string();

        for genome_alignment in &tx.genome_alignments {
            if genome_alignment.contig == db_contig {
                let alt_strand =
                    match Strand::try_from(genome_alignment.strand).expect("invalid strand") {
                        Strand::Plus => 1,
                        Strand::Minus => -1,
                        _ => unreachable!("invalid strand {}", &genome_alignment.strand),
                    };

                let mut exons = Vec::with_capacity(genome_alignment.exons.len());

                for exon in &genome_alignment.exons {
                    exons.push(TxExonsRecord {
                        hgnc: hgnc.clone(),
                        tx_ac: tx_ac_str.clone(),
                        alt_ac: alt_ac_str.clone(),
                        alt_aln_method: alt_aln_method_str.clone(),
                        alt_strand,
                        ord: exon.ord,
                        tx_start_i: exon.alt_cds_start_i.map(|val| val - 1).unwrap_or(-1),
                        tx_end_i: exon.alt_cds_end_i.unwrap_or(-1),
                        alt_start_i: exon.alt_start_i,
                        alt_end_i: exon.alt_end_i,
                        cigar: exon.cigar.clone(),
                        tx_aseq: None,
                        alt_aseq: None,
                        tx_exon_set_id: i32::MAX,
                        alt_exon_set_id: i32::MAX,
                        tx_exon_id: i32::MAX,
                        alt_exon_id: i32::MAX,
                        exon_aln_id: i32::MAX,
                    });
                }

                exons.sort_unstable_by_key(|e| e.alt_start_i);
                return Ok(exons);
            }
        }

        Err(Error::NoAlignmentFound(
            tx_ac.to_string(),
            alt_ac.to_string(),
        ))
    }

    fn get_tx_exon_coords(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        _alt_aln_method: &str,
    ) -> Result<Vec<(i32, i32)>, Error> {
        let db_contig = self
            .contig_alias_map
            .get(alt_ac)
            .map(String::as_str)
            .unwrap_or(alt_ac);
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx as usize];

        for genome_alignment in &tx.genome_alignments {
            if genome_alignment.contig == db_contig {
                let mut coords = Vec::with_capacity(genome_alignment.exons.len());
                for exon in &genome_alignment.exons {
                    let tx_start = exon.alt_cds_start_i.map(|val| val - 1).unwrap_or(-1);
                    let tx_end = exon.alt_cds_end_i.unwrap_or(-1);
                    coords.push((exon.alt_start_i, tx_start, tx_end));
                }
                coords.sort_unstable_by_key(|c| c.0);
                return Ok(coords.into_iter().map(|(_, s, e)| (s, e)).collect());
            }
        }
        Err(Error::NoAlignmentFound(
            tx_ac.to_string(),
            db_contig.to_string(),
        ))
    }

    fn get_cds_start_end(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        _alt_aln_method: &str,
    ) -> Result<(Option<i32>, Option<i32>), Error> {
        let db_contig = self
            .contig_alias_map
            .get(alt_ac)
            .map(String::as_str)
            .unwrap_or(alt_ac);
        let tx_idx = *self
            .tx_map
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        let tx = &self
            .tx_seq_db
            .tx_db
            .as_ref()
            .expect("no tx_db?")
            .transcripts[tx_idx as usize];

        for genome_alignment in &tx.genome_alignments {
            if genome_alignment.contig == db_contig {
                return Ok((tx.start_codon, tx.stop_codon));
            }
        }
        Err(Error::NoAlignmentFound(
            tx_ac.to_string(),
            db_contig.to_string(),
        ))
    }

    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<TxInfoRecord>, Error> {
        let tx_acs = if let Some(tx_acs) = self.get_picked_transcripts(gene) {
            tx_acs
        } else {
            tracing::warn!("no transcripts found for gene: {}", gene);
            return Ok(Vec::default());
        };

        tx_acs
            .iter()
            .filter_map(|tx_ac| -> Option<Result<TxInfoRecord, Error>> {
                if let Some(tx_idx) = self.tx_map.get(tx_ac) {
                    let tx_idx = *tx_idx as usize;
                    let tx = &self
                        .tx_seq_db
                        .tx_db
                        .as_ref()
                        .expect("no tx_db?")
                        .transcripts[tx_idx];

                    if let Some(genome_alignment) = tx.genome_alignments.first() {
                        Some(Ok(TxInfoRecord {
                            hgnc: tx.hgnc_id.clone().unwrap_or_default(),
                            cds_start_i: tx.start_codon,
                            cds_end_i: tx.stop_codon,
                            tx_ac: tx.id.clone(),
                            alt_ac: genome_alignment.contig.to_string(),
                            alt_aln_method: "splign".into(),
                        }))
                    } else {
                        Some(Err(Error::NoAlignmentFound(
                            tx_ac.to_string(),
                            format!("{:?}", self.assembly()),
                        )))
                    }
                } else {
                    tracing::warn!("transcript ID not found {} for gene {}", tx_ac, gene);
                    None
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
        let db_contig = self
            .contig_alias_map
            .get(alt_ac)
            .map(String::as_str)
            .unwrap_or(alt_ac);
        let contig_idx = match self.tx_trees.contig_to_idx.get(db_contig) {
            Some(idx) => *idx,
            None => return Ok(Vec::new()),
        };
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
        let hgnc = tx.hgnc_id.clone().unwrap_or_default();

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
            cds_start_i: tx.start_codon,
            cds_end_i: tx.stop_codon,
            lengths,
            hgnc,
            translation_table: self.translation_table(tx),
        })
    }

    fn get_tx_translation_exceptions(
        &self,
        tx_ac: &str,
    ) -> Result<Vec<TranslationException>, Error> {
        let tx = self
            .get_tx(tx_ac)
            .ok_or_else(|| Error::NoTranscriptFound(tx_ac.to_string()))?;
        let mut exceptions = tx
            .translation_exceptions
            .iter()
            .map(|exception| {
                let amino_acid = exception.amino_acid.parse().map_err(|_| {
                    Error::UnsupportedTranslationException(
                        tx_ac.to_string(),
                        exception.amino_acid.clone(),
                    )
                })?;
                Ok(TranslationException {
                    position: exception.position,
                    amino_acid,
                })
            })
            .collect::<Result<Vec<_>, Error>>()?;
        // A UGA at the CDS end is the stop codon, also where the code reads UGA as
        // selenocysteine, e.g., in SEPHS2.
        if self.translation_table(tx) == TranslationTable::Selenocysteine
            && let (Some(start), Some(stop)) = (tx.start_codon, tx.stop_codon)
            && let (Ok(cds_len), Ok(stop)) = (u32::try_from(stop - start), usize::try_from(stop))
            && cds_len >= 3
            && cds_len % 3 == 0
            && self
                .get_seq_part(tx_ac, Some(stop - 3), Some(stop))?
                .eq_ignore_ascii_case("TGA")
        {
            exceptions.push(TranslationException {
                position: cds_len / 3,
                amino_acid: '*',
            });
        }
        Ok(exceptions)
    }

    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        _alt_aln_method: &str,
    ) -> Result<TxInfoRecord, Error> {
        let db_contig = self
            .contig_alias_map
            .get(alt_ac)
            .map(String::as_str)
            .unwrap_or(alt_ac);
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
            if genome_alignment.contig == db_contig {
                return Ok(TxInfoRecord {
                    hgnc: tx.hgnc_id.clone().unwrap_or_default(),
                    cds_start_i: tx.start_codon,
                    cds_end_i: tx.stop_codon,
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
        let tx = self
            .get_tx(tx_ac)
            .ok_or_else(|| Error::NoTranscriptFound(tx_ac.to_string()))?;

        let current_build = self.assembly();

        let options = tx
            .genome_alignments
            .iter()
            .filter(|aln| aln.genome_build.eq_ignore_ascii_case(&current_build))
            .map(|aln| TxMappingOptionsRecord {
                tx_ac: tx_ac.to_string(),
                alt_ac: aln.contig.clone(),
                alt_aln_method: NCBI_ALN_METHOD.to_string(),
            })
            .collect::<Vec<_>>();

        if options.is_empty() {
            return Err(Error::NoAlignmentFound(tx_ac.to_string(), current_build));
        }

        Ok(options)
    }
}

/// Small transcript databases on chr1 (GRCh38) for tests.
#[cfg(test)]
pub(crate) mod test_db {
    use crate::pbs::txs::{
        ExonAlignment, GeneToTxId, GenomeAlignment, SequenceDb, SourceVersion, Strand, Transcript,
        TranscriptBiotype, TranscriptDb, TxSeqDatabase,
    };

    /// A coding transcript `id` with one exon of 100 bases on the plus strand of chr1. The exon
    /// starts at the 0-based position `start`.
    pub(crate) fn transcript(id: &str, start: i32) -> Transcript {
        Transcript {
            id: id.into(),
            gene_symbol: "GENE1".into(),
            biotype: TranscriptBiotype::Coding.into(),
            genome_alignments: vec![GenomeAlignment {
                genome_build: "grch38".into(),
                contig: "NC_000001.11".into(),
                strand: Strand::Plus.into(),
                exons: vec![ExonAlignment {
                    alt_start_i: start,
                    alt_end_i: start + 100,
                    alt_cds_start_i: Some(1),
                    alt_cds_end_i: Some(100),
                    cigar: "100M".into(),
                    ..Default::default()
                }],
                ..Default::default()
            }],
            ..Default::default()
        }
    }

    /// The gene entry of `tx` with the gene fields of `tx`.
    pub(crate) fn gene(tx: &Transcript) -> GeneToTxId {
        GeneToTxId {
            tx_ids: vec![tx.id.clone()],
            filtered: Some(false),
            gene_id: tx.gene_id.clone(),
            hgnc_id: tx.hgnc_id.clone(),
            ..Default::default()
        }
    }

    /// A database with `transcripts` and the gene entries `gene_to_tx`, but no sequences.
    pub(crate) fn tx_seq_db(
        transcripts: Vec<Transcript>,
        gene_to_tx: Vec<GeneToTxId>,
    ) -> TxSeqDatabase {
        TxSeqDatabase {
            tx_db: Some(TranscriptDb {
                transcripts,
                gene_to_tx,
            }),
            seq_db: Some(SequenceDb::default()),
            version: Some("test".into()),
            source_version: vec![SourceVersion {
                assembly: "GRCh38".into(),
                ..Default::default()
            }],
        }
    }

    /// A database with the transcript `NM_000001.1` of the gene `11` (`HGNC:1`) and its gene
    /// entry. They have the gene fields that `db create` writes for GFF3: `gene_id` `11` and
    /// `hgnc_id` `HGNC:1`.
    pub(crate) fn one_gene_db() -> TxSeqDatabase {
        let tx = Transcript {
            gene_id: "11".into(),
            hgnc_id: Some("HGNC:1".into()),
            ..transcript("NM_000001.1", 1000)
        };
        tx_seq_db(vec![tx.clone()], vec![gene(&tx)])
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn test_sync() {
        fn is_sync<T: Sync>() {}
        is_sync::<super::Provider>();
    }

    #[test]
    fn tx_len_counts_all_exon_bases() {
        use super::PbsTranscriptExt;
        use crate::pbs::txs::{ExonAlignment, GenomeAlignment, Transcript};

        // Exons of 100 and 150 bases, in 1-based inclusive transcript coordinates.
        let exon = |ord, start, end| ExonAlignment {
            ord,
            alt_cds_start_i: Some(start),
            alt_cds_end_i: Some(end),
            ..Default::default()
        };
        let tx = Transcript {
            genome_alignments: vec![GenomeAlignment {
                exons: vec![exon(0, 1, 100), exon(1, 101, 250)],
                ..Default::default()
            }],
            ..Default::default()
        };

        assert_eq!(tx.tx_len(), 250);
    }

    /// The lookup by HGNC ID uses `hgnc_id`, not `gene_id`.
    #[rstest::rstest]
    #[case::hgnc_id("HGNC:1", true)]
    #[case::gene_id("11", false)]
    fn lookup_by_hgnc_id(#[case] hgnc_id: &str, #[case] found: bool) -> Result<(), anyhow::Error> {
        use super::{Provider, test_db};
        use hgvs::data::interface::Provider as _;

        let provider = Provider::new(test_db::one_gene_db(), None::<&str>, true);

        let txs: Vec<_> = provider
            .get_tx_for_gene(hgnc_id)?
            .into_iter()
            .map(|tx| (tx.tx_ac, tx.hgnc))
            .collect();
        let expected = if found {
            vec![("NM_000001.1".to_string(), "HGNC:1".to_string())]
        } else {
            vec![]
        };
        assert_eq!(txs, expected);

        Ok(())
    }

    /// Two genes can have the same HGNC ID. The lookup by HGNC ID returns the transcripts of
    /// both.
    #[test]
    fn transcripts_of_all_genes_with_the_hgnc_id() -> Result<(), anyhow::Error> {
        use super::Provider;
        use crate::pbs::txs::{
            ExonAlignment, GeneToTxId, GenomeAlignment, SequenceDb, SourceVersion, Transcript,
            TranscriptDb, TxSeqDatabase,
        };
        use hgvs::data::interface::Provider as _;

        let tx = |id: &str, gene_id: &str| Transcript {
            id: id.into(),
            gene_symbol: "GENE1".into(),
            gene_id: gene_id.into(),
            hgnc_id: Some("HGNC:1".into()),
            genome_alignments: vec![GenomeAlignment {
                genome_build: "grch38".into(),
                contig: "NC_000001.11".into(),
                exons: vec![ExonAlignment {
                    alt_start_i: 1000,
                    alt_end_i: 1100,
                    alt_cds_start_i: Some(1),
                    alt_cds_end_i: Some(100),
                    cigar: "100M".into(),
                    ..Default::default()
                }],
                ..Default::default()
            }],
            ..Default::default()
        };
        let gene = |tx_id: &str, gene_id: &str| GeneToTxId {
            tx_ids: vec![tx_id.into()],
            filtered: Some(false),
            gene_id: gene_id.into(),
            hgnc_id: Some("HGNC:1".into()),
            ..Default::default()
        };
        let tx_seq_db = TxSeqDatabase {
            tx_db: Some(TranscriptDb {
                transcripts: vec![tx("NM_000001.1", "11"), tx("NM_000002.1", "22")],
                gene_to_tx: vec![gene("NM_000001.1", "11"), gene("NM_000002.1", "22")],
            }),
            seq_db: Some(SequenceDb::default()),
            version: Some("test".into()),
            source_version: vec![SourceVersion {
                assembly: "GRCh38".into(),
                ..Default::default()
            }],
        };
        let provider = Provider::new(tx_seq_db, None::<&str>, true);

        let tx_acs: Vec<_> = provider
            .get_tx_for_gene("HGNC:1")?
            .into_iter()
            .map(|tx| tx.tx_ac)
            .collect();
        assert_eq!(tx_acs, ["NM_000001.1", "NM_000002.1"]);

        Ok(())
    }
}
