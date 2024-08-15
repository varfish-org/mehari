//! Compute molecular consequence of variants.
use std::{collections::HashMap, sync::Arc};

use crate::pbs::txs::{Strand, TranscriptBiotype, TranscriptTag};
use biocommons_bioutils::assemblies::Assembly;
use enumflags2::BitFlags;
use hgvs::parser::NoRef;
use hgvs::{
    data::interface::{Provider, TxForRegionRecord},
    mapper::{assembly, Error},
    parser::{
        Accession, CdsFrom, GenomeInterval, GenomeLocEdit, HgvsVariant, Mu, NaEdit, ProtLocEdit,
    },
};
use itertools::Itertools;
use rustc_hash::FxHashMap;

use super::{
    ann::{Allele, AnnField, Consequence, FeatureBiotype, FeatureType, Pos, Rank, SoFeature},
    provider::Provider as MehariProvider,
};

/// A variant description how VCF would do it.
#[derive(Debug, PartialEq, Eq, Clone, Default)]
pub struct VcfVariant {
    /// Chromosome name.
    pub chromosome: String,
    /// 1-based position on the chromosome of first base of `reference`.
    pub position: i32,
    /// Reference bases.
    pub reference: String,
    /// Alternative bases.
    pub alternative: String,
}

/// Enum that allows to select the transcript source.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Default,
    serde::Deserialize,
    serde::Serialize,
    clap::ValueEnum,
)]
pub enum TranscriptSource {
    /// ENSEMBL
    Ensembl,
    /// RefSeq
    RefSeq,
    /// Both
    #[default]
    Both,
}

/// Configuration for consequence prediction.
#[derive(Debug, Clone, derive_builder::Builder)]
#[builder(pattern = "immutable")]
pub struct Config {
    /// The transcript source to use.
    #[builder(default = "TranscriptSource::Both")]
    pub transcript_source: TranscriptSource,

    /// Whether to report only the worst consequence for each picked transcript.
    #[builder(default = "false")]
    pub report_worst_consequence_only: bool,

    /// Whether to discard intergenic variants.
    #[builder(default = "true")]
    pub discard_intergenic: bool,
}

impl Default for Config {
    fn default() -> Self {
        ConfigBuilder::default().build().unwrap()
    }
}

/// Wrap mapper, provider, and map for consequence prediction.
#[derive(derivative::Derivative)]
#[derivative(Debug)]
pub struct ConsequencePredictor {
    /// The internal transcript provider for locating transcripts.
    #[derivative(Debug = "ignore")]
    provider: Arc<MehariProvider>,
    /// Assembly mapper for variant consequence prediction.
    #[derivative(Debug = "ignore")]
    mapper: assembly::Mapper,
    /// Mapping from chromosome name to accession.
    #[derivative(Debug = "ignore")]
    chrom_to_acc: HashMap<String, String>,
    /// Configuration for the predictor.
    #[derivative(Debug = "ignore")]
    config: Config,
}

/// Padding to look for genes upstream/downstream.
pub const PADDING: i32 = 5_000;
/// Generally used alternative alignment method.
pub const ALT_ALN_METHOD: &str = "splign";

impl ConsequencePredictor {
    pub fn new(provider: Arc<MehariProvider>, assembly: Assembly, config: Config) -> Self {
        tracing::info!("Building transcript interval trees ...");
        let acc_to_chrom: indexmap::IndexMap<String, String> = provider.get_assembly_map(assembly);
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

        let mapper_config = assembly::Config {
            replace_reference: false,
            strict_bounds: false,
            renormalize_g: false,
            genome_seq_available: false,
            ..Default::default()
        };
        let mapper = assembly::Mapper::new(mapper_config, provider.clone());
        tracing::info!("... done building transcript interval trees");

        ConsequencePredictor {
            provider,
            mapper,
            chrom_to_acc,
            config,
        }
    }

    /// Predict the consequences of a variant.
    ///
    /// Note that the predictions will be affected by whether transcript picking has been
    /// enabled in the data provider and the configuration of the predictor, in particular
    /// `Config::report_all_transcripts`.
    ///
    /// # Args
    ///
    /// * `var`: The variant to predict consequences for.
    ///
    /// # Returns
    ///
    /// A list of `AnnField` records, one for each transcript affected by the variant
    /// sorted lexicographically by transcript accession.
    ///
    /// If the accessio is not valid, then `None` will be returned.
    ///
    /// # Errors
    ///
    /// If there was any error during the prediction.
    pub fn predict(&self, var: &VcfVariant) -> Result<Option<Vec<AnnField>>, anyhow::Error> {
        // Normalize variant by stripping common prefix and suffix.
        let norm_var = self.normalize_variant(var);

        // Obtain accession from chromosome name.
        let chrom_acc = self.chrom_to_acc.get(&norm_var.chromosome);
        let chrom_acc = if let Some(chrom_acc) = chrom_acc {
            chrom_acc
        } else {
            tracing::warn!(
                "Could not determine chromosome accession for {:?}; giving up on annotation",
                &norm_var
            );
            return Ok(None);
        };

        // Get all affected transcripts.
        let (var_start, var_end) = if norm_var.reference.is_empty() {
            (norm_var.position - 1, norm_var.position - 1)
        } else {
            (
                norm_var.position - 1,
                norm_var.position + norm_var.reference.len() as i32 - 1,
            )
        };
        let qry_start = var_start - PADDING;
        let qry_end = var_end + PADDING;
        let txs = {
            let mut txs =
                self.provider
                    .get_tx_for_region(chrom_acc, ALT_ALN_METHOD, qry_start, qry_end)?;
            txs.sort_by(|a, b| a.tx_ac.cmp(&b.tx_ac));
            // Filter transcripts to the picked ones from the selected
            // transcript source.
            self.filter_picked_sourced_txs(txs)
        };

        // Handle case of no overlapping transcripts -> intergenic.
        if txs.is_empty() {
            return Ok(Some(self.filter_ann_fields(vec![AnnField {
                allele: Allele::Alt {
                    alternative: var.alternative.clone(),
                },
                gene_id: "".to_string(),
                consequences: vec![Consequence::IntergenicVariant],
                putative_impact: crate::annotate::seqvars::ann::PutativeImpact::Modifier,
                feature_type: FeatureType::Custom {
                    value: "Intergenic".to_string(),
                },
                feature_id: "".to_string(),
                feature_biotype: vec![FeatureBiotype::Noncoding],
                rank: None,
                distance: None,
                strand: 0,
                hgvs_t: None,
                hgvs_p: None,
                tx_pos: None,
                cds_pos: None,
                protein_pos: None,
                gene_symbol: "".to_string(),
                messages: None,
            }])));
        }

        // Compute annotations for all (picked) transcripts first, skipping `None`` results.
        let anns_all_txs = txs
            .into_iter()
            .map(|tx| {
                self.build_ann_field(var, &norm_var, tx, chrom_acc.clone(), var_start, var_end)
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        // Return all or worst annotation only.
        Ok(Some(self.filter_ann_fields(anns_all_txs)))
    }

    // Filter transcripts to the picked ones from the selected transcript source.
    fn filter_picked_sourced_txs(&self, txs: Vec<TxForRegionRecord>) -> Vec<TxForRegionRecord> {
        fn is_ensembl(tx: &TxForRegionRecord) -> bool {
            tx.tx_ac.starts_with("ENST")
        }

        let txs = match self.config.transcript_source {
            TranscriptSource::Ensembl => txs.into_iter().filter(is_ensembl).collect::<Vec<_>>(),
            TranscriptSource::RefSeq => txs
                .into_iter()
                .filter(|tx| !is_ensembl(tx))
                .collect::<Vec<_>>(),
            TranscriptSource::Both => txs,
        };

        // Short-circuit if transcript picking has been disabled.
        if !self.provider.transcript_picking() {
            return txs;
        }

        // Get gene ids for all transcripts in `txs`, then obtain the picked transcript
        // identifiers for these genes and limit `txs` to those transcripts.
        let picked_txs = txs
            .iter()
            .flat_map(|tx| self.provider.get_tx(&tx.tx_ac))
            .flat_map(|tx| self.provider.get_picked_transcripts(&tx.gene_id))
            .flatten()
            .collect::<Vec<_>>();
        // tracing::trace!("Picked transcripts: {:?}", &picked_txs);
        txs.into_iter()
            .filter(|tx| picked_txs.contains(&tx.tx_ac))
            .collect::<Vec<_>>()
    }

    /// Filter the ANN fields depending on the configuration.
    ///
    /// If all transcripts are to be reported then return `ann_fields` as is, otherwise
    /// select one worst consequence per gene.
    fn filter_ann_fields(&self, ann_fields: Vec<AnnField>) -> Vec<AnnField> {
        let ann_fields = if self.config.discard_intergenic {
            ann_fields
                .into_iter()
                .filter(|field| field.consequences != [Consequence::IntergenicVariant])
                .collect()
        } else {
            ann_fields
        };

        // Short-circuit if to report all transcript results.
        if !self.config.report_worst_consequence_only {
            return ann_fields;
        }

        // First, split annotations by gene.
        let mut anns_by_gene: FxHashMap<String, Vec<AnnField>> = FxHashMap::default();
        for ann in ann_fields {
            let gene_id = ann.gene_id.clone();
            anns_by_gene.entry(gene_id).or_default().push(ann);
        }

        /// Return sort order for ANN biotype, gives priority to ManeSelect and ManePlusClinical.
        fn biotype_order(biotypes: &[FeatureBiotype]) -> i32 {
            if biotypes.contains(&FeatureBiotype::ManeSelect) {
                0
            } else if biotypes.contains(&FeatureBiotype::ManePlusClinical) {
                1
            } else {
                2
            }
        }

        // Now, sort by consequence, giving priority to ManeSelect and ManePlusClinical.
        //
        // This uses the invariant that the consequences in the ANN fields are sorted already
        // and there is at least one consequence.
        let mut result = Vec::new();
        for anns in anns_by_gene.values_mut() {
            anns.sort_by_key(|ann| (ann.consequences[0], biotype_order(&ann.feature_biotype)));
            result.push(anns.remove(0));
        }
        result
    }

    fn build_ann_field(
        &self,
        orig_var: &VcfVariant,
        var: &VcfVariant,
        tx_record: TxForRegionRecord,
        chrom_acc: String,
        var_start: i32,
        var_end: i32,
    ) -> Result<Option<AnnField>, anyhow::Error> {
        // NB: The coordinates of var_start, var_end, as well as the exon boundaries
        // are 0-based.

        let tx = self.provider.get_tx(&tx_record.tx_ac);
        let tx = if let Some(tx) = tx {
            // Skip transcripts that are protein coding but do not have a CDS.
            // TODO: do not include such transcripts when building the database.
            if TranscriptBiotype::try_from(tx.biotype).expect("invalid tx biotype")
                == TranscriptBiotype::Coding
                && tx.start_codon.is_none()
            {
                return Ok(None);
            }
            tx
        } else {
            tracing::warn!(
                "Requested transcript accession {}, got None (potentially filtered)",
                &tx_record.tx_ac
            );
            return Ok(None);
        };

        assert_eq!(
            tx.genome_alignments.len(),
            1,
            "At this point, only one genome alignment is expected"
        );

        // TODO: Report selenocysteine modifications
        // let is_seleno = tx.tags.contains(&(TranscriptTag::Selenoprotein as i32));

        let alignment = tx.genome_alignments.first().unwrap();
        let strand = Strand::try_from(alignment.strand).expect("invalid strand");

        let mut consequences: BitFlags<Consequence> = BitFlags::empty();

        let mut min_start = None;
        let mut max_end = None;

        // Find first exon that overlaps with variant or intron that contains the variant.
        //
        // Note that exons are stored in genome position order.
        let mut prev_end = None;
        let mut rank = Rank::default();
        let mut is_exonic = false;
        let mut is_intronic = false;
        let mut distance: Option<i32> = None;
        let mut tx_len = 0;

        fn overlaps(
            var_start: i32,
            var_end: i32,
            exon_intron_start: i32,
            exon_intron_end: i32,
        ) -> bool {
            (var_start < exon_intron_end) && (var_end > exon_intron_start)
        }
        let var_overlaps =
            |start: i32, end: i32| -> bool { overlaps(var_start, var_end, start, end) };

        for exon_alignment in &alignment.exons {
            tx_len += exon_alignment.alt_end_i - exon_alignment.alt_start_i;

            let exon_start = exon_alignment.alt_start_i;
            let exon_end = exon_alignment.alt_end_i;
            let intron_start = prev_end;
            let intron_end = exon_start;

            // Check the cases where the variant overlaps with the exon or is contained within an
            // intron.
            if var_overlaps(exon_start, exon_end) {
                // overlaps with exon
                rank = Rank {
                    ord: exon_alignment.ord + 1,
                    total: alignment.exons.len() as i32,
                };
                is_exonic = true;
                distance = Some(0);
            } else if let Some(intron_start) = intron_start {
                // We are in an intron (the first exon does not have an intron left of it
                // which is expressed by `intron_start` being an `Option<i32>` rather than `i32`.
                if var_start >= intron_start && var_end <= intron_end {
                    // Contained within intron: cannot be in next exon.
                    if !is_exonic {
                        rank = Rank {
                            ord: exon_alignment.ord + 1,
                            total: alignment.exons.len() as i32 - 1,
                        };
                        is_intronic = true;

                        // We compute the "distance" with "+1", the first base of the
                        // intron is "+1", the last one is "-1".
                        let dist_start = var_start + 1 - intron_start;
                        let dist_end = -(intron_end + 1 - var_end);
                        let dist_start_end = if dist_start.abs() <= dist_end.abs() {
                            dist_start
                        } else {
                            dist_end
                        };
                        if distance.is_none()
                            || dist_start_end.abs() <= distance.expect("cannot be None").abs()
                        {
                            distance = Some(dist_start_end);
                        }
                    }
                }
            }

            // Check the cases where the variant overlaps with whole exon.
            if var_start <= exon_start && var_end >= exon_end {
                consequences |= Consequence::ExonLossVariant;
                if var_start < exon_start {
                    if strand == Strand::Plus && !rank.is_first() {
                        consequences |= Consequence::SpliceAcceptorVariant;
                    } else if strand == Strand::Minus && !rank.is_last() {
                        consequences |= Consequence::SpliceDonorVariant;
                    }
                }
                if var_end > exon_end {
                    if strand == Strand::Plus && !rank.is_last() {
                        consequences |= Consequence::SpliceDonorVariant;
                    } else if strand == Strand::Minus && !rank.is_last() {
                        consequences |= Consequence::SpliceAcceptorVariant;
                    }
                }
            }
            if let Some(intron_start) = intron_start {
                // For insertions, we need to consider the case of the insertion being right at
                // the exon/intron junction.  We can express this with a shift of 1 for using
                // "< / >" X +/- shift and meaning <= / >= X.
                let ins_shift = if var.reference.is_empty() { 1 } else { 0 };

                // Check the cases where the variant overlaps with the splice acceptor/donor site.
                if var_overlaps(intron_start - ins_shift, intron_start + 2) {
                    // Left side, is acceptor/donor depending on transcript's strand.
                    match strand {
                        Strand::Plus => consequences.insert(Consequence::SpliceDonorVariant),
                        Strand::Minus => consequences.insert(Consequence::SpliceAcceptorVariant),
                        _ => unreachable!("invalid strand: {}", alignment.strand),
                    }
                }
                // Check the case where the variant overlaps with the splice donor site.
                if var_overlaps(intron_end - 2, intron_end + ins_shift) {
                    // Left side, is acceptor/donor depending on transcript's strand.
                    match strand {
                        Strand::Plus => consequences.insert(Consequence::SpliceAcceptorVariant),
                        Strand::Minus => consequences.insert(Consequence::SpliceDonorVariant),
                        _ => unreachable!("invalid strand: {}", alignment.strand),
                    }
                }
            }
            // Check the case where the variant overlaps with the splice region (1-3 bases in exon
            // or 3-8 bases in intron).  We have to check all cases independently and not with `else`
            // because the variant may be larger.
            if let Some(intron_start) = intron_start {
                if var_overlaps(intron_start + 2, intron_start + 8)
                    || var_overlaps(intron_end - 8, intron_end - 2)
                {
                    consequences |= Consequence::SpliceRegionVariant;
                }
                if var_overlaps(exon_end - 3, exon_end) {
                    if strand == Strand::Plus {
                        if !rank.is_last() {
                            consequences |= Consequence::SpliceRegionVariant;
                        }
                    } else {
                        // alignment.strand == Strand::Minus
                        if !rank.is_first() {
                            consequences |= Consequence::SpliceRegionVariant;
                        }
                    }
                }
                if var_overlaps(exon_start, exon_start + 3) {
                    if strand == Strand::Plus {
                        if !rank.is_first() {
                            consequences |= Consequence::SpliceRegionVariant;
                        }
                    } else {
                        // alignment.strand == Strand::Minus
                        if !rank.is_last() {
                            consequences |= Consequence::SpliceRegionVariant;
                        }
                    }
                }
            }

            if let Some(intron_start) = intron_start {
                // Check the case where the variant overlaps with the polypyrimidine tract.
                // (A sequence variant that falls in the polypyrimidine tract at 3' end of intron between 17 and 3 bases from the end (acceptor -3 to acceptor -17))
                if strand == Strand::Plus && var_overlaps(intron_end - 17, intron_end - 2) {
                    consequences |= Consequence::SplicePolypyrimidineTractVariant;
                }
                if strand == Strand::Minus && var_overlaps(intron_start + 2, intron_start + 17) {
                    consequences |= Consequence::SplicePolypyrimidineTractVariant;
                }
                // Check conditions for splice_donor_region_variant
                // (A sequence variant that falls in the region between the 3rd and 6th base after splice junction (5' end of intron))
                if strand == Strand::Plus && var_overlaps(intron_start + 2, intron_start + 6) {
                    consequences |= Consequence::SpliceDonorRegionVariant;
                }
                if strand == Strand::Minus && var_overlaps(intron_end - 6, intron_end - 2) {
                    consequences |= Consequence::SpliceDonorRegionVariant;
                }
                // Check conditions for splice_donor_5th_base_variant
                // (A sequence variant that causes a change at the 5th base pair after the start of the intron in the orientation of the transcript.)
                if strand == Strand::Plus && var_overlaps(intron_start + 4, intron_start + 5) {
                    consequences |= Consequence::SpliceDonorFifthBaseVariant;
                }
                if strand == Strand::Minus && var_overlaps(intron_end - 5, intron_end - 4) {
                    consequences |= Consequence::SpliceDonorFifthBaseVariant;
                }
            }

            min_start = Some(std::cmp::min(min_start.unwrap_or(exon_start), exon_start));
            max_end = Some(std::cmp::max(max_end.unwrap_or(exon_end), exon_end));

            prev_end = Some(exon_end);
        }

        let min_start = min_start.expect("must have seen exon");
        let max_end = max_end.expect("must have seen exon");

        let transcript_biotype =
            TranscriptBiotype::try_from(tx.biotype).expect("invalid transcript biotype");
        let feature_biotype = {
            let mut feature_biotypes = vec![match transcript_biotype {
                TranscriptBiotype::Coding => FeatureBiotype::Coding,
                TranscriptBiotype::NonCoding => FeatureBiotype::Noncoding,
                _ => unreachable!("invalid biotype: {:?}", transcript_biotype),
            }];

            if tx.tags.contains(&(TranscriptTag::ManeSelect as i32)) {
                feature_biotypes.push(FeatureBiotype::ManeSelect);
            } else if tx.tags.contains(&(TranscriptTag::ManePlusClinical as i32)) {
                feature_biotypes.push(FeatureBiotype::ManePlusClinical);
            }

            feature_biotypes
        };

        let is_upstream = var_end <= min_start;
        let is_downstream = var_start >= max_end;
        if is_exonic {
            if transcript_biotype == TranscriptBiotype::NonCoding {
                consequences |= Consequence::NonCodingTranscriptExonVariant;
            }
        } else if is_intronic {
            if transcript_biotype == TranscriptBiotype::NonCoding {
                consequences |= Consequence::NonCodingTranscriptIntronVariant;
            } else {
                consequences |= Consequence::IntronVariant;
            }
        } else if is_upstream {
            let val = -(min_start + 1 - var_end);
            if val.abs() <= PADDING {
                match strand {
                    Strand::Plus => consequences |= Consequence::UpstreamGeneVariant,
                    Strand::Minus => consequences |= Consequence::DownstreamGeneVariant,
                    _ => unreachable!("invalid strand: {}", alignment.strand),
                }
            } else {
                unreachable!("variant is intergenic but this cannot happen here");
            }
            if distance.is_none() {
                distance = Some(val);
            }
        } else if is_downstream {
            let val = var_start + 1 - max_end;
            if val.abs() <= PADDING {
                match strand {
                    Strand::Plus => consequences |= Consequence::DownstreamGeneVariant,
                    Strand::Minus => consequences |= Consequence::UpstreamGeneVariant,
                    _ => unreachable!("invalid strand: {}", alignment.strand),
                }
            } else {
                unreachable!("variant is intergenic but this cannot happen here");
            }
            if distance.is_none() {
                distance = Some(val);
            }
        }

        let var_g = Self::get_var_g(var, &chrom_acc);

        let (rank, hgvs_t, hgvs_p, tx_pos, cds_pos, protein_pos) = if !is_upstream && !is_downstream
        {
            // TODO: do not include such transcripts when building the tx database.
            let var_n = self.mapper.g_to_n(&var_g, &tx.id).map_or_else(
                |e| match e {
                    Error::NonAdjacentExons(_, _, _, _) => {
                        tracing::warn!("{}, {}: NonAdjacentExons, skipping", &tx.id, &var_g);
                        Ok(None)
                    }
                    _ => Err(e),
                },
                |v| Ok(Some(v)),
            )?;
            if var_n.is_none() {
                return Ok(None);
            }
            let var_n = var_n.unwrap();

            let tx_pos = match &var_n {
                HgvsVariant::TxVariant { loc_edit, .. } => Some(Pos {
                    ord: loc_edit.loc.inner().start.base,
                    total: Some(tx_len),
                }),
                _ => panic!("Invalid tx position: {:?}", &var_n),
            };

            let (var_t, _var_p, hgvs_p, cds_pos, protein_pos) = match transcript_biotype {
                TranscriptBiotype::Coding => {
                    let cds_len = tx.stop_codon.unwrap() - tx.start_codon.unwrap();
                    let prot_len = cds_len / 3;

                    let var_c = self.mapper.n_to_c(&var_n)?;
                    // Gracefully handle the case that the transcript is unsupported because the length
                    // is not a multiple of 3.
                    // TODO: do not include such transcripts when building the tx database.
                    let var_p = self.mapper.c_to_p(&var_c).map_or_else(
                        |e| {
                            if e.to_string()
                                .contains("is not supported because its sequence length of")
                            {
                                Ok(None)
                            } else {
                                Err(e)
                            }
                        },
                        |v| Ok(Some(v)),
                    )?;
                    if var_p.is_none() {
                        return Ok(None);
                    }
                    let var_p = var_p.unwrap();

                    let hgvs_p = format!("{}", &var_p);
                    let hgvs_p = hgvs_p.split(':').nth(1).unwrap().to_owned();
                    let hgvs_p = Some(hgvs_p);
                    let cds_pos = match &var_c {
                        HgvsVariant::CdsVariant { loc_edit, .. } => Some(Pos {
                            ord: loc_edit.loc.inner().start.base,
                            total: Some(cds_len),
                        }),
                        _ => panic!("Invalid CDS position: {:?}", &var_n),
                    };
                    let protein_pos = match &var_p {
                        HgvsVariant::ProtVariant { loc_edit, .. } => match &loc_edit {
                            ProtLocEdit::Ordinary { loc, .. } => Some(Pos {
                                ord: loc.inner().start.number,
                                total: Some(prot_len),
                            }),
                            _ => None,
                        },
                        _ => panic!("Not a protein position: {:?}", &var_n),
                    };

                    let conservative = match &var_c {
                        HgvsVariant::CdsVariant { loc_edit, .. } => {
                            // Handle the cases where the variant touches the start or stop codon based on `var_c`
                            // coordinates.  The cases where the start/stop codon is touched by the variant
                            // directly is handled above based on the `var_p` prediction.
                            let loc = loc_edit.loc.inner();
                            let start_base = loc.start.base;
                            let start_cds_from = loc.start.cds_from;
                            let end_base = loc.end.base;
                            let end_cds_from = loc.end.cds_from;
                            // The variables below mean "VARIANT_{starts,stops}_{left,right}_OF_{start,stop}_CODON".
                            //
                            // start codon
                            let starts_left_of_start =
                                start_cds_from == CdsFrom::Start && start_base < 0;
                            let ends_right_of_start =
                                start_cds_from != CdsFrom::Start || start_base > 0;
                            if starts_left_of_start && ends_right_of_start {
                                consequences |= Consequence::StartLost;
                            }
                            // stop codon
                            let starts_left_of_stop = start_cds_from == CdsFrom::Start;
                            let ends_right_of_stop = end_cds_from == CdsFrom::End;
                            if starts_left_of_stop && ends_right_of_stop {
                                consequences |= Consequence::StopLost;
                            }

                            // Detect variants affecting the 5'/3' UTRs.
                            if start_cds_from == CdsFrom::Start {
                                if start_base < 0 {
                                    if is_intronic {
                                        consequences |= Consequence::FivePrimeUtrIntronVariant;
                                    }
                                    if is_exonic {
                                        consequences |= Consequence::FivePrimeUtrExonVariant;
                                    }
                                }
                            } else if end_cds_from == CdsFrom::End {
                                if is_intronic {
                                    consequences |= Consequence::ThreePrimeUtrIntronVariant;
                                }
                                if is_exonic {
                                    consequences |= Consequence::ThreePrimeUtrExonVariant;
                                }
                            }

                            // The range is "conservative" (regarding deletions and insertions) if
                            // it does not start or end within exons.
                            start_cds_from == CdsFrom::Start
                                && end_cds_from == CdsFrom::Start
                                && start_base % 3 == 1
                                && end_base % 3 == 1
                        }
                        _ => panic!("Must be CDS variant: {}", &var_c),
                    };

                    fn is_stop(s: &str) -> bool {
                        s == "X" || s == "Ter" || s == "*"
                    }

                    // Analyze `var_p` for changes in the protein sequence.
                    match &var_p {
                        HgvsVariant::ProtVariant { loc_edit, .. } => match loc_edit {
                            ProtLocEdit::Ordinary { loc, edit } => {
                                let loc = loc.inner();
                                match edit.inner() {
                                    hgvs::parser::ProteinEdit::Fs { .. } => {
                                        consequences |= Consequence::FrameshiftVariant;
                                    }
                                    hgvs::parser::ProteinEdit::Ext { .. } => {
                                        consequences |= Consequence::StopLost;
                                        consequences |= Consequence::FeatureElongation;
                                    }
                                    hgvs::parser::ProteinEdit::Subst { alternative } => {
                                        if alternative.is_empty() {
                                            consequences |= Consequence::SynonymousVariant;
                                        } else if is_stop(alternative) {
                                            if loc.start == loc.end && is_stop(&loc.start.aa) {
                                                consequences |= Consequence::StopRetainedVariant;
                                            } else {
                                                consequences |= Consequence::StopGained;
                                            }
                                        } else {
                                            consequences |= Consequence::MissenseVariant;
                                        }
                                    }
                                    hgvs::parser::ProteinEdit::DelIns { alternative } => {
                                        if conservative {
                                            consequences |=
                                                Consequence::ConservativeInframeDeletion;
                                        } else {
                                            consequences |= Consequence::DisruptiveInframeDeletion;
                                        }
                                        if alternative.contains('*')
                                            || alternative.contains('X')
                                            || alternative.contains("Ter")
                                        {
                                            consequences |= Consequence::StopGained;
                                        }
                                    }
                                    hgvs::parser::ProteinEdit::Ins { .. }
                                    | hgvs::parser::ProteinEdit::Dup => {
                                        if conservative {
                                            consequences |=
                                                Consequence::ConservativeInframeInsertion;
                                        } else {
                                            consequences |= Consequence::DisruptiveInframeInsertion;
                                        }
                                        consequences |= Consequence::ConservativeInframeInsertion;
                                    }
                                    hgvs::parser::ProteinEdit::Del => {
                                        if conservative {
                                            consequences |=
                                                Consequence::ConservativeInframeDeletion;
                                        } else {
                                            consequences |= Consequence::DisruptiveInframeDeletion;
                                        }
                                    }
                                    hgvs::parser::ProteinEdit::Ident => {
                                        consequences |= Consequence::SynonymousVariant;
                                    }
                                };
                            }
                            ProtLocEdit::NoChange | ProtLocEdit::NoChangeUncertain => {
                                consequences |= Consequence::SynonymousVariant;
                            }
                            ProtLocEdit::InitiationUncertain => {
                                consequences |= Consequence::StartLost;
                            }
                            ProtLocEdit::NoProtein
                            | ProtLocEdit::NoProteinUncertain
                            | ProtLocEdit::Unknown => (),
                        },
                        _ => panic!("Must be protein variant: {}", &var_p),
                    }

                    (var_c, Some(var_p), hgvs_p, cds_pos, protein_pos)
                }
                TranscriptBiotype::NonCoding => (var_n, None, None, None, None),
                _ => unreachable!("invalid transcript biotype: {:?}", transcript_biotype),
            };
            let hgvs_t = format!("{}", &NoRef(&var_t));
            let hgvs_t = hgvs_t.split(':').nth(1).unwrap().to_owned();

            (
                Some(rank),
                Some(hgvs_t),
                hgvs_p,
                tx_pos,
                cds_pos,
                protein_pos,
            )
        } else {
            (None, None, None, None, None, None)
        };

        // Take a highest-ranking consequence and derive putative impact from it.
        if consequences.is_empty() {
            tracing::error!(
                "No consequences for {:?} on {} (hgvs_p={}) - adding `gene_variant`;\
                most likely the transcript has multiple stop codons and the variant \
                lies behind the first.",
                var,
                &tx_record.tx_ac,
                hgvs_p
                    .as_ref()
                    .map(|s| (&s).to_string())
                    .unwrap_or(String::from("None"))
            );
            consequences |= Consequence::GeneVariant;
        }
        let consequences = consequences.iter().collect_vec();
        let putative_impact = (*consequences.first().unwrap()).into();

        let strand = match strand {
            Strand::Unknown => 0,
            Strand::Plus => 1,
            Strand::Minus => -1,
        };

        // Build and return ANN field from the information derived above.
        Ok(Some(AnnField {
            allele: Allele::Alt {
                alternative: orig_var.alternative.clone(),
            },
            consequences,
            putative_impact,
            gene_symbol: tx.gene_symbol,
            gene_id: tx.gene_id,
            feature_type: FeatureType::SoTerm {
                term: SoFeature::Transcript,
            },
            feature_id: tx.id,
            feature_biotype,
            rank,
            hgvs_t,
            hgvs_p,
            tx_pos,
            cds_pos,
            protein_pos,
            strand,
            distance,
            messages: None,
        }))
    }

    fn get_var_g(var: &VcfVariant, chrom_acc: &str) -> HgvsVariant {
        let chrom_acc = chrom_acc.to_string();
        HgvsVariant::GenomeVariant {
            accession: Accession { value: chrom_acc },
            gene_symbol: None,
            loc_edit: if var.reference.is_empty() {
                // insertion
                GenomeLocEdit {
                    loc: Mu::Certain(GenomeInterval {
                        start: Some(var.position - 1),
                        end: Some(var.position),
                    }),
                    edit: Mu::Certain(NaEdit::Ins {
                        alternative: var.alternative.clone(),
                    }),
                }
            } else if var.alternative.is_empty() {
                // deletion
                GenomeLocEdit {
                    loc: Mu::Certain(GenomeInterval {
                        start: Some(var.position),
                        end: Some(var.position + var.reference.len() as i32 - 1),
                    }),
                    edit: Mu::Certain(NaEdit::DelRef {
                        reference: var.reference.clone(),
                    }),
                }
            } else {
                // substitution
                GenomeLocEdit {
                    loc: Mu::Certain(GenomeInterval {
                        start: Some(var.position),
                        end: Some(var.position + var.reference.len() as i32 - 1),
                    }),
                    edit: Mu::Certain(NaEdit::RefAlt {
                        reference: var.reference.clone(),
                        alternative: var.alternative.clone(),
                    }),
                }
            },
        }
    }

    // Normalize variant by stripping common suffixes and prefixes.
    fn normalize_variant(&self, var: &VcfVariant) -> VcfVariant {
        let mut result = var.clone();

        // Strip common suffixes.
        while result.reference.len() > 1 && result.alternative.len() > 1 {
            if result.reference.chars().last().unwrap()
                == result.alternative.chars().last().unwrap()
            {
                result.reference.pop();
                result.alternative.pop();
            } else {
                break;
            }
        }

        // Strip common suffixes.
        while !result.reference.is_empty() && !result.alternative.is_empty() {
            if result.reference.chars().next().unwrap()
                == result.alternative.chars().next().unwrap()
            {
                result.position += 1;
                result.reference.remove(0);
                result.alternative.remove(0);
            } else {
                break;
            }
        }

        result
    }
}

impl ConsequencePredictor {
    /// Return data version string (if set).
    pub fn data_version(&self) -> Option<String> {
        self.provider.as_ref().tx_seq_db.version.clone()
    }
}

#[cfg(test)]
mod test {
    use std::{fs::File, io::BufReader};

    use csv::ReaderBuilder;
    use pretty_assertions::assert_eq;
    use serde::Deserialize;

    use crate::annotate::seqvars::provider::ConfigBuilder as MehariProviderConfigBuilder;
    use crate::annotate::seqvars::{load_tx_db, TranscriptPickType};

    use super::*;

    #[test]
    fn test_sync() {
        fn is_sync<T: Sync>() {}
        is_sync::<super::ConsequencePredictor>();
    }

    #[rstest::rstest]
    #[case("17:41197701:G:C", 0)] // exonic
    #[case("17:41196309:G:C", -3)] // 3bp 3' upstream
    #[case("17:41196310:G:C", -2)] // 2bp 3' upstream
    #[case("17:41196311:G:C", -1)] // 1bp 3' upstream
    #[case("17:41196312:G:C", 0)] // ex. 3' UTR
    #[case("17:41196313:G:C", 0)] // ex. 3' UTR
    #[case("17:41197818:G:C", 0)] // exonic
    #[case("17:41197819:G:C", 0)] // exonic
    #[case("17:41197820:G:C", 1)] // 1bp intronic
    #[case("17:41197821:G:C", 2)] // 2bp intronic
    #[case("17:41197822:G:C", 3)] // 3bp intronic
    #[case("17:41197823:G:C", 4)] // 4bp intronic
    #[case("17:41277379:A:C", 0)] // exonic
    #[case("17:41277380:G:C", 0)] // exonic
    #[case("17:41277381:G:T", 0)] // exonic
    #[case("17:41277382:G:C", 1)] // 1bp upstream
    #[case("17:41277383:A:C", 2)] // 2bp upstream
    #[case("17:41277384:G:C", 3)] // 3bp upstream
    fn annotate_snv_brca1_one_variant(
        #[case] spdi: &str,
        #[case] expected_dist: i32,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", spdi.replace(':', "-"));

        let spdi = spdi.split(':').map(|s| s.to_string()).collect::<Vec<_>>();

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            Assembly::Grch37p10,
            Default::default(),
        ));

        let predictor =
            ConsequencePredictor::new(provider, Assembly::Grch37p10, Default::default());

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].clone(),
                position: spdi[1].parse()?,
                reference: spdi[2].clone(),
                alternative: spdi[3].clone(),
            })?
            .unwrap();

        assert_eq!(res.len(), 5);
        insta::assert_yaml_snapshot!(res);
        assert_eq!(
            res[0].distance,
            Some(expected_dist),
            "spdi = {}",
            spdi.join(":")
        );

        Ok(())
    }

    /// Test some intron specific variants, via the annotated consequences.
    /// GRCh37, BRCA1, NM_007294.4 (MANE, reverse).
    /// The order of the consequences is important: ordered by severity, descending.
    /// cf Consequences enum ordering.
    #[rstest::rstest]
    #[case("17:41197820:G:T", 1, vec![Consequence::SpliceAcceptorVariant, Consequence::IntronVariant]
    )] // 1bp intronic
    #[case("17:41197821:A:C", 2, vec![Consequence::SpliceAcceptorVariant, Consequence::IntronVariant]
    )] // 2bp intronic
    #[case("17:41197822:C:A", 3, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::IntronVariant]
    )] // 3bp intronic
    #[case("17:41197823:C:A", 4, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::IntronVariant]
    )] // 4bp intronic
    #[case("17:41197824:T:G", 5, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::IntronVariant]
    )] // 5bp intronic
    #[case("17:41197825:C:A", 6, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::IntronVariant]
    )] // 6bp intronic
    #[case("17:41197835:T:G", 16, vec![Consequence::SplicePolypyrimidineTractVariant, Consequence::IntronVariant]
    )] // 16bp intronic
    #[case("17:41197836:G:A", 17, vec![Consequence::SplicePolypyrimidineTractVariant, Consequence::IntronVariant]
    )] // 17bp intronic
    #[case("17:41197837:G:A", 18, vec![Consequence::IntronVariant])] // 18bp intronic
    #[case("17:41199660:G:T", 0, vec![Consequence::MissenseVariant, Consequence::SpliceRegionVariant]
    )] // exonic
    #[case("17:41199659:G:T", -1, vec![Consequence::SpliceDonorVariant, Consequence::IntronVariant]
    )] // -1bp intronic
    #[case("17:41199658:T:G", -2, vec![Consequence::SpliceDonorVariant, Consequence::IntronVariant]
    )] // -2bp intronic
    #[case("17:41199657:G:T", -3, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::IntronVariant]
    )] // -3bp intronic
    #[case("17:41199656:A:C", -4, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::IntronVariant]
    )] // -4bp intronic
    #[case("17:41199655:G:T", -5, vec![Consequence::SpliceDonorFifthBaseVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::IntronVariant]
    )] // -5bp intronic
    #[case("17:41199654:G:T", -6, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::IntronVariant]
    )] // -6bp intronic
    #[case("17:41199653:T:G", -7, vec![Consequence::SpliceRegionVariant, Consequence::IntronVariant]
    )] // -7bp intronic
    #[case("17:41199652:G:T", -8, vec![Consequence::SpliceRegionVariant, Consequence::IntronVariant]
    )] // -8bp intronic
    #[case("17:41199651:C:A", -9, vec![Consequence::IntronVariant])] // -9bp intronic
    fn annotate_snv_brca1_csq(
        #[case] spdi: &str,
        #[case] expected_dist: i32,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", spdi.replace(':', "-"));

        let spdi = spdi.split(':').map(|s| s.to_string()).collect::<Vec<_>>();

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            Assembly::Grch37p10,
            MehariProviderConfigBuilder::default()
                .transcript_picking(vec![
                    TranscriptPickType::ManePlusClinical,
                    TranscriptPickType::ManeSelect,
                    TranscriptPickType::Length,
                ])
                .build()?,
        ));

        use crate::annotate::seqvars::ConsequencePredictorConfigBuilder;
        let predictor = ConsequencePredictor::new(
            provider,
            Assembly::Grch37p10,
            ConsequencePredictorConfigBuilder::default()
                .report_worst_consequence_only(true)
                .build()?,
        );

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].clone(),
                position: spdi[1].parse()?,
                reference: spdi[2].clone(),
                alternative: spdi[3].clone(),
            })?
            .unwrap();

        assert_eq!(res.len(), 1);
        assert_eq!(res[0].feature_id, "NM_007294.4");
        assert_eq!(
            res[0].distance,
            Some(expected_dist),
            "spdi = {}",
            spdi.join(":")
        );
        assert_eq!(
            res[0].consequences,
            expected_csqs,
            "spdi = {}",
            spdi.join(":")
        );
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    /// Test some intron specific variants, via the annotated consequences.
    /// GRCh37, OPA1, NM_130837.3 (MANE, forward).
    /// The order of the consequences is important: ordered by severity, descending.
    /// cf Consequences enum ordering.
    #[rstest::rstest]
    #[case("3:193332512:T:G", 0, vec![Consequence::MissenseVariant, Consequence::SpliceRegionVariant]
    )] // exonic
    #[case("3:193332511:G:T", -1, vec![Consequence::SpliceAcceptorVariant, Consequence::IntronVariant]
    )] // -1bp intronic
    #[case("3:193332510:A:G", -2, vec![Consequence::SpliceAcceptorVariant, Consequence::IntronVariant]
    )] // -2bp intronic
    #[case("3:193332509:C:T", -3, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant,  Consequence::IntronVariant]
    )] // -3bp intronic
    #[case("3:193332508:T:C", -4, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant,  Consequence::IntronVariant]
    )] // -4bp intronic
    #[case("3:193332507:T:C", -5, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant,  Consequence::IntronVariant]
    )] // -5bp intronic
    #[case("3:193332506:T:C", -6, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::IntronVariant]
    )] // -6bp intronic
    #[case("3:193332505:C:G", -7, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::IntronVariant]
    )] // -7bp intronic
    #[case("3:193332504:T:C", -8, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::IntronVariant]
    )] // -8bp intronic
    #[case("3:193332503:T:A", -9, vec![Consequence::SplicePolypyrimidineTractVariant, Consequence::IntronVariant]
    )] // -9bp intronic
    #[case("3:193332831:G:T", 1, vec![Consequence::SpliceDonorVariant, Consequence::IntronVariant]
    )] // 1bp intronic
    #[case("3:193332832:T:C", 2, vec![Consequence::SpliceDonorVariant, Consequence::IntronVariant]
    )] // 2bp intronic
    #[case("3:193332833:G:A", 3, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::IntronVariant]
    )] // 3bp intronic
    #[case("3:193332834:A:C", 4, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::IntronVariant]
    )] // 4bp intronic
    #[case("3:193332835:A:T", 5, vec![Consequence::SpliceDonorFifthBaseVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::IntronVariant]
    )] // 5bp intronic
    #[case("3:193332836:C:A", 6, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::IntronVariant]
    )] // 6bp intronic
    #[case("3:193332837:T:G", 7, vec![Consequence::SpliceRegionVariant, Consequence::IntronVariant]
    )] // 7bp intronic
    #[case("3:193332838:T:G", 8, vec![Consequence::SpliceRegionVariant, Consequence::IntronVariant]
    )] // 8bp intronic
    #[case("3:193332839:G:A", 9, vec![Consequence::IntronVariant])] // 9bp intronic
    #[case("3:193332846:A:G", 16, vec![Consequence::IntronVariant])] // 16bp intronic
    #[case("3:193332847:G:A", 17, vec![Consequence::IntronVariant])] // 17bp intronic
    #[case("3:193332848:T:A", 18, vec![Consequence::IntronVariant])] // 18bp intronic
    fn annotate_snv_opa1_csq(
        #[case] spdi: &str,
        #[case] expected_dist: i32,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", spdi.replace(':', "-"));

        let spdi = spdi.split(':').map(|s| s.to_string()).collect::<Vec<_>>();

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            Assembly::Grch37p10,
            MehariProviderConfigBuilder::default()
                .transcript_picking(vec![
                    TranscriptPickType::ManePlusClinical,
                    TranscriptPickType::ManeSelect,
                    TranscriptPickType::Length,
                ])
                .build()?,
        ));

        let predictor =
            ConsequencePredictor::new(provider, Assembly::Grch37p10, Default::default());

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].clone(),
                position: spdi[1].parse()?,
                reference: spdi[2].clone(),
                alternative: spdi[3].clone(),
            })?
            .unwrap();

        assert_eq!(res.len(), 1);
        assert_eq!(res[0].feature_id, "NM_130837.3");
        assert_eq!(
            res[0].distance,
            Some(expected_dist),
            "spdi = {}",
            spdi.join(":")
        );
        assert_eq!(
            res[0].consequences,
            expected_csqs,
            "spdi = {}",
            spdi.join(":")
        );
        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case("17:41197701:G:C", false, true)] // don't pick transcripts, report worst
    #[case("17:41197701:G:C", false, false)] // don't pick transcripts, report all
    #[case("17:41197701:G:C", true, true)] // pick transcripts, report worst
    #[case("17:41197701:G:C", true, false)] // pick transcripts, report all
    fn annotate_snv_brca1_transcript_picking_reporting(
        #[case] spdi: &str,
        #[case] pick_transcripts: bool,
        #[case] report_worst_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!(
            "{}-{}-{}",
            spdi.replace(':', "-"),
            pick_transcripts,
            !report_worst_consequence_only
        );

        let spdi = spdi.split(':').map(|s| s.to_string()).collect::<Vec<_>>();

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let picks = if pick_transcripts {
            vec![
                TranscriptPickType::ManePlusClinical,
                TranscriptPickType::ManeSelect,
                TranscriptPickType::Length,
            ]
        } else {
            vec![]
        };
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            Assembly::Grch37p10,
            MehariProviderConfigBuilder::default()
                .transcript_picking(picks)
                .build()
                .unwrap(),
        ));

        let predictor = ConsequencePredictor::new(
            provider,
            Assembly::Grch37p10,
            ConfigBuilder::default()
                .report_worst_consequence_only(report_worst_consequence_only)
                .build()
                .unwrap(),
        );

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].clone(),
                position: spdi[1].parse()?,
                reference: spdi[2].clone(),
                alternative: spdi[3].clone(),
            })?
            .unwrap();

        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    // Test predictions on TTN where we have a ManeSelect and a ManePlusClinical
    // transcript.
    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case("2:179631246:G:A", false, true)] // don't pick transcripts, report worst
    #[case("2:179631246:G:A", false, false)] // don't pick transcripts, report all
    #[case("2:179631246:G:A", true, true)] // pick transcripts, report worst
    #[case("2:179631246:G:A", true, false)] // pick transcripts, report all
    fn annotate_snv_ttn_transcript_picking_reporting(
        #[case] spdi: &str,
        #[case] pick_transcripts: bool,
        #[case] report_worst_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!(
            "{}-{}-{}",
            spdi.replace(':', "-"),
            pick_transcripts,
            !report_worst_consequence_only
        );

        let spdi = spdi.split(':').map(|s| s.to_string()).collect::<Vec<_>>();

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;

        let picks = if pick_transcripts {
            vec![
                TranscriptPickType::ManePlusClinical,
                TranscriptPickType::ManeSelect,
                TranscriptPickType::Length,
            ]
        } else {
            vec![]
        };

        let provider = Arc::new(MehariProvider::new(
            tx_db,
            Assembly::Grch37p10,
            MehariProviderConfigBuilder::default()
                .transcript_picking(picks)
                .build()
                .unwrap(),
        ));

        let predictor = ConsequencePredictor::new(
            provider,
            Assembly::Grch37p10,
            ConfigBuilder::default()
                .report_worst_consequence_only(report_worst_consequence_only)
                .build()
                .unwrap(),
        );

        let res = predictor
            .predict(&VcfVariant {
                chromosome: spdi[0].clone(),
                position: spdi[1].parse()?,
                reference: spdi[2].clone(),
                alternative: spdi[3].clone(),
            })?
            .unwrap();

        insta::assert_yaml_snapshot!(res);

        Ok(())
    }

    #[derive(Debug, Deserialize)]
    struct Record {
        pub var: String,
        pub tx: String,
        pub csq: String,
    }

    // Compare to SnpEff annotated variants for OPA1, touching special cases.
    #[test]
    fn annotate_opa1_hand_picked_vars() -> Result<(), anyhow::Error> {
        annotate_opa1_vars("tests/data/annotate/seqvars/opa1.hand_picked.tsv", false)
    }

    // Compare to SnpEff annotated ClinVar variants for OPA1 (slow).
    #[test]
    fn annotate_opa1_clinvar_vars_snpeff() -> Result<(), anyhow::Error> {
        annotate_opa1_vars(
            "tests/data/annotate/seqvars/clinvar.excerpt.snpeff.opa1.tsv",
            false,
        )
    }

    // Compare to SnpEff annotated ClinVar variants for OPA1 (slow).
    #[test]
    fn annotate_opa1_clinvar_vars_vep() -> Result<(), anyhow::Error> {
        annotate_opa1_vars(
            "tests/data/annotate/seqvars/clinvar.excerpt.vep.opa1.tsv",
            false,
        )
    }

    fn annotate_opa1_vars(
        path_tsv: &str,
        report_worst_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        let txs = vec![
            String::from("NM_001354663.2"),
            String::from("NM_001354664.2"),
            String::from("NM_015560.3"),
            String::from("NM_130831.3"),
            String::from("NM_130832.3"),
            String::from("NM_130837.3"),
        ];

        annotate_vars(path_tsv, &txs, report_worst_consequence_only)
    }

    // Compare to SnpEff annotated variants for BRCA1, touching special cases.
    #[test]
    fn annotate_brca1_hand_picked_vars() -> Result<(), anyhow::Error> {
        annotate_brca1_vars("tests/data/annotate/seqvars/brca1.hand_picked.tsv", false)
    }

    // Compare to SnpEff annotated ClinVar variants for BRCA1 (slow).
    #[test]
    fn annotate_brca1_clinvar_vars_snpeff() -> Result<(), anyhow::Error> {
        annotate_brca1_vars(
            "tests/data/annotate/seqvars/clinvar.excerpt.snpeff.brca1.tsv",
            false,
        )
    }

    // Compare to SnpEff annotated ClinVar variants for BRCA1 (slow).
    #[test]
    fn annotate_brca1_clinvar_vars_vep() -> Result<(), anyhow::Error> {
        annotate_brca1_vars(
            "tests/data/annotate/seqvars/clinvar.excerpt.vep.brca1.tsv",
            false,
        )
    }

    fn annotate_brca1_vars(
        path_tsv: &str,
        report_worst_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        let txs = vec![
            String::from("NM_007294.4"),
            String::from("NM_007297.4"),
            String::from("NM_007298.3"),
            String::from("NM_007299.4"),
            String::from("NM_007300.4"),
        ];

        annotate_vars(path_tsv, &txs, report_worst_consequence_only)
    }

    fn annotate_vars(
        path_tsv: &str,
        txs: &[String],
        report_worst_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            Assembly::Grch37p10,
            Default::default(),
        ));
        let predictor = ConsequencePredictor::new(
            provider,
            Assembly::Grch37p10,
            ConfigBuilder::default()
                .report_worst_consequence_only(report_worst_consequence_only)
                .build()
                .unwrap(),
        );

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .comment(Some(b'#'))
            .from_reader(File::open(path_tsv).map(BufReader::new)?);

        // Read record with variant, transcript, and predicted consequences.
        //
        // We only have a limited set of transcripts for BRCA1 and this may not be the
        // same as the one in the file.  We thus limit ourselves to the transcripts
        // in `txs`.
        //
        // Also, the predicted ontology terms may not be the same as the ones in the
        // file.  We thus only check that we have a match in at least one term
        // of highest impact.
        let mut lineno = 0;
        for record in reader.deserialize() {
            lineno += 1;

            let record: Record = record?;
            // let mut printed = false;
            if txs.contains(&record.tx) {
                // "Parse" out the variant.
                let arr = record.var.split('-').collect::<Vec<_>>();
                // Predict consequences for the variant using Mehari.
                let anns = predictor
                    .predict(&VcfVariant {
                        chromosome: arr[0].to_string(),
                        position: arr[1].parse::<i32>()?,
                        reference: arr[2].to_string(),
                        alternative: arr[3].to_string(),
                    })?
                    .unwrap();
                // Now, for the overlapping transcripts, check that we have a match with the
                // consequences in the highest impact category.
                for ann in anns.iter().filter(|ann| ann.feature_id == record.tx) {
                    // We perform a comparison based on strings because we may not be able to parse out
                    // all consequences from the other tool.
                    let record_csqs = record
                        .csq
                        .split('&')
                        .map(|s| s.to_string())
                        .collect::<Vec<_>>();

                    let highest_impact = ann.consequences.first().unwrap().impact();
                    let mut expected_one_of = ann
                        .consequences
                        .iter()
                        .filter(|csq| csq.impact() == highest_impact)
                        .map(|csq| csq.to_string())
                        .collect::<Vec<_>>();

                    // Map effects a bit for VEP.
                    if path_tsv.ends_with(".vep.tsv")
                        && (expected_one_of.contains(&String::from("disruptive_inframe_deletion"))
                            || expected_one_of
                                .contains(&String::from("conservative_inframe_deletion")))
                    {
                        expected_one_of.push(String::from("inframe_deletion"));
                    }

                    // Try to find a direct match.
                    let found_one = record_csqs.iter().any(|csq| expected_one_of.contains(csq));
                    // It is common that the other tool predicts a frameshift variant while the actual prediction
                    // is stop_gained or stop_lost.  We thus also check for this case and allow it.
                    let found_one = found_one
                        || (record_csqs.contains(&String::from("frameshift_variant"))
                            && (expected_one_of.contains(&String::from("stop_gained")))
                            || expected_one_of.contains(&String::from("stop_lost")));
                    // VEP does not differentiate between disruptive and conservative inframe deletions and insertions.
                    let found_one = found_one
                        || (record_csqs.contains(&String::from("inframe_deletion"))
                            && (expected_one_of
                                .contains(&String::from("disruptive_inframe_deletion")))
                            || expected_one_of
                                .contains(&String::from("conservative_inframe_deletion")))
                        || (record_csqs.contains(&String::from("inframe_insertion"))
                            && (expected_one_of
                                .contains(&String::from("disruptive_inframe_insertion")))
                            || expected_one_of
                                .contains(&String::from("conservative_inframe_insertion")));
                    // NB: We cannot predict 5_prime_UTR_premature_start_codon_gain_variant yet. For now, we
                    // also accept 5_prime_UTR_variant.
                    let found_one = found_one
                        || ((expected_one_of.contains(&String::from("5_prime_UTR_exon_variant"))
                            || expected_one_of
                                .contains(&String::from("5_prime_UTR_intron_variant")))
                            && (record_csqs.contains(&String::from(
                                "5_prime_UTR_premature_start_codon_gain_variant",
                            ))));
                    // VEP predicts `splice_donor_5th_base_variant` rather than `splice_region_variant`.
                    // Same for `splice_donor_region_variant`.
                    let found_one = found_one
                        || (expected_one_of.contains(&String::from("splice_region_variant"))
                            && (record_csqs
                                .contains(&String::from("splice_donor_5th_base_variant"))
                                || record_csqs
                                    .contains(&String::from("splice_donor_region_variant"))));
                    // In the case of insertions at the end of an exon, VEP predicts `splice_region_variant`
                    // while we predict `splice_donor_variant`, same for start.
                    let found_one = found_one
                        || (expected_one_of.contains(&String::from("splice_donor_variant"))
                            || expected_one_of.contains(&String::from("splice_acceptor_variant")))
                            && (record_csqs.contains(&String::from("splice_region_variant")));
                    // VEP sometimes mispredicts disruptive inframe deletion as missense...
                    // cf. https://github.com/Ensembl/ensembl-vep/issues/1388
                    let found_one = found_one
                        || expected_one_of.contains(&String::from("disruptive_inframe_deletion"))
                            && (record_csqs.contains(&String::from("missense_variant")));
                    // VEP does not provide `exon_loss_variant`, so we also accept `inframe_deletion` and
                    // `splice_region_variant` (BRA1 test case).
                    let found_one = found_one
                        || expected_one_of.contains(&String::from("exon_loss_variant"))
                            && (record_csqs.contains(&String::from("inframe_deletion"))
                                || record_csqs.contains(&String::from("splice_region_variant")));
                    // On BRCA1, there is a case where VEP predicts `protein_altering_variant` rather than
                    // `disruptive_inframe_deletion`.  We accept this as well.
                    let found_one = found_one
                        || (expected_one_of.contains(&String::from("disruptive_inframe_deletion"))
                            || expected_one_of.contains(&String::from("inframe_indel")))
                            && (record_csqs.contains(&String::from("protein_altering_variant")));
                    // In the case of `GRCh37:17:41258543:T:TA`, the `hgvs` prediction is `c.-1_1insT` and
                    // `p.Met1?` which leads to `start_lost` while VEP predicts `5_prime_UTR_variant`.
                    // This may be a bug in `hgvs` and we don't change this for now.  We accept the call
                    // by VEP, of course.
                    let found_one = found_one
                        || expected_one_of.contains(&String::from("start_lost"))
                            && (record_csqs.contains(&String::from("5_prime_UTR_variant")));
                    // We have specialized {5,3}_prime_UTR_{exon,intron}_variant handling, while
                    // vep and snpEff do not
                    let found_one = found_one
                        || record_csqs.contains(&String::from("5_prime_UTR_variant"))
                            && (expected_one_of
                                .contains(&String::from("5_prime_UTR_exon_variant"))
                                || expected_one_of
                                    .contains(&String::from("5_prime_UTR_intron_variant")));
                    let found_one = found_one
                        || record_csqs.contains(&String::from("3_prime_UTR_variant"))
                            && (expected_one_of
                                .contains(&String::from("3_prime_UTR_exon_variant"))
                                || expected_one_of
                                    .contains(&String::from("3_prime_UTR_intron_variant")));
                    // an inframe_indel can be a missense_variant if it is an MNV (which we do not explicitly check here)
                    let found_one = found_one
                        || expected_one_of.contains(&String::from("inframe_indel"))
                            && (record_csqs.contains(&String::from("missense_variant")));
                    // inframe_indel also is a superclass of *_inframe_{deletion, insertion}
                    let found_one = found_one
                        || expected_one_of.contains(&String::from("inframe_indel"))
                            && [
                                "disruptive_inframe_deletion",
                                "conservative_inframe_deletion",
                                "inframe_deletion",
                                "disruptive_inframe_insertion",
                                "conservative_inframe_insertion",
                                "inframe_insertion",
                            ]
                            .iter()
                            .any(|c| record_csqs.contains(&String::from(*c)));
                    // SnpEff has a different interpretation of disruptive/conservative inframe deletions.
                    // We thus allow both.
                    let found_one = found_one
                        || expected_one_of.contains(&String::from("disruptive_inframe_deletion"))
                            && (record_csqs
                                .contains(&String::from("conservative_inframe_deletion")))
                        || expected_one_of.contains(&String::from("disruptive_inframe_insertion"))
                            && (record_csqs
                                .contains(&String::from("conservative_inframe_insertion")));
                    // SnpEff may not predict `splice_region_variant` for 5' UTR correctly, so we
                    // allow this.
                    let found_one = found_one
                        || expected_one_of.contains(&String::from("splice_region_variant"))
                            && (record_csqs.contains(&String::from("5_prime_UTR_variant")));
                    // SnpEff does not predict `splice_polypyrimidine_tract_variant`
                    let found_one = found_one
                        || expected_one_of
                            .contains(&String::from("splice_polypyrimidine_tract_variant"))
                            && (record_csqs.contains(&String::from("splice_region_variant"))
                                || record_csqs.contains(&String::from("intron_variant")));
                    // For `GRCh37:3:193366573:A:ATATTGCCTAGAATGAACT`, SnpEff predicts
                    // `stop_gained` while this rather is a intron variant.  We skip this variant.
                    let found_one = found_one
                        || record_csqs.contains(&String::from("stop_gained"))
                            && record.var == "3-193366573-A-ATATTGCCTAGAATGAACT";
                    // For `GRCh37:3:193409913:ATAAAT:A`, there appears to be a model error
                    // in SnpEff as it predicts `exon_loss`.  We skip this variant.
                    let found_one = found_one
                        || record_csqs.contains(&String::from("exon_loss_variant"))
                            && record.var == "3-193409913-ATAAAT-A";
                    // SnpEff may predict `pMet1.?` as `initiator_codon_variant` rather than `start_lost`.
                    let found_one = found_one
                        || expected_one_of.contains(&String::from("start_lost"))
                            && (record_csqs.contains(&String::from("initiator_codon_variant")));
                    // Similarly, SnpEff may predict `c.-1_1` as `start_retained` rather than `start_lost`.
                    let found_one = found_one
                        || expected_one_of.contains(&String::from("start_lost"))
                            && (record_csqs.contains(&String::from("start_retained_variant")));

                    assert!(
                        found_one,
                        "line no. {}, variant: {}, tx: {}, hgvs_c: {:?}, hgvs_p: {:?}, \
                        record_csqs: {:?}, expected_one_of: {:?}, ann_csqs: {:?}",
                        lineno,
                        record.var,
                        record.tx,
                        ann.hgvs_t.as_ref(),
                        ann.hgvs_p.as_ref(),
                        &record_csqs,
                        &expected_one_of,
                        &ann.consequences,
                    );
                }
            }
        }

        Ok(())
    }
}
