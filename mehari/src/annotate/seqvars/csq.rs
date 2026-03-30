//! Compute molecular consequence of variants.
use super::{
    ann::{Allele, AnnField, Consequence, FeatureBiotype, FeatureType, Pos, Rank, SoFeature},
    provider::Provider as MehariProvider,
};
use crate::annotate::cli::{ConsequenceBy, TranscriptSource};
use crate::annotate::seqvars::ann::{
    ANN_AA_SEQ_ALT, ANN_AA_SEQ_REF, ANN_TX_SEQ_ALT, ANN_TX_SEQ_REF, FeatureTag, GroupedAlleles,
};
use crate::annotate::seqvars::provider::PbsTranscriptExt;
use crate::errors::{GroupValidationError, SeqvarsError};
use crate::pbs::txs::{GenomeAlignment, Strand, Transcript, TranscriptBiotype, TranscriptTag};
use enumflags2::BitFlags;
use hgvs::mapper::altseq::AltSeqBuilder;
use hgvs::parser::{NoRef, ProteinEdit};
use hgvs::{
    data::interface::{Provider, TxForRegionRecord},
    mapper::{Error, assembly},
    parser::{
        Accession, CdsFrom, GenomeInterval, GenomeLocEdit, HgvsVariant, Mu, NaEdit, ProtLocEdit,
    },
};
use itertools::Itertools;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::{collections::HashMap, sync::Arc};

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

/// Configuration for consequence prediction.
#[derive(Debug, Clone, derive_builder::Builder)]
#[builder(pattern = "immutable")]
pub struct Config {
    /// The transcript source to use.
    #[builder(default = "TranscriptSource::Both")]
    pub transcript_source: TranscriptSource,

    /// Whether to report only the worst consequence for each picked transcript.
    #[builder(default)]
    pub report_most_severe_consequence_by: Option<ConsequenceBy>,

    /// Whether to keep intergenic variants.
    #[builder(default = "false")]
    pub keep_intergenic: bool,

    /// Whether to report splice variants in UTRs.
    #[builder(default = "false")]
    pub discard_utr_splice_variants: bool,

    /// Whether to normalize HGVS variants.
    #[builder(default = "true")]
    pub normalize: bool,

    /// Whether re-normalize genomic variants.
    #[builder(default = "true")]
    pub renormalize_g: bool,

    /// Whether to report VEP consequence terms.
    #[builder(default = "false")]
    pub vep_consequence_terms: bool,

    /// Whether to report cDNA sequence.
    #[builder(default)]
    pub report_cdna_sequence: SequenceReporting,

    /// Whether to report protein sequence.
    #[builder(default)]
    pub report_protein_sequence: SequenceReporting,

    /// Ordered list of extra columns registered by plugins/features.
    #[builder(default)]
    pub custom_columns: Vec<String>,
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
    pub(crate) provider: Arc<MehariProvider>,

    /// Assembly mapper for variant consequence prediction.
    #[derivative(Debug = "ignore")]
    mapper: assembly::Mapper,

    /// Configuration for the predictor.
    #[derivative(Debug = "ignore")]
    pub(crate) config: Config,
}

/// Padding to look for genes upstream/downstream.
pub const PADDING: i32 = 5_000;
/// Generally used alternative alignment method.
pub const ALT_ALN_METHOD: &str = "splign";

pub type Consequences = BitFlags<Consequence>;

#[derive(Debug, Clone)]
struct HgvsProjectionContext {
    #[allow(dead_code)]
    g: HgvsVariant,
    n: Option<HgvsVariant>,
    c: Option<HgvsVariant>,
    p: Option<HgvsVariant>,
}

impl HgvsProjectionContext {
    /// Check if the variant falls within the boundaries of the CDS (start to stop),
    /// including intronic parts.
    fn is_within_cds_bounds(&self) -> bool {
        if let Some(HgvsVariant::CdsVariant { loc_edit, .. }) = &self.c {
            let loc = loc_edit.loc.inner();
            let start = &loc.start;
            let end = &loc.end;

            let is_5_prime = start.cds_from == CdsFrom::Start
                && start.base < 0
                && end.cds_from == CdsFrom::Start
                && end.base < 0;

            let is_3_prime = start.cds_from == CdsFrom::End && end.cds_from == CdsFrom::End;

            !is_5_prime && !is_3_prime
        } else {
            false
        }
    }

    #[allow(dead_code)]
    /// Check if the variant is strictly within the coding sequence, i.e., not in UTRs or intronic.
    fn is_within_coding_sequence(&self) -> bool {
        if let Some(HgvsVariant::CdsVariant { loc_edit, .. }) = &self.c {
            let loc = loc_edit.loc.inner();

            if !self.is_within_cds_bounds() {
                return false;
            }

            let start_offset = loc.start.offset.unwrap_or(0);
            let end_offset = loc.end.offset.unwrap_or(0);

            start_offset == 0 && end_offset == 0
        } else {
            false
        }
    }
}

#[derive(Debug)]
struct TranscriptLocationContext {
    rank: Rank,
    distance: Option<i32>,
    is_exonic: bool,
    is_intronic: bool,
    is_upstream: bool,
    is_downstream: bool,
}

#[derive(Debug)]
struct ConsequenceContext {
    cds_consequences: Consequences,
    protein_consequences: Consequences,
    cdna_pos: Option<Pos>,
    cds_pos: Option<Pos>,
    protein_pos: Option<Pos>,
}

impl ConsequencePredictor {
    pub fn new(provider: Arc<MehariProvider>, config: Config) -> Self {
        tracing::info!("Building transcript interval trees ...");

        let reference_available = provider.reference_available();

        let mapper_config = assembly::Config {
            assembly: provider.assembly(),
            replace_reference: reference_available,
            strict_bounds: false,
            renormalize_g: reference_available && config.renormalize_g,
            genome_seq_available: reference_available,
            normalize: config.normalize,
            ..Default::default()
        };
        let mapper = assembly::Mapper::new(mapper_config, provider.clone());
        tracing::info!("... done building transcript interval trees");

        ConsequencePredictor {
            provider,
            mapper,
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
    pub fn predict(&self, var: &VcfVariant) -> Result<Option<Vec<AnnField>>, SeqvarsError> {
        // Normalize variant by stripping common prefix and suffix.
        let mut norm_var = self.normalize_variant(var);

        // TODO check for VCF specification version.
        // According to VCF specification (>=4.1), an alternative of "N" means REF=ALT
        // Prior to 4.1, it indicated a deletion.
        if norm_var.alternative == "N" {
            norm_var.alternative = norm_var.reference.clone();
        }

        // Obtain accession from chromosome name.
        let chrom_acc = self
            .provider
            .contig_manager
            .get_accession(&norm_var.chromosome);
        let chrom_acc = if let Some(chrom_acc) = chrom_acc {
            chrom_acc
        } else {
            tracing::warn!(
                "Could not determine chromosome accession for {:?}; giving up on annotation",
                &norm_var
            );
            return Ok(None);
        };

        // We follow hgvs conventions and therefore normalize input variants
        let var_g = Self::get_var_g(&norm_var, chrom_acc);
        let (var_g_fwd, var_g_rev) = if self.mapper.config.renormalize_g {
            let right = self
                .mapper
                .variant_mapper()
                .right_normalizer()?
                .normalize(&var_g)?;
            let left = self
                .mapper
                .variant_mapper()
                .left_normalizer()?
                .normalize(&var_g)?;
            (right, left)
        } else {
            (var_g.clone(), var_g.clone())
        };

        // Get all affected transcripts.
        let (var_start_fwd, var_end_fwd) = Self::get_var_start_end(&var_g_fwd);
        let (var_start_rev, var_end_rev) = Self::get_var_start_end(&var_g_rev);

        let qry_start = var_start_fwd.min(var_start_rev) - PADDING;
        let qry_end = var_end_fwd.max(var_end_rev) + PADDING;

        let txs = {
            let mut txs = self
                .provider
                .get_tx_for_region(chrom_acc, ALT_ALN_METHOD, qry_start, qry_end)
                .map_err(|e| SeqvarsError::Provider(e.to_string()))?;
            txs.sort_by(|a, b| a.tx_ac.cmp(&b.tx_ac));
            // Filter transcripts to the picked ones from the selected
            // transcript source.
            self.filter_picked_sourced_txs(txs)
        };

        // Handle case of no overlapping transcripts -> intergenic.
        if txs.is_empty() {
            let hgvs_g = format!("{}", &NoRef(&var_g.clone()));
            let hgvs_g = Some(hgvs_g.split(':').nth(1).unwrap().to_owned());

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
                feature_tags: vec![],
                rank: None,
                distance: None,
                strand: 0,
                hgvs_g,
                hgvs_n: None,
                hgvs_c: None,
                hgvs_p: None,
                cdna_pos: None,
                cds_pos: None,
                protein_pos: None,
                gene_symbol: "".to_string(),
                messages: None,
                custom_fields: BTreeMap::new(),
            }])));
        }

        // Compute annotations for all (picked) transcripts first, skipping `None`` results.
        let anns_all_txs = txs
            .into_iter()
            .map(|tx| {
                if tx.alt_strand == -1 {
                    self.build_ann_field(var, var_g_rev.clone(), tx, var_start_rev, var_end_rev)
                } else {
                    self.build_ann_field(var, var_g_fwd.clone(), tx, var_start_fwd, var_end_fwd)
                }
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        // Return all or worst annotation only.
        Ok(Some(self.filter_ann_fields(anns_all_txs)))
    }

    fn get_var_start_end(var_g: &HgvsVariant) -> (i32, i32) {
        match &var_g {
            HgvsVariant::GenomeVariant { loc_edit, .. } => {
                let loc = loc_edit.loc.inner();
                let edit = loc_edit.edit.inner();
                let start = loc
                    .start
                    .map(|s| s - 1)
                    .expect("Failed to get start position");
                let end = loc.end.expect("Failed to get end position");
                // In insertion / duplication cases, range end is exclusive.
                // See https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/
                let end = if edit.is_ins() || edit.is_dup() {
                    end - 1
                } else {
                    end
                };
                (start, end)
            }
            _ => unreachable!(),
        }
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
        let ann_fields = if !self.config.keep_intergenic {
            ann_fields
                .into_iter()
                .filter(|field| field.consequences != [Consequence::IntergenicVariant])
                .collect()
        } else {
            ann_fields
        };

        /// Return sort order for ANN biotype, gives priority to ManeSelect and ManePlusClinical.
        fn tag_order(tags: &[FeatureTag]) -> i32 {
            if tags.contains(&FeatureTag::ManeSelect) {
                0
            } else if tags.contains(&FeatureTag::ManePlusClinical) {
                1
            } else {
                2
            }
        }

        /// Extract the first (i.e. most severe) consequence per group.
        fn first_csq_per_group(grouping: HashMap<String, Vec<AnnField>>) -> Vec<AnnField> {
            // Sort by `Consequence` and `FeatureBiotype`
            //
            // This uses the invariant that the consequences in the ANN fields are sorted already
            // and there is at least one consequence.
            grouping
                .into_values()
                .map(|mut anns| {
                    anns.sort_by_key(|ann| (ann.consequences[0], tag_order(&ann.feature_tags)));
                    anns.into_iter().next().unwrap()
                })
                .collect()
        }

        /// Group ANN fields by a key function.
        fn group_annotations_by<F: Fn(&AnnField) -> String>(
            ann_fields: Vec<AnnField>,
            key_fn: F,
        ) -> HashMap<String, Vec<AnnField>> {
            ann_fields.into_iter().into_group_map_by(key_fn)
        }

        match self.config.report_most_severe_consequence_by {
            // Short-circuit if to report all transcript results.
            None => ann_fields,
            Some(group) => {
                let key = match group {
                    ConsequenceBy::Gene => |ann: &AnnField| ann.gene_id.clone(),
                    ConsequenceBy::Transcript => |ann: &AnnField| ann.feature_id.clone(),
                    ConsequenceBy::Allele => |ann: &AnnField| ann.allele.to_string(),
                };
                first_csq_per_group(group_annotations_by(ann_fields, key))
            }
        }
    }

    fn determine_transcript_context(
        &self,
        alignment: &GenomeAlignment,
        strand: Strand,
        var_g: &HgvsVariant,
        var_start: i32,
        var_end: i32,
    ) -> (TranscriptLocationContext, Consequences, i32) {
        let mut consequences = Consequences::empty();
        let mut rank = Rank::default();
        let mut is_exonic = false;
        let mut is_intronic = false;
        let mut distance: Option<i32> = None;
        let mut tx_len = 0;

        let var_overlaps =
            |start: i32, end: i32| -> bool { overlaps(var_start, var_end, start, end) };

        let cds_start = alignment.cds_start.unwrap_or(-1);
        let cds_end = alignment.cds_end.unwrap_or(-1);

        // Find first exon that overlaps with variant or intron that contains the variant.
        //
        // Note that exons are stored in genome position order.
        let mut prev_end = None;
        let mut min_start = None;
        let mut max_end = None;

        for exon_alignment in &alignment.exons {
            tx_len += exon_alignment.alt_end_i - exon_alignment.alt_start_i;

            let exon_start = exon_alignment.alt_start_i;
            let exon_end = exon_alignment.alt_end_i;

            let is_utr = exon_start < cds_start && !var_overlaps(cds_start, cds_end)
                || exon_end > cds_end && !var_overlaps(cds_start, cds_end);

            // Check the cases where the variant overlaps with the exon or is contained within an
            // intron.
            if var_overlaps(exon_start, exon_end) {
                rank = Rank {
                    ord: exon_alignment.ord + 1,
                    total: alignment.exons.len() as i32,
                };
                is_exonic = true;
                distance = Some(0);
                consequences |= Self::analyze_exonic_variant(
                    strand, var_start, var_end, exon_start, exon_end, &rank, is_utr,
                );
            } else if let Some(intron_start) = prev_end
                && var_start >= intron_start
                && var_end <= exon_end
                && !is_exonic
            {
                rank = Rank {
                    ord: exon_alignment.ord + 1,
                    total: alignment.exons.len() as i32 - 1,
                };
                is_intronic = true;

                // We compute the "distance" with "+1", the first base of the
                // intron is "+1", the last one is "-1".
                let dist_start: i32 = var_start + 1 - intron_start;
                let dist_end: i32 = -(exon_start + 1 - var_end);
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

            if let Some(intron_start) = prev_end {
                consequences |= Self::analyze_intronic_variant(
                    var_g,
                    alignment,
                    strand,
                    var_start,
                    var_end,
                    intron_start,
                    exon_start,
                    is_utr,
                );
            }

            min_start = Some(std::cmp::min(min_start.unwrap_or(exon_start), exon_start));
            max_end = Some(std::cmp::max(max_end.unwrap_or(exon_end), exon_end));
            prev_end = Some(exon_end);
        }

        let min_start = min_start.expect("must have seen exon");
        let max_end = max_end.expect("must have seen exon");
        let is_upstream = var_end <= min_start;
        let is_downstream = var_start >= max_end;

        if !is_exonic && !is_intronic {
            if is_upstream {
                let val = -(min_start + 1 - var_end);
                if val.abs() <= PADDING {
                    consequences |= match strand {
                        Strand::Plus => Consequence::UpstreamGeneVariant,
                        Strand::Minus => Consequence::DownstreamGeneVariant,
                        _ => unreachable!("invalid strand: {}", alignment.strand),
                    };
                }
                if distance.is_none() {
                    distance = Some(val);
                }
            } else if is_downstream {
                let val = var_start + 1 - max_end;
                if val.abs() <= PADDING {
                    consequences |= match strand {
                        Strand::Plus => Consequence::DownstreamGeneVariant,
                        Strand::Minus => Consequence::UpstreamGeneVariant,
                        _ => unreachable!("invalid strand: {}", alignment.strand),
                    };
                }
                if distance.is_none() {
                    distance = Some(val);
                }
            }
        }

        (
            TranscriptLocationContext {
                rank,
                distance,
                is_exonic,
                is_intronic,
                is_upstream,
                is_downstream,
            },
            consequences,
            tx_len,
        )
    }

    fn project_hgvs(
        &self,
        var_g: &HgvsVariant,
        tx: &Transcript,
        transcript_biotype: TranscriptBiotype,
    ) -> Result<HgvsProjectionContext, SeqvarsError> {
        let mut projection = HgvsProjectionContext {
            g: var_g.clone(),
            n: None,
            c: None,
            p: None,
        };

        projection.n = self.mapper.g_to_n(var_g, &tx.id).map_or_else(
            |e| match e {
                Error::NonAdjacentExons(_, _, _, _) => {
                    tracing::warn!("{}, {}: NonAdjacentExons, skipping", &tx.id, var_g);
                    Ok(None)
                }
                _ => Err(SeqvarsError::from(e)),
            },
            |v| Ok(Some(v)),
        )?;

        if let Some(var_n) = &projection.n {
            projection.c = match transcript_biotype {
                TranscriptBiotype::Coding => self
                    .mapper
                    .n_to_c(var_n)
                    .map(Some)
                    .map_err(|e| SeqvarsError::HgvsProjection(e.to_string()))?,
                TranscriptBiotype::NonCoding => Some(var_n.clone()),
                _ => None,
            };

            if let Some(var_c) = &projection.c
                && transcript_biotype == TranscriptBiotype::Coding
            {
                projection.p = self.safe_project_c_to_p(var_c)?;
            }
        }

        Ok(projection)
    }

    fn analyze_transcript_consequences(
        &self,
        projection: &HgvsProjectionContext,
        tx: &Transcript,
        tx_record: &TxForRegionRecord,
        transcript_location: &TranscriptLocationContext,
        tx_len: i32,
        transcript_biotype: TranscriptBiotype,
    ) -> Result<ConsequenceContext, SeqvarsError> {
        let mut context = ConsequenceContext {
            cds_consequences: Consequences::empty(),
            protein_consequences: Consequences::empty(),
            cdna_pos: None,
            cds_pos: None,
            protein_pos: None,
        };

        if let Some(var_n) = &projection.n {
            context.cdna_pos = transcript_location.is_exonic.then_some(match var_n {
                HgvsVariant::TxVariant { loc_edit, .. } => Pos {
                    ord: loc_edit.loc.inner().start.base,
                    total: Some(tx_len),
                },
                _ => panic!("Invalid tx position: {:?}", var_n),
            });
        }

        if let Some(var_c) = &projection.c
            && transcript_biotype == TranscriptBiotype::Coding
        {
            // If there's no stop codon, we can't compute the cds_len.
            // We can, however, still analyze any consequences before the (missing) stop codon.
            let cds_len = tx.stop_codon.map(|stop| stop - tx.start_codon.unwrap());
            context.cds_pos = transcript_location.is_exonic.then_some(match var_c {
                HgvsVariant::CdsVariant { loc_edit, .. } => Pos {
                    ord: loc_edit.loc.inner().start.base,
                    total: cds_len,
                },
                _ => panic!("Invalid CDS position: {:?}", var_c),
            });

            let conservative = is_conservative_cds_variant(var_c);
            let incomplete_3p = tx.is_incomplete_3p();
            let available_cds_len = tx.available_cds_len(tx_len);
            context.cds_consequences = Self::analyze_cds_variant(
                var_c,
                transcript_location.is_exonic,
                conservative,
                incomplete_3p,
                available_cds_len,
            );

            if let Some(var_p) = &projection.p {
                let prot_len = cds_len
                    .expect("cds_len cannot be None if hgvs.p projection has been successful")
                    / 3;
                context.protein_pos = match var_p {
                    HgvsVariant::ProtVariant { loc_edit, .. } => match loc_edit {
                        ProtLocEdit::Ordinary { loc, .. } => Some(Pos {
                            ord: loc.inner().start.number,
                            total: Some(prot_len),
                        }),
                        _ => None,
                    },
                    _ => panic!("Not a protein position: {:?}", var_p),
                };

                context.protein_consequences = self.analyze_protein_variant(
                    var_c,
                    var_p,
                    &context.protein_pos,
                    conservative,
                    &tx_record.tx_ac,
                    incomplete_3p,
                );
            }
        }

        Ok(context)
    }

    fn build_ann_field(
        &self,
        orig_var: &VcfVariant,
        var_g: HgvsVariant,
        tx_record: TxForRegionRecord,
        var_start: i32,
        var_end: i32,
    ) -> Result<Option<AnnField>, SeqvarsError> {
        let tx = match self.provider.get_tx(&tx_record.tx_ac) {
            Some(tx) => {
                if TranscriptBiotype::try_from(tx.biotype).expect("invalid tx biotype")
                    == TranscriptBiotype::Coding
                    && tx.start_codon.is_none()
                {
                    tracing::debug!(
                        "Skipping transcript {} because it is coding but has no known start codon",
                        &tx_record.tx_ac
                    );
                    return Ok(None);
                }
                tx
            }
            None => {
                tracing::warn!(
                    "Requested transcript accession {}, got None (potentially filtered)",
                    &tx_record.tx_ac
                );
                return Ok(None);
            }
        };

        assert_eq!(
            tx.genome_alignments.len(),
            1,
            "At this point, only one genome alignment is expected"
        );

        let alignment = tx.genome_alignments.first().unwrap();
        let strand = Strand::try_from(alignment.strand).expect("invalid strand");
        let transcript_biotype =
            TranscriptBiotype::try_from(tx.biotype).expect("invalid transcript biotype");

        let (transcript_location, transcript_consequences, tx_len) =
            self.determine_transcript_context(alignment, strand, &var_g, var_start, var_end);

        let mut consequences = transcript_consequences;

        if transcript_location.is_exonic {
            if transcript_biotype == TranscriptBiotype::NonCoding {
                consequences |= Consequence::NonCodingTranscriptExonVariant;
            }
        } else if transcript_location.is_intronic {
            if transcript_biotype == TranscriptBiotype::NonCoding {
                consequences |= Consequence::NonCodingTranscriptIntronVariant;
            } else {
                consequences |= Consequence::CodingTranscriptIntronVariant;
            }
        }

        let (rank, projection, cdna_pos, cds_pos, protein_pos) = if !transcript_location.is_upstream
            && !transcript_location.is_downstream
        {
            let projection = self.project_hgvs(&var_g, tx, transcript_biotype)?;
            if projection.n.is_none() {
                return Ok(None); // Early exit if g->n projection failed.
            }

            let consequence_ctx = self.analyze_transcript_consequences(
                &projection,
                tx,
                &tx_record,
                &transcript_location,
                tx_len,
                transcript_biotype,
            )?;

            consequences |= consequence_ctx.cds_consequences | consequence_ctx.protein_consequences;

            self.consequences_fix_special_cases(
                &mut consequences,
                // exon_alignment_consequences, // TODO include these as well
                consequence_ctx.cds_consequences,
                consequence_ctx.protein_consequences,
                &projection,
            );

            (
                Some(transcript_location.rank),
                Some(projection),
                consequence_ctx.cdna_pos,
                consequence_ctx.cds_pos,
                consequence_ctx.protein_pos,
            )
        } else {
            (None, None, None, None, None)
        };

        let mut custom_fields = BTreeMap::new();

        let c_ref = self.config.report_cdna_sequence.includes_ref();
        let c_alt = self.config.report_cdna_sequence.includes_alt();
        let p_ref = self.config.report_protein_sequence.includes_ref();
        let p_alt = self.config.report_protein_sequence.includes_alt();

        if (c_ref || c_alt || p_ref || p_alt)
            && let Some(var_c) = projection.as_ref().and_then(|p| p.c.as_ref())
            && let Ok(ref_data) = hgvs::mapper::altseq::ref_transcript_data_cached(
                self.provider.clone(),
                &tx.id,
                None,
            )
        {
            if c_ref {
                custom_fields.insert(
                    ANN_TX_SEQ_REF.into(),
                    Some(ref_data.transcript_sequence.clone()),
                );
            }
            if p_ref {
                custom_fields.insert(ANN_AA_SEQ_REF.into(), Some(ref_data.aa_sequence.clone()));
            }

            if (c_alt || p_alt)
                && matches!(var_c, HgvsVariant::CdsVariant { .. })
                && let Ok(alt_data_vec) = AltSeqBuilder::new(var_c.clone(), ref_data).build_altseq()
                && let Some(alt_data) = alt_data_vec.into_iter().next()
            {
                if c_alt {
                    custom_fields.insert(ANN_TX_SEQ_ALT.into(), Some(alt_data.transcript_sequence));
                }
                if p_alt {
                    custom_fields.insert(ANN_AA_SEQ_ALT.into(), Some(alt_data.aa_sequence));
                }
            }
        }

        let hgvs_n = projection.as_ref().and_then(|p| p.n.as_ref()).map(|var_n| {
            format!("{}", &NoRef(var_n))
                .split(':')
                .nth(1)
                .unwrap()
                .to_owned()
        });

        let hgvs_c = projection.as_ref().and_then(|p| p.c.as_ref()).map(|var_c| {
            format!("{}", &NoRef(var_c))
                .split(':')
                .nth(1)
                .unwrap()
                .to_owned()
        });

        let hgvs_p = projection
            .as_ref()
            .and_then(|p| p.p.as_ref())
            .map(|var_p| format!("{}", var_p).split(':').nth(1).unwrap().to_owned());

        let feature_biotype = vec![match transcript_biotype {
            TranscriptBiotype::Coding => FeatureBiotype::Coding,
            TranscriptBiotype::NonCoding => FeatureBiotype::Noncoding,
            _ => unreachable!("invalid biotype: {:?}", transcript_biotype),
        }];
        let feature_tags = tx
            .tags
            .iter()
            .map(|tag| TranscriptTag::try_from(*tag).expect("invalid transcript tag"))
            .filter(|tag| !matches!(tag, TranscriptTag::EnsemblGraft))
            .filter_map(|transcript_tag| {
                let tag = FeatureTag::from(transcript_tag);
                if !matches!(tag, FeatureTag::Other(_)) {
                    Some(tag)
                } else {
                    None
                }
            })
            .collect_vec();

        if consequences.is_empty() {
            tracing::error!(
                "No consequences for {:?} on {} (hgvs_n={}, hgvs_c={}, hgvs_p={}) - adding `gene_variant`; \
                most likely the transcript has multiple stop codons and the variant \
                lies behind the first.",
                orig_var,
                &tx_record.tx_ac,
                hgvs_n.as_deref().unwrap_or("None"),
                hgvs_c.as_deref().unwrap_or("None"),
                hgvs_p.as_deref().unwrap_or("None")
            );
            consequences |= Consequence::GeneVariant;
        }

        if self.config.vep_consequence_terms {
            self.adjust_vep_terms(&mut consequences, projection.as_ref());
        }

        let consequences = consequences.iter().collect_vec();
        let putative_impact = (*consequences.first().unwrap()).into();

        let strand = match strand {
            Strand::Unknown => 0,
            Strand::Plus => 1,
            Strand::Minus => -1,
        };

        let hgvs_g = Some(
            format!("{}", &NoRef(&var_g))
                .split(':')
                .nth(1)
                .unwrap()
                .to_owned(),
        );

        Ok(Some(AnnField {
            allele: Allele::Alt {
                alternative: orig_var.alternative.clone(),
            },
            consequences,
            putative_impact,
            gene_symbol: tx.gene_symbol.clone(),
            gene_id: tx.gene_id.clone(),
            feature_type: FeatureType::SoTerm {
                term: SoFeature::Transcript,
            },
            feature_id: tx.id.clone(),
            feature_biotype,
            feature_tags,
            rank,
            hgvs_g,
            hgvs_n,
            hgvs_c,
            hgvs_p,
            cdna_pos,
            cds_pos,
            protein_pos,
            strand,
            distance: transcript_location.distance,
            messages: None,
            custom_fields,
        }))
    }

    fn adjust_vep_terms(
        &self,
        consequences: &mut Consequences,
        projection_context: Option<&HgvsProjectionContext>,
    ) {
        use crate::annotate::seqvars::ann::Consequence::*;

        // vep reports the umbrella intron variant term.
        if consequences.contains(CodingTranscriptIntronVariant) {
            consequences.remove(CodingTranscriptIntronVariant);
            consequences.insert(IntronVariant);
        }
        if consequences.contains(NonCodingTranscriptIntronVariant) {
            consequences.remove(NonCodingTranscriptIntronVariant);
            consequences.insert(IntronVariant);
            consequences.insert(NonCodingTranscriptVariant);
        }
        if consequences.contains(FivePrimeUtrIntronVariant) {
            consequences.remove(FivePrimeUtrIntronVariant);
            consequences.insert(IntronVariant);
        }
        if consequences.contains(ThreePrimeUtrIntronVariant) {
            consequences.remove(ThreePrimeUtrIntronVariant);
            consequences.insert(IntronVariant);
        }

        if consequences.contains(FivePrimeUtrExonVariant) {
            consequences.remove(FivePrimeUtrExonVariant);
            consequences.insert(FivePrimeUtrVariant);
        }
        if consequences.contains(ThreePrimeUtrExonVariant) {
            consequences.remove(ThreePrimeUtrExonVariant);
            consequences.insert(ThreePrimeUtrVariant);
        }

        if consequences.contains(SelenocysteineGain | SelenocysteineLoss) {
            consequences.remove(SelenocysteineGain | SelenocysteineLoss);
            consequences.insert(MissenseVariant);
        }

        // VEP marks some intronic variants as coding_sequence_variant,
        // which is wrong, but we will mimic here for compatibility
        let within_cds_bounds = projection_context
            .map(|c| c.is_within_cds_bounds())
            .unwrap_or(false);
        if consequences.contains(ExonicSpliceRegionVariant) {
            consequences.remove(ExonicSpliceRegionVariant);
            consequences.insert(SpliceRegionVariant);
            if within_cds_bounds {
                consequences.insert(CodingSequenceVariant);
            }
        }

        if consequences.intersects(FrameshiftElongation | FrameshiftTruncation) {
            consequences.remove(FrameshiftElongation | FrameshiftTruncation);
            consequences.insert(FrameshiftVariant);
        }

        let inframe_specifics = DisruptiveInframeDeletion
            | DisruptiveInframeInsertion
            | ConservativeInframeDeletion
            | ConservativeInframeInsertion;

        if consequences.intersects(inframe_specifics) {
            if consequences.intersects(DisruptiveInframeDeletion | ConservativeInframeDeletion) {
                consequences.insert(InframeDeletion);
            }
            if consequences.intersects(DisruptiveInframeInsertion | ConservativeInframeInsertion) {
                consequences.insert(InframeInsertion);
            }
            consequences.remove(inframe_specifics);
        }

        let is_inframe = consequences.intersects(InframeDeletion | InframeInsertion);
        let essential_splice = SpliceDonorVariant | SpliceAcceptorVariant;

        if is_inframe {
            if consequences.intersects(essential_splice) {
                consequences.remove(InframeDeletion | InframeInsertion);
                if within_cds_bounds {
                    consequences.insert(CodingSequenceVariant);
                }
            } else if consequences.intersects(
                SpliceRegionVariant
                    | SpliceDonorFifthBaseVariant
                    | SpliceDonorRegionVariant
                    | SplicePolypyrimidineTractVariant,
            ) && within_cds_bounds
            {
                consequences.insert(CodingSequenceVariant);
            }
        }

        if consequences.contains(ExonLossVariant)
            && *consequences != Into::<Consequences>::into(ExonLossVariant)
        {
            consequences.remove(ExonLossVariant);
        }

        if consequences.contains(FrameshiftVariant) {
            consequences.remove(StopGained | StopLost);
        }

        let suppress_splice_region = SpliceDonorVariant
            | SpliceAcceptorVariant
            | SpliceDonorFifthBaseVariant
            | SpliceDonorRegionVariant;

        if consequences.intersects(suppress_splice_region) {
            consequences.remove(SpliceRegionVariant);
        }

        if consequences.contains(SpliceDonorFifthBaseVariant) {
            consequences.remove(SpliceDonorRegionVariant);
        }

        if consequences.contains(SpliceAcceptorVariant) {
            consequences.remove(SplicePolypyrimidineTractVariant);
        }

        let suppress_intron =
            SpliceDonorVariant | SpliceAcceptorVariant | SpliceDonorFifthBaseVariant;

        if consequences.intersects(suppress_intron) {
            consequences.remove(IntronVariant);
        }
    }

    #[allow(clippy::too_many_arguments, unused_variables)]
    fn consequences_fix_special_cases(
        &self,
        consequences: &mut Consequences,
        consequences_cds: Consequences,
        consequences_protein: Consequences,
        projection: &HgvsProjectionContext,
    ) {
        // If we have a transcript_ablation, we can remove all other consequences
        if consequences.contains(Consequence::TranscriptAblation) {
            *consequences = Consequence::TranscriptAblation.into();
            return;
        }

        if let Some(var_c) = projection.c.as_ref() {
            // Do not report splice variants in UTRs.
            if self.config.discard_utr_splice_variants {
                let splice_variants = Consequence::ExonicSpliceRegionVariant
                    | Consequence::SpliceDonorVariant
                    | Consequence::SpliceAcceptorVariant
                    | Consequence::SpliceRegionVariant
                    | Consequence::SpliceDonorFifthBaseVariant
                    | Consequence::SpliceDonorRegionVariant
                    | Consequence::SplicePolypyrimidineTractVariant;
                let utr_intron_variants = Consequence::FivePrimeUtrIntronVariant
                    | Consequence::ThreePrimeUtrIntronVariant;
                let utr_exon_variants =
                    Consequence::FivePrimeUtrExonVariant | Consequence::ThreePrimeUtrExonVariant;
                let is_utr = match var_c {
                    HgvsVariant::CdsVariant { loc_edit, .. } => {
                        let loc = loc_edit.loc.inner();
                        loc.start.base < 0
                            && loc.end.base < 0
                            && loc.start.cds_from == CdsFrom::Start
                            && loc.end.cds_from == CdsFrom::Start
                            || loc.start.base > 0
                                && loc.end.base > 0
                                && loc.start.cds_from == CdsFrom::End
                                && loc.end.cds_from == CdsFrom::End
                    }
                    _ => false, // Not a CDS variant, can't be UTR in this context
                };
                if is_utr
                    && consequences.intersects(splice_variants)
                    && consequences.intersects(utr_intron_variants | utr_exon_variants)
                {
                    *consequences &= !splice_variants;
                }
            }
        }

        // vep simply reports the frameshift on the protein level, irrespective of the
        // actual outcome
        // i.e., this depends on if you want to have the mechanism or the outcome described
        if !self.config.vep_consequence_terms {
            // If a frameshift/ins/del was predicted on the CDS level,
            // but any relevant consequence (i.e. not just GeneVariant) was produced on the protein level,
            // then it is likely that the frameshift induced a more specific consequence.
            let check_cds_csqs: Consequences = Consequence::FrameshiftVariant.into();
            let checked = consequences_cds & check_cds_csqs;
            if checked != Consequences::empty()
                // if the protein consequence is not effectively empty, we remove the CDS frameshift consequence
                && !(consequences_protein.eq(&Consequence::GeneVariant) || consequences_protein.is_empty())
                // if the protein consequence also includes a frameshift, then we keep it
                && !consequences_protein.intersects(
                Consequence::FrameshiftElongation
                    | Consequence::FrameshiftTruncation
                    | Consequence::FrameshiftVariant,
            ) {
                *consequences &= !checked;
            }
        }

        // In some cases, we predict a stop lost based on the cds variant
        // but the protein translation does not confirm this.
        //
        // e.g.:
        // 20:35511609:CAAGCCGCCTCCAGGTAGCAGCCACAGCCAGGAGCACACAGACAGAAGACTGTGTCATGGGTCATGGCCCCTCCGCACACCTACAGGTTTGCCAAAGGAA:C
        if consequences_cds.contains(Consequence::StopLost)
            && !consequences_protein
                .intersects(Consequence::StopLost | Consequence::ProteinAlteringVariant)
            && projection.p.as_ref().is_some_and(|p| {
                !matches!(
                    p,
                    HgvsVariant::ProtVariant {
                        loc_edit: ProtLocEdit::NoProteinUncertain,
                        ..
                    }
                )
            })
        {
            *consequences &= !Consequence::StopLost;
        }

        // Similarly, for the start lost case
        //
        // e.g.:
        // 13:32316456:TA:T
        // (This case just shortens a poly-A from which the start codon starts)
        if consequences_cds.contains(Consequence::StartLost)
            && !consequences_protein.contains(Consequence::StartLost)
            && !consequences_protein.is_empty()
        {
            *consequences &= !Consequence::StartLost;
        }

        if consequences.contains(Consequence::StartLost)
            && let (
                Some(HgvsVariant::TxVariant {
                    loc_edit: n_loc_edit,
                    accession,
                    ..
                }),
                Some(HgvsVariant::CdsVariant {
                    loc_edit: c_loc_edit,
                    ..
                }),
            ) = (&projection.n, &projection.c)
        {
            let n_loc = n_loc_edit.loc.inner();
            let c_edit = c_loc_edit.edit.inner();
            let c_loc = c_loc_edit.loc.inner();

            // If edit occurs within the first 3 bases of the CDS,
            let (start, end) = (c_loc.start.base, c_loc.end.base);
            if start >= 1
                && end <= 3
                && c_loc.start.cds_from == CdsFrom::Start
                && c_loc.end.cds_from == CdsFrom::Start
            {
                // … then we need to check whether this is a start lost or a start retained.
                // To that end, extract the first 3 bases plus/minus 3 bases …
                if let Ok(first_codon_pm1) = self.provider.get_seq_part(
                    &accession.value,
                    Some(
                        usize::try_from(n_loc.start.base - start + 1)
                            .unwrap()
                            .saturating_sub(4),
                    ),
                    Some(usize::try_from(n_loc.end.base - start + 1).unwrap() + 5),
                ) {
                    // … and introduce the change into the sequence.
                    let mut first_codon = first_codon_pm1.clone();
                    let (start, end) = (start as usize, end as usize);
                    let start_retained = match c_edit {
                        NaEdit::DelRef { .. } => {
                            first_codon.replace_range(3 + start - 1..=3 + end - 1, "");
                            // If the first codon is still a start codon, then it is a start retained.
                            first_codon[2..5].contains("ATG")
                        }
                        // TODO: handle other cases
                        _ => false,
                    };
                    if start_retained {
                        tracing::trace!("Fixing StartLost → StartRetained for {:?}", &projection,);
                        *consequences &= !Consequence::StartLost;
                        *consequences |= Consequence::StartRetainedVariant;
                    }
                }
            }
        }

        if let Some(HgvsVariant::CdsVariant { loc_edit, .. }) = projection.c.as_ref() {
            let loc = loc_edit.loc.inner();
            let start_base = loc.start.base;
            let start_cds_from = loc.start.cds_from;
            let end_base = loc.end.base;
            let end_cds_from = loc.end.cds_from;

            let starts_left_of_start = start_cds_from == CdsFrom::Start && start_base < 0;
            let ends_left_of_start = end_cds_from == CdsFrom::Start && end_base < 0;

            if consequences.contains(Consequence::ExonLossVariant)
                && starts_left_of_start
                && ends_left_of_start
            {
                *consequences &= !Consequence::ExonLossVariant;
            }
        }

        if let Some(HgvsVariant::ProtVariant {
            loc_edit: ProtLocEdit::Unknown,
            ..
        }) = projection.p.as_ref()
            && consequences.is_empty()
            && projection.is_within_coding_sequence()
        {
            *consequences |= Consequence::CodingSequenceVariant;
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn analyze_exonic_variant(
        strand: Strand,
        var_start: i32,
        var_end: i32,
        exon_start: i32,
        exon_end: i32,
        rank: &Rank,
        _is_utr: bool,
    ) -> Consequences {
        let mut consequences: Consequences = Consequences::empty();

        let var_overlaps =
            |start: i32, end: i32| -> bool { overlaps(var_start, var_end, start, end) };

        // Check the cases where the variant overlaps with whole exon.
        if var_start <= exon_start && var_end >= exon_end {
            // FIXME: this is not true if the var_c variant is effectively pre start completely
            // we address that in fix_special_cases
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

        // Check splice region variants
        if var_overlaps(exon_end - 3, exon_end) {
            if strand == Strand::Plus {
                if !rank.is_last() {
                    consequences |= Consequence::ExonicSpliceRegionVariant;
                }
            } else {
                // alignment.strand == Strand::Minus
                if !rank.is_first() {
                    consequences |= Consequence::ExonicSpliceRegionVariant;
                }
            }
        }
        if var_overlaps(exon_start, exon_start + 3) {
            if strand == Strand::Plus {
                if !rank.is_first() {
                    consequences |= Consequence::ExonicSpliceRegionVariant;
                }
            } else {
                // alignment.strand == Strand::Minus
                if !rank.is_last() {
                    consequences |= Consequence::ExonicSpliceRegionVariant;
                }
            }
        }
        consequences
    }

    #[allow(clippy::too_many_arguments)]
    fn analyze_intronic_variant(
        var_g: &HgvsVariant,
        alignment: &GenomeAlignment,
        strand: Strand,
        var_start: i32,
        var_end: i32,
        intron_start: i32,
        intron_end: i32,
        _is_utr: bool,
    ) -> Consequences {
        let mut consequences: Consequences = Consequences::empty();

        let var_overlaps =
            |start: i32, end: i32| -> bool { overlaps(var_start, var_end, start, end) };

        // Insertion ranges are end-exclusive, so subtract 1 from the end, where applicable.
        let ins_shift = Self::ins_shift(var_g);

        // Check the cases where the variant overlaps with the splice acceptor/donor site.
        if var_overlaps(intron_start - ins_shift, intron_start + 2) {
            // Left side, is acceptor/donor depending on transcript's strand.
            match strand {
                Strand::Plus => {
                    consequences |= Consequence::SpliceDonorVariant;
                }
                Strand::Minus => {
                    consequences |= Consequence::SpliceAcceptorVariant;
                }
                _ => unreachable!("invalid strand: {}", alignment.strand),
            }
        }

        // Check the case where the variant overlaps with the splice donor site.
        if var_overlaps(intron_end - 2, intron_end + ins_shift) {
            // Left side, is acceptor/donor depending on transcript's strand.
            match strand {
                Strand::Plus => {
                    consequences |= Consequence::SpliceAcceptorVariant;
                }
                Strand::Minus => {
                    consequences |= Consequence::SpliceDonorVariant;
                }
                _ => unreachable!("invalid strand: {}", alignment.strand),
            }
        }

        // Check the case where the variant overlaps with the splice region (1-3 bases in exon
        // or 3-8 bases in intron).
        // n.b. the 1-3 bases in exon check is already done within `analyze_exonic_variant`.
        // We have to check all cases independently and not with `else`
        // because the variant may be larger.
        if var_overlaps(intron_start + 2, intron_start + 8)
            || var_overlaps(intron_end - 8, intron_end - 2)
        {
            consequences |= Consequence::SpliceRegionVariant;
        }

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
        // Note that this is two bases short of the intronic part of splice_region_variant
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

        consequences
    }

    fn ins_shift(var_g: &HgvsVariant) -> i32 {
        // For insertions, we need to consider the case of the insertion being right at
        // the exon/intron junction.  We can express this with a shift of 1 for using
        // "< / >" X +/- shift and meaning <= / >= X.

        match var_g {
            HgvsVariant::GenomeVariant {
                loc_edit: GenomeLocEdit { edit, .. },
                ..
            } => {
                let edit = edit.inner();
                if edit.is_ins() || edit.is_dup() { 1 } else { 0 }
            }
            _ => unreachable!(),
        }
    }

    fn analyze_cds_variant(
        var_c: &HgvsVariant,
        is_exonic: bool,
        conservative: bool,
        incomplete_3p: bool,
        available_cds_len: Option<i32>,
    ) -> Consequences {
        let mut consequences: Consequences = Consequences::empty();

        if let HgvsVariant::CdsVariant { loc_edit, .. } = &var_c {
            // Handle the cases where the variant touches the start or stop codon based on `var_c`
            // coordinates.  The cases where the start/stop codon is touched by the variant
            // directly is handled above based on the `var_p` prediction.
            let loc = loc_edit.loc.inner();
            let edit = loc_edit.edit.inner();
            let start_base = loc.start.base;
            let start_cds_from = loc.start.cds_from;
            let end_base = loc.end.base;
            let end_cds_from = loc.end.cds_from;
            let loc_start_offset = loc.start.offset.unwrap_or(0);
            let loc_end_offset = loc.end.offset.unwrap_or(0);

            // Update is_intronic flag with information from var_c.
            // From hgvs spec:
            // > Base-Offset coordinates use a base position,
            //   which is an index in the specified sequence,
            //   and an optional offset from that base position.
            //   Non-zero offsets refer to non-coding sequence,
            //   such as 5’ UTR, 3’ UTR, or intronic position.
            let is_intronic_or_utr = loc_start_offset != 0 && loc_end_offset != 0;

            // The variables below mean "VARIANT_{starts,stops}_{left,right}_OF_{start,stop}_CODON".
            //
            // start codon
            let starts_left_of_start = start_cds_from == CdsFrom::Start && start_base < 0;
            let ends_right_of_start = start_cds_from != CdsFrom::Start || start_base > 0;
            if starts_left_of_start && ends_right_of_start {
                consequences |= Consequence::StartLost;
            }
            // stop codon
            let starts_left_of_stop = start_cds_from == CdsFrom::Start;
            let ends_right_of_stop = end_cds_from == CdsFrom::End;
            if starts_left_of_stop && ends_right_of_stop && !incomplete_3p {
                consequences |= Consequence::StopLost;
            }

            if incomplete_3p
                && let Some(cds_len) = available_cds_len
                && start_base <= cds_len
                && end_base >= cds_len - 2
            {
                consequences |= Consequence::IncompleteTerminalCodonVariant;
            }

            if (start_cds_from == CdsFrom::Start && start_base <= 0)
                && ends_right_of_stop
                && matches!(edit, NaEdit::DelNum { .. } | NaEdit::DelRef { .. })
            {
                consequences |= Consequence::TranscriptAblation;
            }

            // Detect variants affecting the 5'/3' UTRs.
            if starts_left_of_start && start_base < 0 {
                if is_intronic_or_utr {
                    consequences |= Consequence::FivePrimeUtrIntronVariant;
                } else if is_exonic {
                    consequences |= Consequence::FivePrimeUtrExonVariant;
                }
            }
            if ends_right_of_stop {
                if is_intronic_or_utr {
                    consequences |= Consequence::ThreePrimeUtrIntronVariant;
                } else if is_exonic {
                    consequences |= Consequence::ThreePrimeUtrExonVariant;
                }
            }

            if matches!(edit, NaEdit::DelNum { .. } | NaEdit::DelRef { .. })
                && end_base > start_base
                && end_base > 0
            {
                if start_base <= 3 && start_cds_from == CdsFrom::Start {
                    consequences |= Consequence::StartLost;
                }

                if loc_start_offset < 0 && loc_end_offset > 0 {
                    consequences |= Consequence::ExonLossVariant;
                }
            }

            // Make sure not to report frameshift variants
            // that occur completely within intronic sequence
            // i.e. not within the CDS, as the definition is
            // "A sequence variant which causes a disruption of the translational reading frame,
            // because the number of nucleotides inserted or deleted is not a multiple of three."
            let within_exonic_sequence = loc_start_offset == 0 && loc_end_offset == 0;
            let _crosses_boundary = (loc_start_offset != 0) ^ (loc_end_offset != 0);

            if !(ends_right_of_stop || starts_left_of_start || is_intronic_or_utr) {
                match edit {
                    NaEdit::RefAlt {
                        reference,
                        alternative,
                    } => {
                        if reference.len().abs_diff(alternative.len()) % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else {
                            // Check for inframe insertions/deletions (that are not delins)
                            match (
                                reference.len().cmp(&alternative.len()),
                                alternative.is_empty() ^ reference.is_empty(),
                            ) {
                                (Ordering::Less, true) => {
                                    if conservative {
                                        consequences |= Consequence::ConservativeInframeInsertion;
                                    } else {
                                        consequences |= Consequence::DisruptiveInframeInsertion;
                                    }
                                }
                                (Ordering::Greater, true) => {
                                    if conservative {
                                        consequences |= Consequence::ConservativeInframeDeletion;
                                    } else {
                                        consequences |= Consequence::DisruptiveInframeDeletion;
                                    }
                                }
                                _ => {}
                            }
                        }
                    }
                    NaEdit::DelRef { reference } => {
                        if reference.len() % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else if conservative {
                            consequences |= Consequence::ConservativeInframeDeletion;
                        } else {
                            consequences |= Consequence::DisruptiveInframeDeletion;
                        }
                    }
                    NaEdit::DelNum { count } => {
                        if count % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else if conservative {
                            consequences |= Consequence::ConservativeInframeDeletion;
                        } else {
                            consequences |= Consequence::DisruptiveInframeDeletion;
                        }
                    }
                    NaEdit::Ins { alternative } => {
                        if alternative.len() % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else if conservative {
                            consequences |= Consequence::ConservativeInframeInsertion;
                        } else {
                            consequences |= Consequence::DisruptiveInframeInsertion;
                        }
                    }
                    NaEdit::Dup { reference } => {
                        if reference.len() % 3 != 0 {
                            if within_exonic_sequence {
                                consequences |= Consequence::FrameshiftVariant;
                            }
                        } else if conservative {
                            consequences |= Consequence::ConservativeInframeInsertion;
                        } else {
                            consequences |= Consequence::DisruptiveInframeInsertion;
                        }
                    }
                    _ => {}
                }
            }
        } else {
            panic!("Must be CDS variant: {}", &var_c)
        };
        consequences
    }

    fn analyze_protein_variant(
        &self,
        var_c: &HgvsVariant,
        var_p: &HgvsVariant,
        protein_pos: &Option<Pos>,
        conservative: bool,
        tx_accession: &str,
        incomplete_3p: bool,
    ) -> Consequences {
        let mut consequences: Consequences = Consequences::empty();

        // TODO move to hgvs-rs library as method of `ProtPos` or similar
        fn is_stop(s: &str) -> bool {
            s == "X" || s == "Ter" || s == "*"
        }

        fn has_stop(s: &str) -> bool {
            s.contains('*') || s.contains('X') || s.contains("Ter")
        }

        match var_p {
            HgvsVariant::ProtVariant { loc_edit, .. } => match loc_edit {
                ProtLocEdit::Ordinary { loc, edit } => {
                    let loc = loc.inner();
                    match edit.inner() {
                        ProteinEdit::Fs { .. } => {
                            consequences |= Consequence::FrameshiftVariant;

                            // in the case of frameshifts, we will get the altered protein sequence
                            // in order to compare it with the unaltered one

                            if let Ok(reference_data) =
                                hgvs::mapper::altseq::ref_transcript_data_cached(
                                    self.provider.clone(),
                                    tx_accession,
                                    None,
                                )
                            {
                                let original_sequence_len = reference_data.aa_sequence.len();
                                if let Ok(alt_data) =
                                    AltSeqBuilder::new(var_c.clone(), reference_data).build_altseq()
                                    && let Some(alt_data) = alt_data.first()
                                {
                                    let altered_sequence = &alt_data.aa_sequence;

                                    // trim altered sequence to the first stop encountered
                                    let altered_sequence =
                                        if let Some(pos) = altered_sequence.find('*') {
                                            // do not use the 'X' fallback here,
                                            // as that is _usually_ only added
                                            // when the number of bases is not divisible by 3.
                                            // We only want to identify cases where a new/later
                                            // stop codon is encountered
                                            // .or_else(|| altered_sequence.find('X'))
                                            &altered_sequence[..=pos]
                                        } else {
                                            altered_sequence
                                        };

                                    match altered_sequence.len().cmp(&original_sequence_len) {
                                        Ordering::Less => {
                                            consequences |= Consequence::FrameshiftTruncation;
                                        }
                                        Ordering::Equal => {
                                            if !self.config.vep_consequence_terms {
                                                consequences |= Consequence::MissenseVariant;
                                                // TODO: discuss stop_retained
                                                // consequences |= Consequence::StopRetainedVariant;
                                                consequences &= !Consequence::FrameshiftVariant;
                                            }
                                        }
                                        Ordering::Greater => {
                                            consequences |= Consequence::FrameshiftElongation;
                                        }
                                    }
                                }
                            }
                        }
                        ProteinEdit::Ext { .. } => {
                            if !incomplete_3p {
                                consequences |= Consequence::StopLost;
                            }
                            consequences |= Consequence::FeatureElongation;
                        }
                        ProteinEdit::Subst { alternative } => {
                            if alternative.is_empty() {
                                consequences |= Consequence::SynonymousVariant;
                            } else if is_stop(alternative) {
                                if loc.start == loc.end && is_stop(&loc.start.aa) {
                                    consequences |= Consequence::StopRetainedVariant;
                                } else {
                                    consequences |= Consequence::StopGained;
                                    // if the substitution happens right before the stop codon
                                    // and if it is a conservative change
                                    // then it is not a stop gained
                                    // cf. 1:43450470:GCCT:G, ENST00000634258.3:c.10294_10296del/p.Leu3432Ter
                                    if let Some(p) = protein_pos
                                        && p.total.is_some_and(|t| p.ord == t - 1)
                                        && conservative
                                    {
                                        consequences &= !Consequence::StopGained;
                                        consequences |= Consequence::ConservativeInframeDeletion;
                                    }
                                }
                            } else {
                                consequences |= Consequence::MissenseVariant;
                                // Missense variants that affect selenocysteine are marked
                                // as rare amino acid variants / selenocysteine gain/loss variants.
                                let alt_has_selenocysteine = alternative.contains('U');
                                let ref_has_selenocysteine =
                                    (loc.start == loc.end) && loc.start.aa == "U";
                                match (ref_has_selenocysteine, alt_has_selenocysteine) {
                                    (true, false) => {
                                        consequences |= Consequence::SelenocysteineLoss;
                                    }
                                    (false, true) => {
                                        consequences |= Consequence::SelenocysteineGain;
                                    }
                                    (true, true) => {
                                        consequences |= Consequence::RareAminoAcidVariant;
                                    }
                                    _ => {}
                                }
                            }
                        }
                        ProteinEdit::DelIns { alternative } => {
                            consequences |= Consequence::ProteinAlteringVariant;
                            if alternative
                                .len()
                                .cmp(&(loc.start.number.abs_diff(loc.end.number) as usize + 1))
                                == Ordering::Equal
                            {
                                // When the delins does not change the CDS length,
                                // it is a missense variant, not an inframe deletion
                                // cf https://github.com/Ensembl/ensembl-vep/issues/1388
                                consequences |= Consequence::MissenseVariant;
                                consequences &= !Consequence::ProteinAlteringVariant;
                            }

                            if (is_stop(&loc.start.aa) || is_stop(&loc.end.aa))
                                && !has_stop(alternative)
                                && !incomplete_3p
                            {
                                consequences |= Consequence::StopLost;
                            }

                            if has_stop(alternative) {
                                consequences |= Consequence::StopGained;
                            }
                        }
                        ProteinEdit::Ins { .. } | ProteinEdit::Dup => {
                            if conservative {
                                consequences |= Consequence::ConservativeInframeInsertion;
                            } else {
                                consequences |= Consequence::DisruptiveInframeInsertion;
                            }
                        }
                        ProteinEdit::Del => {
                            if conservative {
                                consequences |= Consequence::ConservativeInframeDeletion;
                            } else {
                                consequences |= Consequence::DisruptiveInframeDeletion;
                            }
                        }
                        ProteinEdit::Ident => {
                            if loc.start == loc.end && is_stop(&loc.start.aa) {
                                consequences |= Consequence::StopRetainedVariant;
                            } else {
                                consequences |= Consequence::SynonymousVariant;
                            }
                        }
                    };
                }
                ProtLocEdit::NoChange | ProtLocEdit::NoChangeUncertain => {
                    consequences |= Consequence::SynonymousVariant;
                }
                ProtLocEdit::InitiationUncertain => {
                    consequences |= Consequence::StartLost;
                }
                ProtLocEdit::NoProtein | ProtLocEdit::NoProteinUncertain | ProtLocEdit::Unknown => {
                }
            },
            _ => panic!("Must be protein variant: {}", &var_p),
        }
        consequences
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

    /// Predict the combined consequence of multiple (phased) variants.
    /// This is an _experimental_ feature.
    pub fn predict_multiple(
        &self,
        variants: &[VcfVariant],
    ) -> Result<Option<Vec<AnnField>>, SeqvarsError> {
        // Run each input variant through single-variant normalization first.
        let normalized_variants: Vec<VcfVariant> = variants
            .iter()
            .map(|v| {
                let (_, r1, a1) =
                    hgvs::sequences::trim_common_suffixes(&v.reference, &v.alternative);
                let (prefix_trim, r2, a2) = hgvs::sequences::trim_common_prefixes(&r1, &a1);
                VcfVariant {
                    chromosome: v.chromosome.clone(),
                    position: v.position + prefix_trim as i32,
                    reference: r2,
                    alternative: a2,
                }
            })
            .collect();

        let mut paired_vars: Vec<_> = variants.iter().cloned().zip(normalized_variants).collect();
        paired_vars.sort_by_key(|(_, norm)| norm.position);

        let sorted_originals: Vec<VcfVariant> =
            paired_vars.iter().map(|(orig, _)| orig.clone()).collect();
        let sorted_normalized: Vec<VcfVariant> =
            paired_vars.into_iter().map(|(_, norm)| norm).collect();

        let sorted_vars = Self::validate_and_sort_variant_group(&sorted_normalized)?;
        if sorted_vars.is_empty() {
            return Ok(None);
        }

        let chrom_acc = self
            .provider
            .contig_manager
            .get_accession(&sorted_vars[0].chromosome)
            .ok_or_else(|| SeqvarsError::UnknownChromosomeAccession)?;

        // build a pseudo hgvs.g description
        let min_pos = sorted_vars.first().unwrap().position;
        let max_var = sorted_vars.last().unwrap();
        let max_pos = max_var.position + max_var.reference.len() as i32 - 1;

        let mut hgvs_g = None;
        if let Ok(ref_seq_g) = self.provider.get_seq_part(
            chrom_acc,
            Some((min_pos as usize).saturating_sub(1)),
            Some(max_pos as usize),
        ) {
            let mut alt_seq_g = ref_seq_g.clone();

            let mut g_edits: Vec<_> = sorted_vars
                .iter()
                .map(|var| {
                    let start = (var.position - min_pos) as usize;
                    let end = start + var.reference.len();
                    (start, end, var.alternative.clone())
                })
                .collect();

            g_edits.sort_by(|a, b| b.0.cmp(&a.0));
            for (start, end, alt) in g_edits {
                alt_seq_g.replace_range(start..end, &alt);
            }

            let compound_var_g = HgvsVariant::GenomeVariant {
                accession: Accession::new(chrom_acc),
                gene_symbol: None,
                loc_edit: GenomeLocEdit {
                    loc: Mu::Certain(GenomeInterval {
                        start: Some(min_pos),
                        end: Some(max_pos),
                    }),
                    edit: Mu::Certain(NaEdit::RefAlt {
                        reference: ref_seq_g,
                        alternative: alt_seq_g,
                    }),
                },
            };
            hgvs_g = Some(
                format!("{}", &NoRef(&compound_var_g))
                    .split(':')
                    .nth(1)
                    .unwrap()
                    .to_owned(),
            );
        }

        let txs = self.get_transcripts_for_variant_group(&sorted_vars, chrom_acc)?;
        if txs.is_empty() {
            return Ok(None);
        }

        let mut multi_anns = Vec::new();

        for tx_record in txs {
            if let Some(ann) = self.predict_multiple_for_transcript(
                &sorted_vars,
                &sorted_originals,
                &tx_record,
                chrom_acc,
                hgvs_g.clone(),
            )? {
                multi_anns.push(ann);
            }
        }

        Ok(Some(self.filter_ann_fields(multi_anns)))
    }

    /// Ensures variants don't overlap and are in the correct order.
    fn validate_and_sort_variant_group(
        variants: &[VcfVariant],
    ) -> Result<Vec<VcfVariant>, SeqvarsError> {
        if variants.is_empty() {
            return Ok(Vec::new());
        }

        let chrom = &variants[0].chromosome;
        if !variants.iter().all(|v| v.chromosome == *chrom) {
            return Err(SeqvarsError::GroupValidation(
                GroupValidationError::DifferentChromosomes,
            ));
        }

        let mut sorted_vars = variants.to_vec();
        sorted_vars.sort_by_key(|v| v.position);

        for window in sorted_vars.windows(2) {
            let v1 = &window[0];
            let v2 = &window[1];
            let v1_end = v1.position + v1.reference.len() as i32 - 1;

            if v1_end >= v2.position {
                return Err(SeqvarsError::GroupValidation(
                    GroupValidationError::OverlappingVariants(v1.position, v2.position),
                ));
            }
        }

        Ok(sorted_vars)
    }

    /// Fetches all relevant transcripts for the group of variants.
    fn get_transcripts_for_variant_group(
        &self,
        sorted_vars: &[VcfVariant],
        chrom_acc: &str,
    ) -> Result<Vec<TxForRegionRecord>, SeqvarsError> {
        let min_pos = sorted_vars.first().unwrap().position;
        let max_pos = sorted_vars.last().unwrap().position
            + sorted_vars.last().unwrap().reference.len() as i32
            - 1;

        let mut txs = self
            .provider
            .get_tx_for_region(
                chrom_acc,
                ALT_ALN_METHOD,
                min_pos - PADDING,
                max_pos + PADDING,
            )
            .map_err(|e| SeqvarsError::Provider(e.to_string()))?;
        txs.sort_by(|a, b| a.tx_ac.cmp(&b.tx_ac));
        Ok(self.filter_picked_sourced_txs(txs))
    }

    /// Evaluates a group of variants for a specific transcript.
    /// Returns a tuple of projections for exonic variants (used for cDNA sequence assembly)
    /// and the combined baseline consequences of all variants in the group.
    fn get_group_projections(
        &self,
        sorted_vars: &[VcfVariant],
        tx: &Transcript,
        chrom_acc: &str,
        transcript_biotype: TranscriptBiotype,
    ) -> Result<Option<(Vec<HgvsProjectionContext>, Consequences)>, SeqvarsError> {
        let alignment = tx.genome_alignments.first().unwrap();
        let strand = Strand::try_from(alignment.strand).expect("invalid strand");

        let splice_csqs = Consequence::SpliceAcceptorVariant
            | Consequence::SpliceDonorVariant
            | Consequence::SpliceRegionVariant
            | Consequence::SplicePolypyrimidineTractVariant
            | Consequence::SpliceDonorRegionVariant
            | Consequence::SpliceDonorFifthBaseVariant
            | Consequence::ExonicSpliceRegionVariant;

        let mut projections = Vec::new();
        let mut group_consequences = Consequences::empty();

        for var in sorted_vars {
            let var_g = Self::get_var_g(var, chrom_acc);
            let (var_start, var_end) = Self::get_var_start_end(&var_g);
            let (tx_loc, tx_csqs, _) =
                self.determine_transcript_context(alignment, strand, &var_g, var_start, var_end);

            group_consequences |= tx_csqs;

            if tx_csqs.intersects(splice_csqs) {
                tracing::warn!(
                    "Phased variant {:?} affects splicing on {}. Skipping multi-variant assembly.",
                    var,
                    tx.id
                );
                return Ok(None);
            }

            if tx_loc.is_intronic && !tx_loc.is_exonic {
                tracing::debug!(
                    "Skipping purely intronic variant {:?} for sequence assembly on {}",
                    var,
                    tx.id
                );
                continue;
            }

            if !tx_loc.is_exonic {
                return Ok(None);
            }

            let proj = self.project_hgvs(&var_g, tx, transcript_biotype)?;
            if proj.n.is_none() || proj.c.is_none() || !proj.is_within_cds_bounds() {
                return Ok(None);
            }

            projections.push(proj);
        }

        Ok(Some((projections, group_consequences)))
    }

    /// Assembles the compound `delins` sequence and constructs the final annotation field.
    fn predict_multiple_for_transcript(
        &self,
        sorted_vars: &[VcfVariant],
        sorted_originals: &[VcfVariant],
        tx_record: &TxForRegionRecord,
        chrom_acc: &str,
        hgvs_g: Option<String>,
    ) -> Result<Option<AnnField>, SeqvarsError> {
        let tx = match self.provider.get_tx(&tx_record.tx_ac) {
            Some(t) => t,
            None => return Ok(None),
        };

        let transcript_biotype = TranscriptBiotype::try_from(tx.biotype).unwrap();
        if transcript_biotype != TranscriptBiotype::Coding || tx.start_codon.is_none() {
            return Ok(None);
        }

        let (projections, base_group_consequences) =
            match self.get_group_projections(sorted_vars, tx, chrom_acc, transcript_biotype)? {
                Some((p, csqs)) => (p, csqs),
                None => return Ok(None),
            };

        if projections.is_empty() {
            // All variants in this group were purely intronic and thus skipped.
            return Ok(None);
        }

        let ref_data = match hgvs::mapper::altseq::ref_transcript_data_cached(
            self.provider.clone(),
            &tx.id,
            None,
        ) {
            Ok(r) => r,
            Err(_) => return Ok(None),
        };

        struct NEdit {
            n_loc_start: i32,
            n_loc_end: i32,
            replace_start: usize,
            replace_end: usize,
            alt: String,
        }
        let mut n_edits = Vec::new();

        for proj in &projections {
            if let Some(HgvsVariant::TxVariant { loc_edit, .. }) = &proj.n {
                let loc = loc_edit.loc.inner();
                let edit = loc_edit.edit.inner();

                let n_loc_start = loc.start.base;
                let n_loc_end = loc.end.base;

                if n_loc_start > n_loc_end {
                    tracing::warn!(
                        "Invalid transcript coordinates ({} > {}) after HGVS projection. \
                        Skipping compound prediction for transcript {}.",
                        n_loc_start,
                        n_loc_end,
                        tx.id
                    );
                    return Ok(None);
                }

                let replace_start;
                let replace_end;
                let alt;

                match edit {
                    NaEdit::RefAlt { alternative, .. } | NaEdit::NumAlt { alternative, .. } => {
                        replace_start = (n_loc_start - 1) as usize;
                        replace_end = n_loc_end as usize;
                        alt = alternative.clone();
                    }
                    NaEdit::DelRef { .. } | NaEdit::DelNum { .. } => {
                        replace_start = (n_loc_start - 1) as usize;
                        replace_end = n_loc_end as usize;
                        alt = "".to_string();
                    }
                    NaEdit::Dup { reference } => {
                        replace_start = (n_loc_start - 1) as usize;
                        replace_end = n_loc_end as usize;
                        alt = format!("{}{}", reference, reference);
                    }
                    NaEdit::Ins { alternative } => {
                        replace_start = n_loc_start as usize;
                        replace_end = n_loc_start as usize;
                        alt = alternative.clone();
                    }
                    _ => {
                        tracing::warn!(
                            "Unsupported NaEdit type {:?} in multi-assembly. Skipping.",
                            edit
                        );
                        return Ok(None);
                    }
                }

                n_edits.push(NEdit {
                    n_loc_start,
                    n_loc_end,
                    replace_start,
                    replace_end,
                    alt,
                });
            }
        }

        n_edits.sort_by(|a, b| b.replace_start.cmp(&a.replace_start));

        let tx_len = ref_data.transcript_sequence.len() as i32;
        let mut alt_seq = ref_data.transcript_sequence.clone();
        let n_min = n_edits.iter().map(|e| e.n_loc_start).min().unwrap();
        let n_max = n_edits.iter().map(|e| e.n_loc_end).max().unwrap();

        let mut total_delta = 0i32;
        for edit in &n_edits {
            if edit.replace_start > alt_seq.len() || edit.replace_end > alt_seq.len() {
                tracing::warn!(
                    "Edit range out of bounds: {}..{} exceeds sequence length {}. Cannot assemble variant.",
                    edit.replace_start,
                    edit.replace_end,
                    alt_seq.len()
                );
                return Ok(None);
            }

            if edit.replace_start > edit.replace_end {
                tracing::warn!(
                    "Invalid edit range: start {} > end {}. Cannot assemble variant.",
                    edit.replace_start,
                    edit.replace_end
                );
                return Ok(None);
            }

            alt_seq.replace_range(edit.replace_start..edit.replace_end, &edit.alt);
            let orig_len = edit.replace_end - edit.replace_start;
            total_delta += edit.alt.len() as i32 - orig_len as i32;
        }

        let new_length = (n_max - n_min + 1) + total_delta;
        if new_length < 0 {
            tracing::warn!("Calculated new_length {} is negative.", new_length);
            return Ok(None);
        }
        let start_idx = (n_min - 1) as usize;
        let end_idx = start_idx + new_length as usize;

        if start_idx > alt_seq.len() || end_idx > alt_seq.len() {
            tracing::warn!(
                "Slice index out of bounds: start={} end={} > {} (alt_seq length). Cannot assemble variant.",
                start_idx,
                end_idx,
                alt_seq.len()
            );
            return Ok(None);
        }

        let new_substring = &alt_seq[start_idx..end_idx];
        let ref_substring = &ref_data.transcript_sequence[(n_min - 1) as usize..n_max as usize];

        let compound_var_n = HgvsVariant::TxVariant {
            accession: projections[0].n.as_ref().unwrap().accession().clone(),
            gene_symbol: projections[0].n.as_ref().unwrap().gene_symbol().clone(),
            loc_edit: hgvs::parser::TxLocEdit {
                loc: Mu::Certain(hgvs::parser::TxInterval {
                    start: hgvs::parser::TxPos {
                        base: n_min,
                        offset: None,
                    },
                    end: hgvs::parser::TxPos {
                        base: n_max,
                        offset: None,
                    },
                }),
                edit: Mu::Certain(NaEdit::RefAlt {
                    reference: ref_substring.to_string(),
                    alternative: new_substring.to_string(),
                }),
            },
        };

        let compound_var_c = match self.mapper.n_to_c(&compound_var_n) {
            Ok(c) => c,
            Err(e) => {
                tracing::debug!("Failed to map compound n_loc to c_loc: {}", e);
                return Ok(None);
            }
        };

        let compound_var_n = HgvsVariant::TxVariant {
            accession: projections[0].n.as_ref().unwrap().accession().clone(),
            gene_symbol: projections[0].n.as_ref().unwrap().gene_symbol().clone(),
            loc_edit: hgvs::parser::TxLocEdit {
                loc: Mu::Certain(hgvs::parser::TxInterval {
                    start: hgvs::parser::TxPos {
                        base: n_min,
                        offset: None,
                    },
                    end: hgvs::parser::TxPos {
                        base: n_max,
                        offset: None,
                    },
                }),
                edit: Mu::Certain(NaEdit::RefAlt {
                    reference: ref_substring.to_string(),
                    alternative: new_substring.to_string(),
                }),
            },
        };

        let compound_var_p = self.safe_project_c_to_p(&compound_var_c)?;

        let compound_proj = HgvsProjectionContext {
            g: projections[0].g.clone(),
            n: Some(compound_var_n),
            c: Some(compound_var_c.clone()),
            p: compound_var_p,
        };

        let mut custom_fields = BTreeMap::new();
        let c_ref = self.config.report_cdna_sequence.includes_ref();
        let c_alt = self.config.report_cdna_sequence.includes_alt();
        let p_ref = self.config.report_protein_sequence.includes_ref();
        let p_alt = self.config.report_protein_sequence.includes_alt();

        if c_ref || c_alt || p_ref || p_alt {
            if c_ref {
                custom_fields.insert(
                    ANN_TX_SEQ_REF.into(),
                    Some(ref_data.transcript_sequence.clone()),
                );
            }
            if p_ref {
                custom_fields.insert(ANN_AA_SEQ_REF.into(), Some(ref_data.aa_sequence.clone()));
            }
            if (c_alt || p_alt)
                && let Ok(alt_data_vec) =
                    AltSeqBuilder::new(compound_var_c.clone(), ref_data).build_altseq()
                && let Some(alt_data) = alt_data_vec.into_iter().next()
            {
                if c_alt {
                    custom_fields.insert(ANN_TX_SEQ_ALT.into(), Some(alt_data.transcript_sequence));
                }
                if p_alt {
                    custom_fields.insert(ANN_AA_SEQ_ALT.into(), Some(alt_data.aa_sequence));
                }
            }
        }

        let tlc = TranscriptLocationContext {
            rank: Rank { ord: 1, total: 1 }, // dummy rank since we potentially span multiple exons
            distance: Some(0),
            is_exonic: true,
            is_intronic: false,
            is_upstream: false,
            is_downstream: false,
        };

        let c_ctx = self.analyze_transcript_consequences(
            &compound_proj,
            tx,
            tx_record,
            &tlc,
            tx_len,
            transcript_biotype,
        )?;

        let mut consequences =
            c_ctx.cds_consequences | c_ctx.protein_consequences | base_group_consequences;
        self.consequences_fix_special_cases(
            &mut consequences,
            c_ctx.cds_consequences,
            c_ctx.protein_consequences,
            &compound_proj,
        );

        if self.config.vep_consequence_terms {
            self.adjust_vep_terms(&mut consequences, Some(&compound_proj));
        }
        if consequences.is_empty() {
            consequences |= Consequence::GeneVariant;
        }

        let consequences_vec = consequences.iter().collect_vec();
        let putative_impact = (*consequences_vec.first().unwrap()).into();

        let hgvs_n = compound_proj.n.as_ref().map(|n| {
            format!("{}", &NoRef(n))
                .split(':')
                .nth(1)
                .unwrap()
                .to_owned()
        });
        let hgvs_c = Some(
            format!("{}", &NoRef(compound_proj.c.as_ref().unwrap()))
                .split(':')
                .nth(1)
                .unwrap()
                .to_owned(),
        );
        let hgvs_p = compound_proj
            .p
            .as_ref()
            .map(|p| format!("{}", p).split(':').nth(1).unwrap().to_owned());

        let ref_alts = sorted_originals
            .iter()
            .map(|v| v.reference.clone())
            .collect();
        let alt_alts = sorted_originals
            .iter()
            .map(|v| v.alternative.clone())
            .collect();

        let strand = match tx.genome_alignments.first().unwrap().strand {
            1 => 1,
            -1 => -1,
            _ => 0,
        };

        let feature_tags = tx
            .tags
            .iter()
            .map(|tag| TranscriptTag::try_from(*tag).expect("invalid transcript tag"))
            .filter(|tag| !matches!(tag, TranscriptTag::EnsemblGraft))
            .filter_map(|t| {
                if !matches!(FeatureTag::from(t), FeatureTag::Other(_)) {
                    Some(FeatureTag::from(t))
                } else {
                    None
                }
            })
            .collect_vec();

        Ok(Some(AnnField {
            allele: Allele::Grouped(GroupedAlleles {
                references: ref_alts,
                alternatives: alt_alts,
            }),
            consequences: consequences_vec,
            putative_impact,
            gene_symbol: tx.gene_symbol.clone(),
            gene_id: tx.gene_id.clone(),
            feature_type: FeatureType::SoTerm {
                term: SoFeature::Transcript,
            },
            feature_id: tx.id.clone(),
            feature_biotype: vec![FeatureBiotype::Coding],
            feature_tags,
            rank: None,
            hgvs_g,
            hgvs_n,
            hgvs_c,
            hgvs_p,
            cdna_pos: c_ctx.cdna_pos,
            cds_pos: c_ctx.cds_pos,
            protein_pos: c_ctx.protein_pos,
            strand,
            distance: Some(0),
            messages: None,
            custom_fields,
        }))
    }

    /// Safely projects a CDS variant to a Protein variant, gracefully catching
    /// and swallowing expected incomplete-transcript errors as `None`.
    fn safe_project_c_to_p(
        &self,
        var_c: &HgvsVariant,
    ) -> Result<Option<HgvsVariant>, SeqvarsError> {
        self.mapper.c_to_p(var_c).map_or_else(
            |e| {
                if matches!(
                    e,
                    Error::TranscriptLengthInvalid(_, _)
                        | Error::CannotConvertIntervalEnd(_)
                        | Error::MultipleAAVariants
                ) {
                    tracing::debug!("c_to_p failed gracefully (typed error): {}", e);
                    return Ok(None);
                }

                let err_str = e.to_string();
                if err_str.contains("does not contain a stop codon")
                    || err_str.contains("multiple of 3")
                    || err_str.contains("multiple of three")
                    || err_str.contains("out of bound")
                    || err_str.contains("outside of sequence bounds")
                {
                    tracing::debug!(
                        "c_to_p failed gracefully (nested error string): {}",
                        err_str
                    );
                    Ok(None)
                } else {
                    Err(SeqvarsError::HgvsProjection(format!(
                        "c_to_p mapping failed: {}",
                        e
                    )))
                }
            },
            |v| Ok(Some(v)),
        )
    }
}

fn is_conservative_cds_variant(var_c: &HgvsVariant) -> bool {
    match var_c {
        HgvsVariant::CdsVariant { loc_edit, .. } => {
            // Handle the cases where the variant touches the start or stop codon based on `var_c`
            // coordinates. The cases where the start/stop codon is touched by the variant
            // directly is handled elsewhere based on the `var_p` prediction.
            let loc = loc_edit.loc.inner();
            let start_base = loc.start.base;
            let start_cds_from = loc.start.cds_from;
            let end_base = loc.end.base;
            let end_cds_from = loc.end.cds_from;
            // The range is "conservative" (regarding deletions and insertions) if
            // it does not start or end within codons.
            start_cds_from == CdsFrom::Start
                && end_cds_from == CdsFrom::Start
                && start_base % 3 == 1
                && (end_base + 1) % 3 == 1
        }
        _ => panic!("Expected CdsVariant, got {:#?}", var_c),
    }
}

#[inline]
fn overlaps(start_a: i32, end_a: i32, start_b: i32, end_b: i32) -> bool {
    (start_a < end_b) && (end_a > start_b)
}

impl ConsequencePredictor {
    /// Return data version string (if set).
    pub fn data_version(&self) -> Option<String> {
        self.provider.as_ref().tx_seq_db.version.clone()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::annotate::cli::{PredictorSettings, TranscriptPickType, TranscriptSettings};
    use crate::annotate::seqvars::provider::ConfigBuilder as MehariProviderConfigBuilder;
    use crate::annotate::seqvars::{
        Args, AsyncAnnotatedVariantWriter, PathOutput, load_tx_db, run_with_writer,
    };
    use crate::common::TsvContigStyle;
    use crate::common::noodles::{NoodlesVariantReader, open_variant_reader, open_variant_writer};
    use csv::ReaderBuilder;
    use futures::TryStreamExt;
    use insta::assert_yaml_snapshot;
    use noodles::vcf::variant::Record as NoodlesRecord;
    use noodles::vcf::variant::record_buf::info::field::Value;
    use noodles::vcf::variant::record_buf::info::field::value::Array;
    use pretty_assertions::assert_eq;
    use serde::Deserialize;
    use std::collections::BTreeMap;
    use std::path::{Path, PathBuf};
    use std::{fs::File, io::BufReader};
    use tempfile::NamedTempFile;

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
            None::<PathBuf>,
            true,
            Default::default(),
        ));

        let predictor = ConsequencePredictor::new(provider, Default::default());

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
    #[case("17:41197820:G:T", 1, vec![Consequence::SpliceAcceptorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 1bp intronic
    #[case("17:41197821:A:C", 2, vec![Consequence::SpliceAcceptorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 2bp intronic
    #[case("17:41197822:C:A", 3, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 3bp intronic
    #[case("17:41197823:C:A", 4, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 4bp intronic
    #[case("17:41197824:T:G", 5, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 5bp intronic
    #[case("17:41197825:C:A", 6, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 6bp intronic
    #[case("17:41197835:T:G", 16, vec![Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 16bp intronic
    #[case("17:41197836:G:A", 17, vec![Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 17bp intronic
    #[case("17:41197837:G:A", 18, vec![Consequence::CodingTranscriptIntronVariant]
    )] // 18bp intronic
    #[case("17:41199660:G:T", 0, vec![Consequence::MissenseVariant, Consequence::ExonicSpliceRegionVariant]
    )] // exonic
    #[case("17:41199659:G:T", -1, vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -1bp intronic
    #[case("17:41199658:T:G", -2, vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -2bp intronic
    #[case("17:41199657:G:T", -3, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -3bp intronic
    #[case("17:41199656:A:C", -4, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -4bp intronic
    #[case("17:41199655:G:T", -5, vec![Consequence::SpliceDonorFifthBaseVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -5bp intronic
    #[case("17:41199654:G:T", -6, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -6bp intronic
    #[case("17:41199653:T:G", -7, vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -7bp intronic
    #[case("17:41199652:G:T", -8, vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -8bp intronic
    #[case("17:41199651:C:A", -9, vec![Consequence::CodingTranscriptIntronVariant]
    )] // -9bp intronic
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
            None::<PathBuf>,
            true,
            MehariProviderConfigBuilder::default()
                .pick_transcript(vec![
                    TranscriptPickType::ManePlusClinical,
                    TranscriptPickType::ManeSelect,
                    TranscriptPickType::Length,
                ])
                .build()?,
        ));

        use crate::annotate::seqvars::ConsequencePredictorConfigBuilder;
        let predictor = ConsequencePredictor::new(
            provider,
            ConsequencePredictorConfigBuilder::default()
                .report_most_severe_consequence_by(Some(ConsequenceBy::Gene))
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
    #[case("3:193332512:T:G", 0, vec![Consequence::MissenseVariant, Consequence::ExonicSpliceRegionVariant]
    )] // exonic
    #[case("3:193332511:G:T", -1, vec![Consequence::SpliceAcceptorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -1bp intronic
    #[case("3:193332510:A:G", -2, vec![Consequence::SpliceAcceptorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -2bp intronic
    #[case("3:193332509:C:T", -3, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant,  Consequence::CodingTranscriptIntronVariant]
    )] // -3bp intronic
    #[case("3:193332508:T:C", -4, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant,  Consequence::CodingTranscriptIntronVariant]
    )] // -4bp intronic
    #[case("3:193332507:T:C", -5, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant,  Consequence::CodingTranscriptIntronVariant]
    )] // -5bp intronic
    #[case("3:193332506:T:C", -6, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -6bp intronic
    #[case("3:193332505:C:G", -7, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -7bp intronic
    #[case("3:193332504:T:C", -8, vec![Consequence::SpliceRegionVariant, Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -8bp intronic
    #[case("3:193332503:T:A", -9, vec![Consequence::SplicePolypyrimidineTractVariant, Consequence::CodingTranscriptIntronVariant]
    )] // -9bp intronic
    #[case("3:193332831:G:T", 1, vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 1bp intronic
    #[case("3:193332832:T:C", 2, vec![Consequence::SpliceDonorVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 2bp intronic
    #[case("3:193332833:G:A", 3, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 3bp intronic
    #[case("3:193332834:A:C", 4, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 4bp intronic
    #[case("3:193332835:A:T", 5, vec![Consequence::SpliceDonorFifthBaseVariant, Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 5bp intronic
    #[case("3:193332836:C:A", 6, vec![Consequence::SpliceRegionVariant, Consequence::SpliceDonorRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 6bp intronic
    #[case("3:193332837:T:G", 7, vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 7bp intronic
    #[case("3:193332838:T:G", 8, vec![Consequence::SpliceRegionVariant, Consequence::CodingTranscriptIntronVariant]
    )] // 8bp intronic
    #[case("3:193332839:G:A", 9, vec![Consequence::CodingTranscriptIntronVariant])] // 9bp intronic
    #[case("3:193332846:A:G", 16, vec![Consequence::CodingTranscriptIntronVariant])] // 16bp intronic
    #[case("3:193332847:G:A", 17, vec![Consequence::CodingTranscriptIntronVariant])] // 17bp intronic
    #[case("3:193332848:T:A", 18, vec![Consequence::CodingTranscriptIntronVariant])] // 18bp intronic
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
            None::<PathBuf>,
            true,
            MehariProviderConfigBuilder::default()
                .pick_transcript(vec![
                    TranscriptPickType::ManePlusClinical,
                    TranscriptPickType::ManeSelect,
                    TranscriptPickType::Length,
                ])
                .build()?,
        ));

        let predictor = ConsequencePredictor::new(provider, Default::default());

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

    #[rstest::rstest]
    #[case("3:193311167:ATGT:T", vec![Consequence::StartLost, Consequence::ConservativeInframeDeletion]
    )]
    #[case("3:193311170:TGGC:C", vec![Consequence::ConservativeInframeDeletion])]
    #[case("3:193311170:TGGCG:G", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation]
    )]
    #[case("3:193311180:GTCG:G", vec![Consequence::DisruptiveInframeDeletion])]
    #[case("3:193409910:GAAA:G", vec![Consequence::ConservativeInframeDeletion])]
    #[case("3:193409913:ATAA:A", vec![Consequence::StopLost, Consequence::FeatureElongation, Consequence::ConservativeInframeDeletion]
    )]
    fn annotate_del_opa1_csqs(
        #[case] spdi: &str,
        #[case] expected_csqs: Vec<Consequence>,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", spdi.replace(':', "-"));

        let spdi = spdi.split(':').map(|s| s.to_string()).collect::<Vec<_>>();

        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            None::<PathBuf>,
            true,
            MehariProviderConfigBuilder::default()
                .pick_transcript(vec![
                    TranscriptPickType::ManePlusClinical,
                    TranscriptPickType::ManeSelect,
                    TranscriptPickType::Length,
                ])
                .build()?,
        ));

        let predictor = ConsequencePredictor::new(provider, Default::default());

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
        #[case] report_most_severe_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!(
            "{}-{}-{}",
            spdi.replace(':', "-"),
            pick_transcripts,
            !report_most_severe_consequence_only
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
            None::<PathBuf>,
            true,
            MehariProviderConfigBuilder::default()
                .pick_transcript(picks)
                .build()
                .unwrap(),
        ));
        let report_most_severe_consequence_by = if report_most_severe_consequence_only {
            Some(ConsequenceBy::Gene)
        } else {
            None
        };

        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .report_most_severe_consequence_by(report_most_severe_consequence_by)
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
        #[case] report_most_severe_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!(
            "{}-{}-{}",
            spdi.replace(':', "-"),
            pick_transcripts,
            !report_most_severe_consequence_only
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
            None::<PathBuf>,
            true,
            MehariProviderConfigBuilder::default()
                .pick_transcript(picks)
                .build()
                .unwrap(),
        ));

        let report_most_severe_consequence_by = if report_most_severe_consequence_only {
            Some(ConsequenceBy::Gene)
        } else {
            None
        };

        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .report_most_severe_consequence_by(report_most_severe_consequence_by)
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

    /// This is a set of variants where VEP and mehari to disagree,
    /// i.e. interesting/edge cases that are not as clear-cut as others.
    ///
    /// This test ensures we do not regress on these cases.
    #[tokio::test]
    async fn annotate_vep_disagreement_cases() -> Result<(), anyhow::Error> {
        let tx_path =
            "tests/data/annotate/db/grch38/GRCh38-ensembl.disagreement-subset.txs.bin.zst";

        let path_input_vcf = "tests/data/annotate/seqvars/vep.disagreement-cases.vcf";
        let output = NamedTempFile::new()?;
        let mut writer = open_variant_writer(output.as_ref()).await?;
        run_with_writer(
            &mut writer,
            &Args {
                threads: 1,
                reference: None,
                in_memory_reference: true,
                genome_release: None,
                path_input_ped: None,
                path_input_vcf: path_input_vcf.into(),
                output: PathOutput {
                    path_output_vcf: Some(output.as_ref().to_str().unwrap().into()),
                    path_output_tsv: None,
                },
                predictor_settings: PredictorSettings {
                    transcript_settings: TranscriptSettings {
                        report_most_severe_consequence_by: Some(ConsequenceBy::Allele),
                        pick_transcript: vec![TranscriptPickType::ManeSelect],
                        ..Default::default()
                    },
                    ..Default::default()
                },
                max_var_count: None,
                hgnc: None,
                sources: crate::annotate::seqvars::Sources {
                    transcripts: Some(vec![tx_path.into()]),
                    frequencies: None,
                    clinvar: None,
                },
                tsv_contig_style: TsvContigStyle::Auto,
            },
        )
        .await?;
        writer.shutdown().await?;

        let records_written = read_vcf(output).await?;

        let mut snapshot_data = BTreeMap::new();
        let header = noodles::vcf::io::reader::Builder::default()
            .build_from_path(path_input_vcf)?
            .read_header()?;

        for record in records_written {
            let key = format!(
                "{}:{}:{}:{}:{}",
                record.reference_sequence_name(),
                record
                    .variant_start()
                    .map_or_else(|| "0".into(), |s| s.to_string()),
                record
                    .variant_end(&header)
                    .map_or_else(|_| "0".into(), |s| s.to_string()),
                record.reference_bases(),
                record.alternate_bases().as_ref().join(",")
            );

            let ann_field = record.info().get("ANN").flatten().map(|v| match v {
                Value::Array(Array::String(inner)) => inner
                    .iter()
                    .map(|s| s.clone().unwrap_or_default())
                    .join("|"),
                _ => "".into(),
            });

            snapshot_data.insert(key, ann_field);
        }

        assert_yaml_snapshot!("vep_disagreement_cases_output", snapshot_data);

        Ok(())
    }

    async fn read_vcf(
        path: impl AsRef<Path>,
    ) -> Result<Vec<noodles::vcf::variant::RecordBuf>, anyhow::Error> {
        let mut output_reader = open_variant_reader(path.as_ref()).await?;
        let header = output_reader.read_header().await?;
        let mut record_iter = output_reader.records(&header).await;
        let mut records = Vec::new();
        while let Some(record) = record_iter.try_next().await? {
            records.push(record);
        }
        Ok(records)
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
        report_most_severe_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        let txs = vec![
            String::from("NM_001354663.2"),
            String::from("NM_001354664.2"),
            String::from("NM_015560.3"),
            String::from("NM_130831.3"),
            String::from("NM_130832.3"),
            String::from("NM_130837.3"),
        ];

        annotate_vars(path_tsv, &txs, report_most_severe_consequence_only, false)
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
        report_most_severe_consequence_only: bool,
    ) -> Result<(), anyhow::Error> {
        let txs = vec![
            String::from("NM_007294.4"),
            String::from("NM_007297.4"),
            String::from("NM_007298.3"),
            String::from("NM_007299.4"),
            String::from("NM_007300.4"),
        ];

        annotate_vars(path_tsv, &txs, report_most_severe_consequence_only, false)
    }

    fn annotate_vars(
        path_tsv: &str,
        txs: &[String],
        report_most_severe_consequence_only: bool,
        with_reference: bool,
    ) -> Result<(), anyhow::Error> {
        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let reference_path = "resources/GCF_000001405.25_GRCh37.p13_genomic.fna";
        let provider = Arc::new(MehariProvider::new(
            tx_db,
            with_reference.then_some(reference_path),
            true,
            Default::default(),
        ));

        let report_most_severe_consequence_by = if report_most_severe_consequence_only {
            Some(ConsequenceBy::Gene)
        } else {
            None
        };

        let predictor = ConsequencePredictor::new(
            provider,
            ConfigBuilder::default()
                .report_most_severe_consequence_by(report_most_severe_consequence_by)
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

            // Because for this variant our highest impact is "HIGH" and vep's is not (and because we filter for the highest impact), skip it.
            // We predict FrameshiftVariant, FrameshiftTruncation, SpliceRegionVariant and ThreePrimeUtrExonVariant,
            // while vep calls "splice_region_variant", "coding_sequence_variant", "3_prime_UTR_variant".
            if record.var
                == "17-41196310-GGTGGAAGTGTTTGCTACCAAGTTTATTTGCAGTGTTAACAGCACAACATTTACAAAACGTATTTTGTACAATCAAGTCTTCACTGCCCTTGCACACTGGGGGGGCTAGGGAAGACCTAGTCCTTCCAACAGCTATAAACAGTCCTGGATAATGGGTTTATGAAAAACACTTTTTCTTCCTTCAGCAAGCAAAATTATTTATGAAGCTGTATGGTTTCAGCAACAGGGAGCAAAGGAAAAAAATCACCTCAAAGAAAGCAACAGCTTCCTTCCTGGTGGGATCTGTCATTTTATAGATATGAAATATTCATGCCAGAGGTCTTATATTTTAAGAGGAATGGATTATATACCAGAGCTACAACAATAAACATTTTACTTATTACTAATGAGGAATTAGAAGACTGTCTTTGGAAACCGGTTCTTGAAAATCTTCTGCTGTTTTAGAACACATTCTTTAGAAATCTAGCAAATATATCTCAGACTTTTAGAAATCTCTTCTAGTTTCATTTTCCTTTTTTTTTTTTTTTTTTTGAGCCACAGTCTCACTGTCACCCAGGCTGGAGTGCCGTGGTATGATCTTGGCTCACTGCAACCTCCACCTCCCGGGCTGAAGTGATTCTCCTGCCTTAGCCACCTGAGTAGCTGGGATTACAGGTGTCCACCACCATGACCGGCTAATTTCTGTATTTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTTTCGAACTCCTGACCTCCAGTGATCTGCCCACCTTGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCATGCCCAGGTTTCAAGTTTCCTTTTCATTTCTAATACCTGCCTCAGAATTTCCTCCCCAATGTTCCACTCCAACATTTGAGAACTGCCCAAGGACTATTCTGACTTTAAGTCACATAATCGATCCCAAGCACTCTCCTTCCATTGAAGGGTCTGACTCTCTGCCTTTGTGAACACAGGGTTTTAGAGAAGTAAACTTAGGGAAACCAGCTATTCTCTTGAGGCCAAGCCACTCTGTGCTTCCAGCCCTAAGCCAACAACAGCCTGAATAGAAAGAATAGGGCTGATAAATAATGAATCAGCATCTTGCTCAATTGGTGGCGTTTAAATGGTTTTAAAATCTTCTCAGGTGAAAAATTACCATAATTTTGTGCTCATGGCAGATTTCCAAGGGAGACTTCAAGCAGAAAATCTTTAAGGGACCCTTGCATAGCCAGAAGTCCTTTTCAGGCTGATGTACATAAAATATTTAGTAGCCAGGACAGTAGAAGGACTGAAGAGTGAGAGGAGCTCCCAGGGCCTGGAAAGGCCACTTTGTAAGCTCATTCTTGGGGTCCTGTGGCTCTGTACCTGTGGCTGGCTGCAGTCAGTAGTGGCTGTGGGGGATCTGGGGTATCAGGTAGGTGTCCAGCTCCTGGCACTGGTAGAGTGCTACACTGTCCAACACCCACTCTCGGGTCACCACAGGTGCCTCACACATCTGCCCAATT-G"
            {
                continue;
            }

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
                    let record_csqs = record.csq.split('&').collect::<Vec<_>>();

                    let highest_impact = ann.consequences.first().unwrap().impact();
                    let expected_one_of = ann
                        .consequences
                        .iter()
                        .filter(|csq| csq.impact() == highest_impact)
                        .map(|csq| csq.to_string())
                        .collect::<Vec<_>>();
                    let mut expected_one_of = expected_one_of
                        .iter()
                        .map(|s| s.as_str())
                        .collect::<Vec<_>>();

                    // Map effects a bit for VEP.
                    if path_tsv.contains(".vep")
                        && (expected_one_of.contains(&"disruptive_inframe_deletion")
                            || expected_one_of.contains(&"conservative_inframe_deletion"))
                    {
                        expected_one_of.push("inframe_deletion");
                    }

                    let found_one = [
                        // Try to find a direct match.
                        record_csqs.iter().any(|csq| expected_one_of.contains(csq)),
                        // vep sometimes only reports a coding_sequence_variant, so we accept anything
                        path_tsv.contains(".vep")
                            && (record_csqs == ["coding_sequence_variant"]
                                && !expected_one_of.is_empty()),
                        // It is common that the other tool predicts a frameshift variant while the actual prediction
                        // is stop_gained or stop_lost.  We thus also check for this case and allow it.
                        (record_csqs.contains(&"frameshift_variant")
                            || record_csqs.contains(&"frameshift_truncation")
                            || record_csqs.contains(&"frameshift_elongation"))
                            && (expected_one_of.contains(&"stop_gained"))
                            || expected_one_of.contains(&"stop_lost"),
                        // … or vice-versa
                        (expected_one_of.contains(&"frameshift_variant")
                            || expected_one_of.contains(&"frameshift_truncation")
                            || expected_one_of.contains(&"frameshift_elongation"))
                            && (record_csqs.contains(&"stop_gained"))
                            || record_csqs.contains(&"stop_lost"),
                        // VEP does not differentiate between disruptive and conservative inframe deletions and insertions.
                        (record_csqs.contains(&"inframe_deletion")
                            && (expected_one_of.contains(&"disruptive_inframe_deletion")
                                || expected_one_of.contains(&"conservative_inframe_deletion")))
                            || (record_csqs.contains(&"inframe_insertion")
                                && (expected_one_of.contains(&"disruptive_inframe_insertion")
                                    || expected_one_of
                                        .contains(&"conservative_inframe_insertion"))),
                        // delins on protein level are sometimes erroneously reported as inframe_insertion/deletion instead
                        (expected_one_of.contains(&"protein_altering_variant")
                            && ([
                                "disruptive_inframe_deletion",
                                "disruptive_inframe_insertion",
                                "conservative_inframe_deletion",
                                "conservative_inframe_insertion",
                            ]
                            .iter()
                            .any(|c| record_csqs.contains(c)))),
                        // NB: We cannot predict 5_prime_UTR_premature_start_codon_gain_variant yet. For now, we
                        // also accept 5_prime_UTR_variant.
                        ((expected_one_of.contains(&"5_prime_UTR_exon_variant")
                            || expected_one_of.contains(&"5_prime_UTR_intron_variant"))
                            && (record_csqs
                                .contains(&"5_prime_UTR_premature_start_codon_gain_variant"))),
                        // We accept 5_prime_UTR_exon_variant and 5_prime_UTR_intron_variant if the
                        // other tool predicts upstream_gene_variant.
                        ((expected_one_of.contains(&"5_prime_UTR_exon_variant")
                            || expected_one_of.contains(&"5_prime_UTR_intron_variant"))
                            && (record_csqs.contains(&"upstream_gene_variant"))),
                        // A coding_transcript_intron_variant is a more specific intron_variant
                        (expected_one_of.contains(&"coding_transcript_intron_variant")
                            && (record_csqs.contains(&"intron_variant"))),
                        // VEP predicts `splice_donor_5th_base_variant` rather than `splice_region_variant`.
                        // Same for `splice_donor_region_variant`.
                        (expected_one_of.contains(&"splice_region_variant")
                            && (record_csqs.contains(&"splice_donor_5th_base_variant")
                                || record_csqs.contains(&"splice_donor_region_variant"))),
                        // In the case of insertions at the end of an exon, VEP predicts `splice_region_variant`
                        // while we predict `splice_donor_variant`, same for start.
                        (expected_one_of.contains(&"splice_donor_variant")
                            || expected_one_of.contains(&"splice_acceptor_variant"))
                            && (record_csqs.contains(&"splice_region_variant")),
                        // VEP sometimes mispredicts disruptive inframe deletion as missense...
                        // cf. https://github.com/Ensembl/ensembl-vep/issues/1388
                        expected_one_of.contains(&"disruptive_inframe_deletion")
                            && (record_csqs.contains(&"missense_variant")),
                        // VEP does not provide `exon_loss_variant`, so we also accept `inframe_deletion` and
                        // `splice_region_variant` (BRA1 test case).
                        expected_one_of.contains(&"exon_loss_variant")
                            && (record_csqs.contains(&"inframe_deletion")
                                || record_csqs.contains(&"splice_region_variant")),
                        // On BRCA1, there is a case where VEP predicts `protein_altering_variant` rather than
                        // `disruptive_inframe_deletion`.  We accept this as well.
                        (expected_one_of.contains(&"disruptive_inframe_deletion")
                            || expected_one_of.contains(&"inframe_indel"))
                            && (record_csqs.contains(&"protein_altering_variant")),
                        // In the case of `GRCh37:17:41258543:T:TA`, the `hgvs` prediction is `c.-1_1insT` and
                        // `p.Met1?` which leads to `start_lost` while VEP predicts `5_prime_UTR_variant`.
                        // This may be a bug in `hgvs` and we don't change this for now.  We accept the call
                        // by VEP, of course.
                        expected_one_of.contains(&"start_lost")
                            && (record_csqs.contains(&"5_prime_UTR_variant")),
                        // We have specialized {5,3}_prime_UTR_{exon,intron}_variant handling, while
                        // vep and snpEff do not
                        record_csqs.contains(&"5_prime_UTR_variant")
                            && (expected_one_of.contains(&"5_prime_UTR_exon_variant")
                                || expected_one_of.contains(&"5_prime_UTR_intron_variant")),
                        record_csqs.contains(&"3_prime_UTR_variant")
                            && (expected_one_of.contains(&"3_prime_UTR_exon_variant")
                                || expected_one_of.contains(&"3_prime_UTR_intron_variant")),
                        // an inframe_indel can be a missense_variant if it is an MNV (which we do not explicitly check here)
                        expected_one_of.contains(&"inframe_indel")
                            && (record_csqs.contains(&"missense_variant")),
                        // inframe_indel also is a superclass of *_inframe_{deletion, insertion}
                        expected_one_of.contains(&"inframe_indel")
                            && [
                                "disruptive_inframe_deletion",
                                "conservative_inframe_deletion",
                                "inframe_deletion",
                                "disruptive_inframe_insertion",
                                "conservative_inframe_insertion",
                                "inframe_insertion",
                            ]
                            .iter()
                            .any(|c| record_csqs.contains(c)),
                        // SnpEff has a different interpretation of disruptive/conservative inframe deletions.
                        // We thus allow both.
                        expected_one_of.contains(&"disruptive_inframe_deletion")
                            && (record_csqs.contains(&"conservative_inframe_deletion"))
                            || expected_one_of.contains(&"disruptive_inframe_insertion")
                                && (record_csqs.contains(&"conservative_inframe_insertion")),
                        // SnpEff may not predict `splice_region_variant` for 5' UTR correctly, so we
                        // allow this.
                        (expected_one_of.contains(&"splice_region_variant")
                            || expected_one_of.contains(&"exonic_splice_region_variant"))
                            && (record_csqs.contains(&"5_prime_UTR_variant")),
                        // SnpEff does not predict `splice_polypyrimidine_tract_variant`
                        expected_one_of.contains(&"splice_polypyrimidine_tract_variant")
                            && (record_csqs.contains(&"splice_region_variant")
                                || record_csqs.contains(&"intron_variant")),
                        // For `GRCh37:3:193366573:A:ATATTGCCTAGAATGAACT`, SnpEff predicts
                        // `stop_gained` while this rather is a intron variant.  We skip this variant.
                        record_csqs.contains(&"stop_gained")
                            && record.var == "3-193366573-A-ATATTGCCTAGAATGAACT",
                        // For `GRCh37:3:193409913:ATAAAT:A`, there appears to be a model error
                        // in SnpEff as it predicts `exon_loss`.  We skip this variant.
                        record_csqs.contains(&"exon_loss_variant")
                            && record.var == "3-193409913-ATAAAT-A",
                        // SnpEff may predict `pMet1.?` as `initiator_codon_variant` rather than `start_lost`.
                        expected_one_of.contains(&"start_lost")
                            && (record_csqs.contains(&"initiator_codon_variant")),
                        // Similarly, SnpEff may predict `c.-1_1` as `start_retained` rather than `start_lost`.
                        expected_one_of.contains(&"start_lost")
                            && (record_csqs.contains(&"start_retained_variant")),
                        // SnpEff calls this insertion at c.5193+2_5193+3insT a splice donor variant
                        // even though the third intronic base is affected, not the first or second
                        record_csqs.contains(&"splice_donor_variant")
                            && expected_one_of.contains(&"splice_region_variant")
                            && [
                                "17-41215347-T-TA",
                                "17-41215888-T-TA",
                                "17-41242958-T-TA",
                                "17-41256882-T-TA",
                                "17-41276031-T-TA",
                                "17-41277285-T-TA",
                            ]
                            .contains(&record.var.as_str()),
                        // we call exonic_splice_region_variant, while the others only call splice_region_variant
                        record_csqs.contains(&"splice_region_variant")
                            && expected_one_of.contains(&"exonic_splice_region_variant"),
                    ]
                    .iter()
                    .any(|b| *b);

                    assert!(
                        found_one,
                        "line no. {}, variant: {}, tx: {}, hgvs_c: {:?}, hgvs_p: {:?}, \
                        their_csqs: {:?}, expected_one_of: {:?}, our_csqs: {:?}",
                        lineno,
                        record.var,
                        record.tx,
                        ann.hgvs_c.as_ref(),
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

    #[test]
    fn test_predict_multiple_bounds_checking() -> Result<(), anyhow::Error> {
        let alt_seq = String::from("ATGCGTACGTAGCTAGCT");
        let n_min = 5;
        let n_max = 7;

        let total_delta_negative = -10;
        let new_length_negative = (n_max - n_min + 1) + total_delta_negative;
        assert!(
            new_length_negative < 0,
            "Length should be negative and caught by the guard"
        );

        let total_delta_past_end = 100;
        let new_length = (n_max - n_min + 1) + total_delta_past_end;
        let start_idx = (n_min - 1) as usize;
        let end_idx = start_idx + new_length as usize;

        assert!(
            end_idx > alt_seq.len(),
            "End index should exceed alt_seq length and be caught by the guard"
        );

        let total_delta_valid = 2; // e.g., del 1 base, ins 3 bases
        let new_length_valid = (n_max - n_min + 1) + total_delta_valid;
        let end_idx_valid = start_idx + new_length_valid as usize;

        assert!(new_length_valid >= 0);
        assert!(end_idx_valid <= alt_seq.len());
        let _new_substring = &alt_seq[start_idx..end_idx_valid];

        Ok(())
    }
}

#[derive(
    Debug,
    Default,
    Clone,
    Copy,
    PartialEq,
    Eq,
    clap::ValueEnum,
    parse_display::FromStr,
    parse_display::Display,
)]
#[display(style = "kebab-case")]
pub enum SequenceReporting {
    #[default]
    None,
    Reference,
    Alternative,
    Both,
}

impl SequenceReporting {
    pub fn includes_ref(&self) -> bool {
        matches!(self, Self::Reference | Self::Both)
    }
    pub fn includes_alt(&self) -> bool {
        matches!(self, Self::Alternative | Self::Both)
    }
}
