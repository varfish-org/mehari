//! Compute molecular consequence of variants.
use super::{
    ann::{Allele, AnnField, Consequence, FeatureBiotype, FeatureType, Pos, Rank, SoFeature},
    provider::Provider as MehariProvider,
};
use crate::annotate::cli::{ConsequenceBy, TranscriptSource};
use crate::pbs::txs::{GenomeAlignment, Strand, TranscriptBiotype, TranscriptTag};
use enumflags2::BitFlags;
use hgvs::mapper::altseq::AltSeqBuilder;
use hgvs::parser::{NoRef, ProteinEdit};
use hgvs::{
    data::interface::{Provider, TxForRegionRecord},
    mapper::{assembly, Error},
    parser::{
        Accession, CdsFrom, GenomeInterval, GenomeLocEdit, HgvsVariant, Mu, NaEdit, ProtLocEdit,
    },
};
use itertools::Itertools;
use std::cmp::Ordering;
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
    pub(crate) provider: Arc<MehariProvider>,
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

pub type Consequences = BitFlags<Consequence>;

impl ConsequencePredictor {
    pub fn new(provider: Arc<MehariProvider>, config: Config) -> Self {
        tracing::info!("Building transcript interval trees ...");
        let chrom_to_acc = provider.build_chrom_to_acc(None);

        let reference_available = provider.reference_available();

        let mapper_config = assembly::Config {
            assembly: provider.assembly(),
            replace_reference: reference_available,
            strict_bounds: false,
            renormalize_g: reference_available,
            genome_seq_available: reference_available,
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
        let mut norm_var = self.normalize_variant(var);

        // TODO check for VCF specification version.
        // According to VCF specification (>=4.1), an alternative of "N" means REF=ALT
        // Prior to 4.1, it indicated a deletion.
        if norm_var.alternative == "N" {
            norm_var.alternative = norm_var.reference.clone();
        }

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

        // We follow hgvs conventions and therefore normalize input variants
        let var_g = Self::get_var_g(&norm_var, chrom_acc);
        let var_g = if self.mapper.config.renormalize_g {
            self.mapper.maybe_normalize(&var_g)?
        } else {
            var_g
        };

        // Get all affected transcripts.
        let (var_start, var_end) = match &var_g {
            HgvsVariant::GenomeVariant { loc_edit, .. } => {
                let loc = loc_edit.loc.inner();
                let edit = loc_edit.edit.inner();
                let start = loc
                    .start
                    .map(|s| s - 1)
                    .expect("Failed to get start position");
                let end = loc.end.map(|s| s - 1).expect("Failed to get end position");
                let ref_len = norm_var.reference.len() as i32;
                let old_end = start + ref_len;
                let end = if edit.is_ins() || edit.is_dup() {
                    start
                } else {
                    end + 1
                };
                assert_eq!(end, old_end, "End position mismatch, start is {}, end is {}, ref_len is {}, var_g is {}, norm_var is {:?}", start, end, ref_len, &var_g, &norm_var);
                (start, end)
            }
            _ => unreachable!(),
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
            let hgvs_g = format!("{}", &NoRef(&var_g));
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
                rank: None,
                distance: None,
                strand: 0,
                hgvs_g,
                hgvs_c: None,
                hgvs_p: None,
                cdna_pos: None,
                cds_pos: None,
                protein_pos: None,
                gene_symbol: "".to_string(),
                messages: None,
            }])));
        }

        // Compute annotations for all (picked) transcripts first, skipping `None`` results.
        let anns_all_txs = txs
            .into_iter()
            .map(|tx| self.build_ann_field(var, var_g.clone(), tx, var_start, var_end))
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

        /// Extract the first (i.e. most severe) consequence per group.
        fn first_csq_per_group(grouping: HashMap<String, Vec<AnnField>>) -> Vec<AnnField> {
            // Sort by `Consequence` and `FeatureBiotype`
            //
            // This uses the invariant that the consequences in the ANN fields are sorted already
            // and there is at least one consequence.
            grouping
                .into_values()
                .map(|mut anns| {
                    anns.sort_by_key(|ann| {
                        (ann.consequences[0], biotype_order(&ann.feature_biotype))
                    });
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

    fn build_ann_field(
        &self,
        orig_var: &VcfVariant,
        var_g: HgvsVariant,
        tx_record: TxForRegionRecord,
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

        // tracing::trace!("Annotating variant {:?} w.r.t. transcript {}", orig_var, &tx_record.tx_ac);

        assert_eq!(
            tx.genome_alignments.len(),
            1,
            "At this point, only one genome alignment is expected"
        );

        // TODO: Report selenocysteine modifications
        // let is_seleno = tx.tags.contains(&(TranscriptTag::Selenoprotein as i32));

        let alignment = tx.genome_alignments.first().unwrap();
        let strand = Strand::try_from(alignment.strand).expect("invalid strand");

        let mut consequences = Consequences::empty();

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
                        let dist_start: i32 = var_start + 1 - intron_start;
                        let dist_end: i32 = -(intron_end + 1 - var_end);
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

            if var_overlaps(exon_start, exon_end) {
                let consequences_exonic = Self::analyze_exonic_variant(
                    strand, var_start, var_end, exon_start, exon_end, &rank,
                );
                consequences |= consequences_exonic;
            }

            if let Some(intron_start) = intron_start {
                let consequences_intronic = Self::analyze_intronic_variant(
                    &var_g,
                    alignment,
                    strand,
                    var_start,
                    var_end,
                    intron_start,
                    intron_end,
                );
                consequences |= consequences_intronic;
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
                // unreachable!("variant is intergenic but this cannot happen here");
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
                // unreachable!("variant is intergenic but this cannot happen here");
            }
            if distance.is_none() {
                distance = Some(val);
            }
        }

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

                    let conservative = is_conservative_cds_variant(&var_c);

                    let consequences_cds =
                        Self::analyze_cds_variant(&var_c, is_exonic, is_intronic, conservative);

                    // Analyze `var_p` for changes in the protein sequence.
                    let consequences_protein = self.analyze_protein_variant(
                        &var_c,
                        &var_p,
                        &protein_pos,
                        conservative,
                        &tx_record.tx_ac,
                    );

                    consequences |= consequences_cds | consequences_protein;

                    self.consequences_fix_special_cases(
                        &mut consequences,
                        consequences_cds,
                        consequences_protein,
                        &var_g,
                        &var_n,
                        &var_c,
                        &var_p,
                    );

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
                "No consequences for {:?} on {} (hgvs_p={}) - adding `gene_variant`; \
                most likely the transcript has multiple stop codons and the variant \
                lies behind the first.",
                orig_var,
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

        let hgvs_g = format!("{}", &NoRef(&var_g));
        let hgvs_g = Some(hgvs_g.split(':').nth(1).unwrap().to_owned());

        // Build and return ANN field from the information derived above.
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
            rank,
            hgvs_g,
            hgvs_c: hgvs_t,
            hgvs_p,
            cdna_pos: tx_pos,
            cds_pos,
            protein_pos,
            strand,
            distance,
            messages: None,
        }))
    }

    #[allow(clippy::too_many_arguments, unused_variables)]
    fn consequences_fix_special_cases(
        &self,
        consequences: &mut Consequences,
        consequences_cds: Consequences,
        consequences_protein: Consequences,
        var_g: &HgvsVariant,
        var_n: &HgvsVariant,
        var_c: &HgvsVariant,
        var_p: &HgvsVariant,
    ) {
        // If we have a transcript_ablation, we can remove all other consequences
        if consequences.contains(Consequence::TranscriptAblation) {
            *consequences = Consequence::TranscriptAblation.into();
            return;
        }

        // If a frameshift/ins/del was predicted on the CDS level,
        // but any relevant consequence (i.e. not just GeneVariant) was produced on the protein level,
        // then it is likely that the frameshift induced a more specific consequence.
        let check_cds_csqs: Consequences = Consequence::FrameshiftVariant.into();
        // | Consequence::DisruptiveInframeDeletion
        // | Consequence::ConservativeInframeDeletion
        // | Consequence::DisruptiveInframeInsertion
        // | Consequence::DisruptiveInframeInsertion;
        let checked = consequences_cds & check_cds_csqs;
        if checked != Consequences::empty()
            // if the protein consequence is not effectively empty, we remove the CDS frameshift consequence
            && !(consequences_protein.eq(&Consequence::GeneVariant)
            || consequences_protein.is_empty())
            // if the protein consequence also includes a frameshift, then we keep it
            && !consequences_protein
            .intersects(Consequence::FrameshiftElongation | Consequence::FrameshiftTruncation | Consequence::FrameshiftVariant)
        {
            *consequences &= !checked;
        }

        // In some cases, we predict a stop lost based on the cds variant
        // but the protein translation does not confirm this.
        //
        // e.g.:
        // 20:35511609:CAAGCCGCCTCCAGGTAGCAGCCACAGCCAGGAGCACACAGACAGAAGACTGTGTCATGGGTCATGGCCCCTCCGCACACCTACAGGTTTGCCAAAGGAA:C
        if consequences_cds.contains(Consequence::StopLost)
            && !consequences_protein
                .intersects(Consequence::StopLost | Consequence::ProteinAlteringVariant)
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
        if consequences.contains(Consequence::StartLost) {
            if let (
                HgvsVariant::TxVariant {
                    loc_edit: n_loc_edit,
                    accession,
                    ..
                },
                HgvsVariant::CdsVariant {
                    loc_edit: c_loc_edit,
                    ..
                },
            ) = (var_n, var_c)
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
                            tracing::trace!(
                                "Fixing StartLost → StartRetained for {}, {}, {}",
                                &var_g,
                                &var_c,
                                &var_p
                            );
                            *consequences &= !Consequence::StartLost;
                            *consequences |= Consequence::StartRetainedVariant;
                        }
                    }
                }
            }
        }
        match var_c {
            HgvsVariant::CdsVariant { loc_edit, .. } => {
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
            _ => unreachable!(),
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
    ) -> Consequences {
        let mut consequences: Consequences = Consequences::empty();

        let var_overlaps =
            |start: i32, end: i32| -> bool { overlaps(var_start, var_end, start, end) };

        // For insertions, we need to consider the case of the insertion being right at
        // the exon/intron junction.  We can express this with a shift of 1 for using
        // "< / >" X +/- shift and meaning <= / >= X.
        let ins_shift = match var_g {
            HgvsVariant::GenomeVariant {
                loc_edit: GenomeLocEdit { edit, .. },
                ..
            } => {
                let edit = edit.inner();
                if edit.is_ins() || edit.is_dup() {
                    1
                } else {
                    0
                }
            }
            _ => unreachable!(),
        };

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

    fn analyze_cds_variant(
        var_c: &HgvsVariant,
        is_exonic: bool,
        is_intronic: bool,
        conservative: bool,
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
            if starts_left_of_stop && ends_right_of_stop {
                consequences |= Consequence::StopLost;
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
                                let original_sequence = reference_data.aa_sequence.clone();
                                let alt_data = AltSeqBuilder::new(var_c.clone(), reference_data)
                                    .build_altseq()
                                    .unwrap();
                                let alt_data = alt_data.first().unwrap();
                                let altered_sequence = Some(alt_data.aa_sequence.clone());
                                if let Some(ref altered_sequence) = altered_sequence {
                                    let original_sequence = &original_sequence;

                                    // trim altered sequence to the first stop encountered
                                    let altered_sequence = if let Some(pos) =
                                        altered_sequence.find('*')
                                    // do not use the 'X' fallback here,
                                    // as that is _usually_ only added
                                    // when the number of bases is not divisible by 3.
                                    // We only want to identify cases where a new/later
                                    // stop codon is encountered
                                    // .or_else(|| altered_sequence.find('X'))
                                    {
                                        &altered_sequence[..=pos]
                                    } else {
                                        altered_sequence
                                    };

                                    match altered_sequence.len().cmp(&original_sequence.len()) {
                                        Ordering::Less => {
                                            consequences |= Consequence::FrameshiftTruncation;
                                        }
                                        Ordering::Equal => {
                                            consequences |= Consequence::MissenseVariant;
                                            // TODO: discuss stop_retained
                                            // consequences |= Consequence::StopRetainedVariant;
                                            consequences &= !Consequence::FrameshiftVariant;
                                        }
                                        Ordering::Greater => {
                                            consequences |= Consequence::FrameshiftElongation;
                                        }
                                    }
                                }
                            }
                        }
                        ProteinEdit::Ext { .. } => {
                            consequences |= Consequence::StopLost;
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
                                    if let Some(p) = protein_pos {
                                        if p.total.is_some_and(|t| p.ord == t - 1) && conservative {
                                            consequences &= !Consequence::StopGained;
                                            consequences |=
                                                Consequence::ConservativeInframeDeletion;
                                        }
                                    }
                                }
                            } else {
                                consequences |= Consequence::MissenseVariant;
                                // Missense variants that affect selenocysteine are marked
                                // as rare amino acid variants.
                                if alternative.contains("U")
                                    || (loc.start == loc.end) && loc.start.aa == "U"
                                {
                                    consequences |= Consequence::RareAminoAcidVariant;
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
    use crate::annotate::cli::{TranscriptPickType, TranscriptSettings};
    use crate::annotate::seqvars::provider::ConfigBuilder as MehariProviderConfigBuilder;
    use crate::annotate::seqvars::{
        load_tx_db, run_with_writer, Args, AsyncAnnotatedVariantWriter, PathOutput,
    };
    use crate::common::noodles::{open_variant_reader, open_variant_writer, NoodlesVariantReader};
    use csv::ReaderBuilder;
    use futures::TryStreamExt;
    use pretty_assertions::assert_eq;
    use serde::Deserialize;
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
    #[case("3:193311167:ATGT:T", vec![Consequence::StartLost, Consequence::ConservativeInframeDeletion])]
    #[case("3:193311170:TGGC:C", vec![Consequence::ConservativeInframeDeletion])]
    #[case("3:193311170:TGGCG:G", vec![Consequence::FrameshiftVariant, Consequence::FrameshiftTruncation])]
    #[case("3:193311180:GTCG:G", vec![Consequence::DisruptiveInframeDeletion])]
    #[case("3:193409910:GAAA:G", vec![Consequence::ConservativeInframeDeletion])]
    #[case("3:193409913:ATAA:A", vec![Consequence::StopLost, Consequence::FeatureElongation, Consequence::ConservativeInframeDeletion])]
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
        let expected_vcf = "tests/data/annotate/seqvars/vep.disagreement-cases.expected.vcf";
        let output = NamedTempFile::new()?;
        let mut writer = open_variant_writer(output.as_ref()).await?;
        run_with_writer(
            &mut writer,
            &Args {
                reference: None,
                in_memory_reference: true,
                genome_release: None,
                path_input_ped: None,
                path_input_vcf: path_input_vcf.into(),
                output: PathOutput {
                    path_output_vcf: Some(output.as_ref().to_str().unwrap().into()),
                    path_output_tsv: None,
                },
                transcript_settings: TranscriptSettings {
                    report_most_severe_consequence_by: Some(ConsequenceBy::Allele),
                    pick_transcript: vec![TranscriptPickType::ManeSelect],
                    ..Default::default()
                },
                max_var_count: None,
                hgnc: None,
                sources: crate::annotate::seqvars::Sources {
                    transcripts: Some(vec![tx_path.into()]),
                    frequencies: None,
                    clinvar: None,
                },
            },
        )
        .await?;
        writer.shutdown().await?;

        let records_written = read_vcf(output).await?;
        let records_expected = read_vcf(expected_vcf).await?;
        assert_eq!(records_expected, records_written);

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
        with_reference: bool
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
            if record.var == "17-41196310-GGTGGAAGTGTTTGCTACCAAGTTTATTTGCAGTGTTAACAGCACAACATTTACAAAACGTATTTTGTACAATCAAGTCTTCACTGCCCTTGCACACTGGGGGGGCTAGGGAAGACCTAGTCCTTCCAACAGCTATAAACAGTCCTGGATAATGGGTTTATGAAAAACACTTTTTCTTCCTTCAGCAAGCAAAATTATTTATGAAGCTGTATGGTTTCAGCAACAGGGAGCAAAGGAAAAAAATCACCTCAAAGAAAGCAACAGCTTCCTTCCTGGTGGGATCTGTCATTTTATAGATATGAAATATTCATGCCAGAGGTCTTATATTTTAAGAGGAATGGATTATATACCAGAGCTACAACAATAAACATTTTACTTATTACTAATGAGGAATTAGAAGACTGTCTTTGGAAACCGGTTCTTGAAAATCTTCTGCTGTTTTAGAACACATTCTTTAGAAATCTAGCAAATATATCTCAGACTTTTAGAAATCTCTTCTAGTTTCATTTTCCTTTTTTTTTTTTTTTTTTTGAGCCACAGTCTCACTGTCACCCAGGCTGGAGTGCCGTGGTATGATCTTGGCTCACTGCAACCTCCACCTCCCGGGCTGAAGTGATTCTCCTGCCTTAGCCACCTGAGTAGCTGGGATTACAGGTGTCCACCACCATGACCGGCTAATTTCTGTATTTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTTTCGAACTCCTGACCTCCAGTGATCTGCCCACCTTGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCATGCCCAGGTTTCAAGTTTCCTTTTCATTTCTAATACCTGCCTCAGAATTTCCTCCCCAATGTTCCACTCCAACATTTGAGAACTGCCCAAGGACTATTCTGACTTTAAGTCACATAATCGATCCCAAGCACTCTCCTTCCATTGAAGGGTCTGACTCTCTGCCTTTGTGAACACAGGGTTTTAGAGAAGTAAACTTAGGGAAACCAGCTATTCTCTTGAGGCCAAGCCACTCTGTGCTTCCAGCCCTAAGCCAACAACAGCCTGAATAGAAAGAATAGGGCTGATAAATAATGAATCAGCATCTTGCTCAATTGGTGGCGTTTAAATGGTTTTAAAATCTTCTCAGGTGAAAAATTACCATAATTTTGTGCTCATGGCAGATTTCCAAGGGAGACTTCAAGCAGAAAATCTTTAAGGGACCCTTGCATAGCCAGAAGTCCTTTTCAGGCTGATGTACATAAAATATTTAGTAGCCAGGACAGTAGAAGGACTGAAGAGTGAGAGGAGCTCCCAGGGCCTGGAAAGGCCACTTTGTAAGCTCATTCTTGGGGTCCTGTGGCTCTGTACCTGTGGCTGGCTGCAGTCAGTAGTGGCTGTGGGGGATCTGGGGTATCAGGTAGGTGTCCAGCTCCTGGCACTGGTAGAGTGCTACACTGTCCAACACCCACTCTCGGGTCACCACAGGTGCCTCACACATCTGCCCAATT-G" {
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
                        expected_one_of.contains(&"splice_region_variant")
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
                    ]
                    .iter()
                    .any(|b| *b);

                    assert!(
                        found_one,
                        "line no. {}, variant: {}, tx: {}, hgvs_c: {:?}, hgvs_p: {:?}, \
                        record_csqs: {:?}, expected_one_of: {:?}, ann_csqs: {:?}",
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
}
