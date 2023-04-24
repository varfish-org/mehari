//! Compute molecular consequence of variants.

use std::{collections::HashMap, rc::Rc};

use hgvs::{
    data::interface::{Provider, TxForRegionRecord},
    mapper::assembly::{Config as AssemblyConfig, Mapper as AssemblyMapper},
    parser::{
        Accession, CdsFrom, GenomeInterval, GenomeLocEdit, HgvsVariant, Mu, NaEdit, ProtLocEdit,
    },
    static_data::Assembly,
};

use crate::db::create::txs::data::{Strand, TranscriptBiotype};

use super::{
    ann::{Allele, AnnField, Consequence, FeatureBiotype, FeatureType, Pos, Rank, SoFeature},
    provider::MehariProvider,
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

pub struct ConsequencePredictor {
    /// The internal transcript provider for locating transcripts.
    provider: Rc<MehariProvider>,
    /// Assembly mapper for variant consequence prediction.
    mapper: AssemblyMapper,
    /// Mapping from chromosome name to accession.
    chrom_to_acc: HashMap<String, String>,
}

/// Padding to look for genes upstream/downstream.
pub const PADDING: i32 = 5_000;
/// Generally used alternative alignment method.
pub const ALT_ALN_METHOD: &str = "splign";

impl ConsequencePredictor {
    pub fn new(provider: Rc<MehariProvider>, assembly: Assembly) -> Self {
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

        let config = AssemblyConfig {
            replace_reference: false,
            strict_bounds: false,
            renormalize_g: false,
            genome_seq_available: false,
            ..Default::default()
        };
        let mapper = AssemblyMapper::new(config, provider.clone());

        ConsequencePredictor {
            provider,
            mapper,
            chrom_to_acc,
        }
    }

    pub fn predict(&self, var: &VcfVariant) -> Result<Option<Vec<AnnField>>, anyhow::Error> {
        // Normalize variant by stripping common prefix and suffix.
        let norm_var = self.normalize_variant(var);

        // Obtain accession from chromosome name.
        let chrom_acc = self.chrom_to_acc.get(&norm_var.chromosome);
        let chrom_acc = if let Some(chrom_acc) = chrom_acc {
            chrom_acc
        } else {
            tracing::debug!(
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
        let mut txs =
            self.provider
                .get_tx_for_region(chrom_acc, ALT_ALN_METHOD, qry_start, qry_end)?;
        txs.sort_by(|a, b| a.tx_ac.cmp(&b.tx_ac));

        // Generate `AnnField` records for each transcript.
        //
        // Skip `None` results.
        Ok(Some(
            txs.into_iter()
                .map(|tx| {
                    self.build_ann_field(var, &norm_var, tx, chrom_acc.clone(), var_start, var_end)
                })
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .flatten()
                .collect::<Vec<_>>(),
        ))
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

        let tx = self.provider.get_tx(&tx_record.tx_ac).unwrap();

        // Skip transcripts that are protein coding but do not have a CDS.
        // TODO: do not include such transcripts when building the database.
        if TranscriptBiotype::from_i32(tx.biotype).expect("invalid tx biotype")
            == TranscriptBiotype::Coding
            && tx.start_codon.is_none()
        {
            return Ok(None);
        }

        let mut consequences: Vec<Consequence> = Vec::new();

        let alignment = tx.genome_alignments.first().unwrap();

        let mut min_start = None;
        let mut max_end = None;

        // Find first exon that overlaps with variant or intron that contains the variant.
        //
        // Note that exons are stored in genome position order.
        let mut prev_end = None;
        let mut rank = Rank::default();
        let mut is_exonic = false;
        let mut is_intronic = false;
        let mut distance = None;
        let mut tx_len = 0;
        for exon_alignment in &alignment.exons {
            tx_len += exon_alignment.alt_end_i - exon_alignment.alt_start_i;

            let exon_start = exon_alignment.alt_start_i;
            let exon_end = exon_alignment.alt_end_i;
            let intron_start = prev_end;
            let intron_end = exon_start;

            // Check the cases where the variant overlaps with the exon or is contained within an
            // intron.
            if var_start < exon_end && exon_start < var_end {
                // overlaps with exon
                rank = Rank {
                    ord: exon_alignment.ord + 1,
                    total: alignment.exons.len() as i32,
                };
                is_exonic = true;

                if var_start <= exon_start || var_end >= exon_end {
                    // Overlaps with exon/intron boundary.
                    distance = Some(0);
                } else {
                    let dist_start = -(var_start - exon_start + 1);
                    let dist_end = exon_end - var_end + 1;
                    if dist_end >= dist_start.abs() {
                        distance = Some(dist_end);
                    } else {
                        distance = Some(dist_start);
                    }
                }
            } else if let Some(intron_start) = intron_start {
                if var_start >= intron_start && var_end <= intron_end {
                    // Contained within intron: cannot be in next exon.
                    if !is_exonic {
                        rank = Rank {
                            ord: exon_alignment.ord + 1,
                            total: alignment.exons.len() as i32 - 1,
                        };
                        is_intronic = true;

                        let dist_start = -(var_start - intron_start);
                        let dist_end = intron_end - var_end;
                        if dist_end <= dist_start.abs() {
                            distance = Some(dist_end);
                        } else {
                            distance = Some(dist_start);
                        }
                    }
                }
            }

            // Check the cases where the variant overlaps with whole exon.
            if var_start <= exon_start && var_end >= exon_end {
                consequences.push(Consequence::ExonLossVariant);
                if var_start < exon_start {
                    if Strand::from_i32(alignment.strand).expect("invalid strand") == Strand::Plus
                        && rank.ord != 1
                    {
                        consequences.push(Consequence::SpliceAcceptorVariant);
                    } else if Strand::from_i32(alignment.strand).expect("invalid strand")
                        == Strand::Minus
                        && rank.ord != rank.total
                    {
                        consequences.push(Consequence::SpliceDonorVariant);
                    }
                }
                if var_end > exon_end {
                    if Strand::from_i32(alignment.strand).expect("invalid strand") == Strand::Plus
                        && rank.ord != rank.total
                    {
                        consequences.push(Consequence::SpliceDonorVariant);
                    } else if Strand::from_i32(alignment.strand).expect("invalid strand")
                        == Strand::Minus
                        && rank.ord != rank.total
                    {
                        consequences.push(Consequence::SpliceAcceptorVariant);
                    }
                }
            }
            if let Some(intron_start) = intron_start {
                // For insertions, we need to consider the case of the insertion being right at
                // the exon/intron junction.  We can express this with a shift of 1 for using
                // "</>" X +/- shift and meaning <=/>= X.
                let ins_shift = if var.reference.is_empty() { 1 } else { 0 };

                // Check the cases where the variant overlaps with the splice acceptor/donor site.
                if var_start < intron_start + 2 && var_end > intron_start - ins_shift {
                    // Left side, is acceptor/donor depending on transcript's strand.
                    match Strand::from_i32(alignment.strand).expect("invalid strand") {
                        Strand::Plus => consequences.push(Consequence::SpliceDonorVariant),
                        Strand::Minus => consequences.push(Consequence::SpliceAcceptorVariant),
                    }
                }
                // Check the case where the variant overlaps with the splice donor site.
                if var_start < intron_end + ins_shift && var_end > intron_end - 2 {
                    // Left side, is acceptor/donor depending on transcript's strand.
                    match Strand::from_i32(alignment.strand).expect("invalid strand") {
                        Strand::Plus => consequences.push(Consequence::SpliceAcceptorVariant),
                        Strand::Minus => consequences.push(Consequence::SpliceDonorVariant),
                    }
                }
            }
            // Check the case where the variant overlaps with the splice region (1-3 bases in exon
            // or 3-8 bases in intron).  We have to check all cases independently and not with `else`
            // because the variant may be larger.
            if let Some(intron_start) = intron_start {
                if (var_start < intron_start + 8 && var_end > intron_start + 2)
                    || (var_start < intron_end - 8 && var_end > intron_end - 2)
                {
                    consequences.push(Consequence::SpliceRegionVariant);
                }
                if var_start < exon_end && var_end > exon_end - 3 {
                    if Strand::from_i32(alignment.strand).expect("invalid strand") == Strand::Plus {
                        if rank.ord != rank.total {
                            consequences.push(Consequence::SpliceRegionVariant);
                        }
                    } else {
                        // alignment.strand == Strand::Minus
                        if rank.ord != 1 {
                            consequences.push(Consequence::SpliceRegionVariant);
                        }
                    }
                }
                if var_start < exon_start + 3 && var_end > exon_start {
                    if Strand::from_i32(alignment.strand).expect("invalid strand") == Strand::Plus {
                        if rank.ord != 1 {
                            consequences.push(Consequence::SpliceRegionVariant);
                        }
                    } else {
                        // alignment.strand == Strand::Minus
                        if rank.ord != rank.total {
                            consequences.push(Consequence::SpliceRegionVariant);
                        }
                    }
                }
            }

            min_start = Some(std::cmp::min(min_start.unwrap_or(exon_start), exon_start));
            max_end = Some(std::cmp::max(max_end.unwrap_or(exon_end), exon_end));

            prev_end = Some(exon_end);
        }

        let min_start = min_start.expect("must have seen exon");
        let max_end = max_end.expect("must have seen exon");

        let feature_biotype =
            match TranscriptBiotype::from_i32(tx.biotype).expect("invalid transcript biotype") {
                TranscriptBiotype::Coding => FeatureBiotype::Coding,
                TranscriptBiotype::NonCoding => FeatureBiotype::Noncoding,
            };

        let is_upstream = var_end <= min_start;
        let is_downstream = var_start >= max_end;
        if is_exonic {
            if !feature_biotype.is_coding() {
                consequences.push(Consequence::NonCodingTranscriptExonVariant);
            }
        } else if is_intronic {
            if !feature_biotype.is_coding() {
                consequences.push(Consequence::NonCodingTranscriptIntronVariant);
            } else {
                consequences.push(Consequence::IntronVariant);
            }
        } else if is_upstream {
            let val = -(min_start - var_end);
            if val.abs() <= 5_000 {
                match Strand::from_i32(alignment.strand).expect("invalid strand") {
                    Strand::Plus => consequences.push(Consequence::UpstreamGeneVariant),
                    Strand::Minus => consequences.push(Consequence::DownstreamGeneVariant),
                }
            }
            distance = Some(val);
        } else if is_downstream {
            let val = var_start - max_end;
            if val.abs() <= 5_000 {
                match Strand::from_i32(alignment.strand).expect("invalid strand") {
                    Strand::Plus => consequences.push(Consequence::DownstreamGeneVariant),
                    Strand::Minus => consequences.push(Consequence::UpstreamGeneVariant),
                }
            }
            distance = Some(val);
        }

        let var_g = HgvsVariant::GenomeVariant {
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
        };

        let (rank, hgvs_t, hgvs_p, tx_pos, cds_pos, protein_pos) = if !is_upstream && !is_downstream
        {
            // Gracefully handle problems in the projection (in this case "Non-adjacent exons for ...").
            // TODO: do not include such transcripts when building the tx database.
            let var_n = self.mapper.g_to_n(&var_g, &tx.id).map_or_else(
                |e| {
                    if e.to_string().contains("Non-adjacent exons for") {
                        Ok(None)
                    } else {
                        Err(e)
                    }
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

            let (var_t, _var_p, hgvs_p, cds_pos, protein_pos) = match feature_biotype {
                FeatureBiotype::Coding => {
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
                    let mut hgvs_p = hgvs_p.split(':').nth(1).unwrap().to_owned();
                    for (aa3, aa1) in hgvs::sequences::AA3_TO_AA1_VEC.iter() {
                        hgvs_p = hgvs_p.replace(aa3, aa1);
                    }
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
                                consequences.push(Consequence::StartLost)
                            }
                            // stop codon
                            let starts_left_of_stop = start_cds_from == CdsFrom::Start;
                            let ends_right_of_stop = end_cds_from == CdsFrom::End;
                            if starts_left_of_stop && ends_right_of_stop {
                                consequences.push(Consequence::StopLost)
                            }

                            // Detect variants affecting the 5'/3' UTRs.
                            if start_cds_from == CdsFrom::Start {
                                if start_base < 0 {
                                    consequences.push(Consequence::FivePrimeUtrVariant);
                                }
                            } else if end_cds_from == CdsFrom::End {
                                consequences.push(Consequence::ThreePrimeUtrVariant);
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
                                        consequences.push(Consequence::FrameshiftVariant);
                                    }
                                    hgvs::parser::ProteinEdit::Ext { .. } => {
                                        consequences.push(Consequence::StopLost);
                                        consequences.push(Consequence::FeatureElongation);
                                    }
                                    hgvs::parser::ProteinEdit::Subst { alternative } => {
                                        if alternative.is_empty() {
                                            consequences.push(Consequence::SynonymousVariant);
                                        } else if is_stop(alternative) {
                                            if loc.start == loc.end && is_stop(&loc.start.aa) {
                                                consequences.push(Consequence::StopRetainedVariant);
                                            } else {
                                                consequences.push(Consequence::StopGained);
                                            }
                                        } else {
                                            consequences.push(Consequence::MissenseVariant);
                                        }
                                    }
                                    hgvs::parser::ProteinEdit::DelIns { alternative } => {
                                        if conservative {
                                            consequences
                                                .push(Consequence::ConservativeInframeDeletion);
                                        } else {
                                            consequences
                                                .push(Consequence::DisruptiveInframeDeletion);
                                        }
                                        if alternative.contains('*')
                                            || alternative.contains('X')
                                            || alternative.contains("Ter")
                                        {
                                            consequences.push(Consequence::StopGained);
                                        }
                                    }
                                    hgvs::parser::ProteinEdit::Ins { .. }
                                    | hgvs::parser::ProteinEdit::Dup => {
                                        if conservative {
                                            consequences
                                                .push(Consequence::ConservativeInframeInsertion);
                                        } else {
                                            consequences
                                                .push(Consequence::DisruptiveInframeInsertion);
                                        }
                                        consequences
                                            .push(Consequence::ConservativeInframeInsertion);
                                    }
                                    hgvs::parser::ProteinEdit::Del => {
                                        if conservative {
                                            consequences
                                                .push(Consequence::ConservativeInframeDeletion);
                                        } else {
                                            consequences
                                                .push(Consequence::DisruptiveInframeDeletion);
                                        }
                                    }
                                    hgvs::parser::ProteinEdit::Ident => {
                                        consequences.push(Consequence::SynonymousVariant)
                                    }
                                };
                            }
                            ProtLocEdit::NoChange | ProtLocEdit::NoChangeUncertain => {
                                consequences.push(Consequence::SynonymousVariant)
                            }
                            ProtLocEdit::InitiationUncertain => {
                                consequences.push(Consequence::StartLost)
                            }
                            ProtLocEdit::NoProtein
                            | ProtLocEdit::NoProteinUncertain
                            | ProtLocEdit::Unknown => (),
                        },
                        _ => panic!("Must be protein variant: {}", &var_p),
                    }

                    (var_c, Some(var_p), hgvs_p, cds_pos, protein_pos)
                }
                FeatureBiotype::Noncoding => (var_n, None, None, None, None),
            };
            let hgvs_t = format!("{}", &var_t);
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
        consequences.sort();
        consequences.dedup();
        if consequences.is_empty() {
            tracing::debug!(
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
            consequences.push(Consequence::GeneVariant);
        }
        let putative_impact = (*consequences.first().unwrap()).into();

        // Build and return ANN field from the information derived above.
        Ok(Some(AnnField {
            allele: Allele::Alt {
                alternative: orig_var.alternative.clone(),
            },
            consequences,
            putative_impact,
            gene_symbol: tx.gene_name.clone(),
            gene_id: format!("HGNC:{}", &tx.gene_id),
            feature_type: FeatureType::SoTerm {
                term: SoFeature::Transcript,
            },
            feature_id: tx.id.clone(),
            feature_biotype,
            rank,
            hgvs_t,
            hgvs_p,
            tx_pos,
            cds_pos,
            protein_pos,
            distance,
            messages: None,
        }))
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

#[cfg(test)]
mod test {
    use std::{fs::File, io::BufReader};

    use csv::ReaderBuilder;
    use pretty_assertions::assert_eq;
    use serde::Deserialize;

    use crate::annotate::seqvars::{ann::PutativeImpact, load_tx_db};

    use super::*;

    #[test]
    fn annotate_snv_brca1_one_variant() -> Result<(), anyhow::Error> {
        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Rc::new(MehariProvider::new(tx_db, Assembly::Grch37p10));

        let predictor = ConsequencePredictor::new(provider, Assembly::Grch37p10);

        let res = predictor
            .predict(&VcfVariant {
                chromosome: String::from("17"),
                position: 41_197_701,
                reference: String::from("G"),
                alternative: String::from("C"),
            })?
            .unwrap();

        assert_eq!(res.len(), 4);
        assert_eq!(
            res[0],
            AnnField {
                allele: Allele::Alt {
                    alternative: String::from("C")
                },
                consequences: vec![Consequence::MissenseVariant,],
                putative_impact: PutativeImpact::Moderate,
                gene_symbol: String::from("BRCA1"),
                gene_id: String::from("HGNC:1100"),
                feature_type: FeatureType::SoTerm {
                    term: SoFeature::Transcript
                },
                feature_id: String::from("NM_007294.4"),
                feature_biotype: FeatureBiotype::Coding,
                rank: Some(Rank { ord: 23, total: 23 }),
                hgvs_t: Some(String::from("c.5586C>G")),
                hgvs_p: Some(String::from("p.H1862Q")),
                tx_pos: Some(Pos {
                    ord: 5699,
                    total: Some(7088)
                }),
                cds_pos: Some(Pos {
                    ord: 5586,
                    total: Some(5592)
                }),
                protein_pos: Some(Pos {
                    ord: 1862,
                    total: Some(1864)
                }),
                distance: Some(-1390),
                messages: None,
            }
        );

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
        annotate_opa1_vars("tests/data/annotate/vars/opa1.hand_picked.tsv")
    }

    // Compare to SnpEff annotated ClinVar variants for OPA1 (slow).
    #[ignore]
    #[test]
    fn annotate_opa1_clinvar_vars_snpeff() -> Result<(), anyhow::Error> {
        annotate_opa1_vars("tests/data/annotate/vars/clinvar.excerpt.snpeff.opa1.tsv")
    }

    // Compare to SnpEff annotated ClinVar variants for OPA1 (slow).
    #[ignore]
    #[test]
    fn annotate_opa1_clinvar_vars_vep() -> Result<(), anyhow::Error> {
        annotate_opa1_vars("tests/data/annotate/vars/clinvar.excerpt.vep.opa1.tsv")
    }

    fn annotate_opa1_vars(path_tsv: &str) -> Result<(), anyhow::Error> {
        let txs = vec![
            String::from("NM_001354663.2"),
            String::from("NM_001354664.2"),
            String::from("NM_015560.3"),
            String::from("NM_130831.3"),
            String::from("NM_130832.3"),
            String::from("NM_130837.3"),
        ];

        annotate_vars(path_tsv, &txs)
    }

    // Compare to SnpEff annotated variants for BRCA1, touching special cases.
    #[test]
    fn annotate_brca1_hand_picked_vars() -> Result<(), anyhow::Error> {
        annotate_brca1_vars("tests/data/annotate/vars/brca1.hand_picked.tsv")
    }

    // Compare to SnpEff annotated ClinVar variants for BRCA1 (slow).
    #[ignore]
    #[test]
    fn annotate_brca1_clinvar_vars_snpeff() -> Result<(), anyhow::Error> {
        annotate_brca1_vars("tests/data/annotate/vars/clinvar.excerpt.snpeff.brca1.tsv")
    }

    // Compare to SnpEff annotated ClinVar variants for BRCA1 (slow).
    #[ignore]
    #[test]
    fn annotate_brca1_clinvar_vars_vep() -> Result<(), anyhow::Error> {
        annotate_brca1_vars("tests/data/annotate/vars/clinvar.excerpt.vep.brca1.tsv")
    }

    fn annotate_brca1_vars(path_tsv: &str) -> Result<(), anyhow::Error> {
        let txs = vec![
            String::from("NM_007294.4"),
            String::from("NM_007297.4"),
            String::from("NM_007298.3"),
            String::from("NM_007299.4"),
            String::from("NM_007300.4"),
        ];

        annotate_vars(path_tsv, &txs)
    }

    fn annotate_vars(path_tsv: &str, txs: &[String]) -> Result<(), anyhow::Error> {
        let tx_path = "tests/data/annotate/db/grch37/txs.bin.zst";
        let tx_db = load_tx_db(tx_path)?;
        let provider = Rc::new(MehariProvider::new(tx_db, Assembly::Grch37p10));
        let predictor = ConsequencePredictor::new(provider, Assembly::Grch37p10);

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
                        || (expected_one_of.contains(&String::from("5_prime_UTR_variant"))
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
                        || expected_one_of.contains(&String::from("disruptive_inframe_deletion"))
                            && (record_csqs.contains(&String::from("protein_altering_variant")));
                    // In the case of `GRCh37:17:41258543:T:TA`, the `hgvs` prediction is `c.-1_1insT` and
                    // `p.Met1?` which leads to `start_lost` while VEP predicts `5_prime_UTR_variant`.
                    // This may be a bug in `hgvs` and we don't change this for now.  We accept the call
                    // by VEP, of course.
                    let found_one = found_one
                        || expected_one_of.contains(&String::from("start_lost"))
                            && (record_csqs.contains(&String::from("5_prime_UTR_variant")));
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
