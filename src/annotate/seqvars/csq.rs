//! Compute molecular consequence of variants.

use std::{collections::HashMap, rc::Rc};

use hgvs::{
    data::interface::{Provider, TxForRegionRecord},
    mapper::assembly::{Config as AssemblyConfig, Mapper as AssemblyMapper},
    parser::{Accession, GenomeInterval, GenomeLocEdit, HgvsVariant, Mu, NaEdit},
    static_data::Assembly,
};

use crate::db::create::txs::data::TranscriptBiotype;

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
        // Obtain accession from chromosome name.
        let chrom_acc = self.chrom_to_acc.get(&var.chromosome);
        let chrom_acc = if let Some(chrom_acc) = chrom_acc {
            chrom_acc
        } else {
            tracing::debug!(
                "Could not determine chromosome accession for {:?}; giving up on annotation",
                &var
            );
            return Ok(None);
        };

        // Get all affected transcripts.
        let start = var.position - PADDING;
        let end = var.position + var.reference.len() as i32 - 1 + PADDING;
        let mut txs = self
            .provider
            .get_tx_for_region(&chrom_acc, ALT_ALN_METHOD, start, end)?;
        txs.sort_by(|a, b| a.tx_ac.cmp(&b.tx_ac));

        // Generate `AnnField` records for each transcript.
        Ok(Some(
            txs.into_iter()
                .map(|tx| self.build_ann_field(var, tx, chrom_acc.clone(), start, end))
                .collect::<Result<Vec<_>, _>>()?,
        ))
    }

    fn build_ann_field(
        &self,
        var: &VcfVariant,
        tx_record: TxForRegionRecord,
        chrom_acc: String,
        var_start: i32,
        var_end: i32,
    ) -> Result<AnnField, anyhow::Error> {
        let tx = self.provider.get_tx(&tx_record.tx_ac).unwrap();
        let mut consequences: Vec<Consequence> = Vec::new();

        let alignment = tx.genome_alignments.first().unwrap();

        let mut min_start = std::i32::MAX;
        let mut max_end = std::i32::MIN;

        // Find first exon that overlaps with variant or intron that contains the variant.
        //
        // Note that exons are stored in genome position order.
        let mut prev_end = std::i32::MAX;
        let mut rank = Rank::default();
        let mut is_exonic = false;
        let mut is_intronic = false;
        let mut distance = None;
        let mut tx_len = 0;
        let mut cds_len = 0;
        for exon_alignment in &alignment.exons {
            tx_len += exon_alignment.alt_end_i - exon_alignment.alt_start_i + 1;
            cds_len += exon_alignment
                .alt_cds_end_i
                .map(|alt_cds_end_i| alt_cds_end_i - exon_alignment.alt_cds_start_i.unwrap() + 1)
                .unwrap_or(0);

            let exon_start = exon_alignment.alt_start_i;
            let exon_end = exon_alignment.alt_end_i;
            let intron_start = prev_end;
            let intron_end = exon_start;

            if var_start <= exon_end && exon_start <= var_end {
                // overlaps with exon
                rank = Rank {
                    ord: exon_alignment.ord + 1,
                    total: alignment.exons.len() as i32,
                };
                is_exonic = true;

                if var_start <= exon_start || var_end >= exon_end {
                    // overlaps with exon/intron boundary
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
            } else if var_start >= intron_start && var_end <= intron_end {
                // contained within intron: cannot be in next exon
                rank = Rank {
                    ord: exon_alignment.ord + 1,
                    total: alignment.exons.len() as i32 - 1,
                };
                if !is_exonic {
                    is_intronic = true;

                    let dist_start = -(var_start - intron_start + 1);
                    let dist_end = intron_end - var_end + 1;
                    if dist_end >= dist_start.abs() {
                        distance = Some(dist_end);
                    } else {
                        distance = Some(dist_start);
                    }
                }
            }

            min_start = std::cmp::min(min_start, exon_start);
            max_end = std::cmp::min(max_end, exon_end);
            prev_end = exon_end;
        }
        let prot_len = cds_len / 3;

        let is_upstream = var_end < min_start;
        let is_downstream = var_start > max_end;
        if is_exonic {
            consequences.push(Consequence::ExonVariant);
        } else if is_intronic {
            consequences.push(Consequence::IntronVariant);
        } else if is_upstream {
            let val = -(min_start - var_end + 1);
            if val.abs() < 5_000 {
                consequences.push(Consequence::UpstreamGeneVariant);
            }
            distance = Some(val);
        } else if is_downstream {
            let val = var_start - max_end + 1;
            if val.abs() < 5_000 {
                consequences.push(Consequence::DownstreamGeneVariant);
            }
            distance = Some(val);
        }

        let feature_biotype = match tx.biotype {
            TranscriptBiotype::Coding => FeatureBiotype::Coding,
            TranscriptBiotype::NonCoding => FeatureBiotype::Noncoding,
        };

        let var_g = HgvsVariant::GenomeVariant {
            accession: Accession { value: chrom_acc },
            gene_symbol: None,
            loc_edit: GenomeLocEdit {
                loc: Mu::Certain(GenomeInterval {
                    start: Some(var.position),
                    end: Some(var.position + var.reference.len() as i32 - 1),
                }),
                edit: Mu::Certain(NaEdit::RefAlt {
                    reference: var.reference.clone(),
                    alternative: var.alternative.clone(),
                }),
            },
        };

        let var_n = self.mapper.g_to_n(&var_g, &tx.id)?;
        let tx_pos = match &var_n {
            HgvsVariant::TxVariant { loc_edit, .. } => Some(Pos {
                ord: loc_edit.loc.inner().start.base,
                total: Some(tx_len),
            }),
            _ => panic!("Invalid tx position: {:?}", &var_n),
        };

        let (var_t, _var_p, hgvs_p, cds_pos, protein_pos) = match feature_biotype {
            FeatureBiotype::Coding => {
                let var_c = self.mapper.n_to_c(&var_n)?;
                let var_p = self.mapper.c_to_p(&var_c)?;
                let hgvs_p = Some(format!("{}", &var_p));
                let cds_pos = Some(Pos {
                    ord: -1,
                    total: Some(cds_len),
                }); // TODO
                let protein_pos = Some(Pos {
                    ord: -1,
                    total: Some(prot_len),
                }); // TODO
                (var_c, Some(var_p), hgvs_p, cds_pos, protein_pos)
            }
            FeatureBiotype::Noncoding => (var_n, None, None, None, None),
        };
        let hgvs_t = Some(format!("{}", &var_t));

        // Take a highest-ranking consequence and derive putative impact from it.
        consequences.sort();
        let putative_impact = (*consequences.first().unwrap()).into();

        // Build and return ANN field from the information derived above.
        Ok(AnnField {
            allele: Allele::Alt {
                alternative: var.alternative.clone(),
            },
            consequences,
            putative_impact,
            gene_id: tx.gene_id.clone(),
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
        })
    }
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use crate::annotate::seqvars::{ann::PutativeImpact, load_tx_db};

    use super::*;

    #[test]
    fn annotate_snvs_brca1() -> Result<(), anyhow::Error> {
        let tx_path = "tests/data/annotate/db/seqvars/grch37/txs.bin";
        let tx_db = load_tx_db(tx_path, 5_000_000)?;
        let provider = Rc::new(MehariProvider::new(tx_db));

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
        println!("{:#?}", &res);
        assert_eq!(
            res[0],
            AnnField {
                allele: Allele::Alt {
                    alternative: String::from("C")
                },
                consequences: vec![Consequence::ExonVariant],
                putative_impact: PutativeImpact::Modifier,
                gene_id: String::from("1100"),
                feature_type: FeatureType::SoTerm {
                    term: SoFeature::Transcript
                },
                feature_id: String::from("NM_007294.4"),
                feature_biotype: FeatureBiotype::Coding,
                rank: Rank { ord: 23, total: 23 },
                hgvs_t: Some(String::from("NM_007294.4:c.5586C>G")),
                hgvs_p: Some(String::from("NP_009225.1:p.His1862Gln")),
                tx_pos: Some(Pos {
                    ord: 0,
                    total: Some(7088)
                }),
                cds_pos: Some(Pos {
                    ord: 0,
                    total: Some(5592)
                }),
                protein_pos: Some(Pos {
                    ord: 1862,
                    total: Some(1864)
                }),
                distance: Some(0),
                messages: None,
            }
        );
        // assert_eq!(
        //     format!("{}", res[1]),
        //     "T|exon_variant|MODIFIER|1100|transcript|ENST00000309486.4|Coding|21/22|\
        //     ENST00000309486.4:c.*6229G>A|MD5_555d63f4613b8654d04949014cc53e9a:p.?|11961|||0|"
        // );
        // assert_eq!(
        //     format!("{}", res[2]),
        //     "T|exon_variant|MODIFIER|1100|transcript|ENST00000346315.3|Coding|18/19|\
        //     ENST00000346315.3:c.*6226G>A|MD5_e000d3f294909461e2bdb311a8ceb781:p.?|11295|||0|"
        // );
        // assert_eq!(
        //     format!("{}", res[3]),
        //     "T|exon_variant|MODIFIER|1100|transcript|ENST00000351666.3|Coding|18/19|\
        //     ENST00000351666.3:c.*6226G>A|MD5_7f256726f0eff1d7208ada2c937df1a9:p.?|8288|||0|"
        // );
        // assert_eq!(
        //     format!("{}", res[4]),
        //     "T|exon_variant|MODIFIER|1100|transcript|ENST00000352993.3|Coding|21/22|\
        //     ENST00000352993.3:c.*6229G>A|MD5_920c89c4cdd507c06cc17b6614230efe:p.?|8627|||0|"
        // );
        // assert_eq!(
        //     format!("{}", res[5]),
        //     "T|exon_variant|MODIFIER|1100|transcript|ENST00000354071.3|Coding|17/18|\
        //     ENST00000354071.3:c.*6225G>A|MD5_316c7cafc42206aa3030d4fba92a5297:p.?|11254|||0|"
        // );

        Ok(())
    }
}
