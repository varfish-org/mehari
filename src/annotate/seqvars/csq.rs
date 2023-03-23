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
        let var_start = var.position;
        let var_end = var.position + var.reference.len() as i32 - 1;
        let qry_start = var_start - PADDING;
        let qry_end = var_end + PADDING;
        let mut txs =
            self.provider
                .get_tx_for_region(chrom_acc, ALT_ALN_METHOD, qry_start, qry_end)?;
        txs.sort_by(|a, b| a.tx_ac.cmp(&b.tx_ac));

        // Generate `AnnField` records for each transcript.
        Ok(Some(
            txs.into_iter()
                .map(|tx| self.build_ann_field(var, tx, chrom_acc.clone(), var_start, var_end))
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
        for exon_alignment in &alignment.exons {
            tx_len += exon_alignment.alt_end_i - exon_alignment.alt_start_i;

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
                if !is_exonic {
                    rank = Rank {
                        ord: exon_alignment.ord + 1,
                        total: alignment.exons.len() as i32 - 1,
                    };
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
        println!(
            "var_n = {}, tx_pos = {:#?}",
            &var_n,
            &tx_pos.as_ref().unwrap()
        );

        let (var_t, _var_p, hgvs_p, cds_pos, protein_pos) = match feature_biotype {
            FeatureBiotype::Coding => {
                let cds_len = tx.stop_codon.unwrap() - tx.start_codon.unwrap();
                let prot_len = cds_len / 3;

                let var_c = self.mapper.n_to_c(&var_n)?;
                let var_p = self.mapper.c_to_p(&var_c)?;
                let hgvs_p = Some(format!("{}", &var_p));
                let cds_pos = match &var_c {
                    HgvsVariant::CdsVariant { loc_edit, .. } => Some(Pos {
                        ord: loc_edit.loc.inner().start.base,
                        total: Some(cds_len),
                    }),
                    _ => panic!("Invalid CDS position: {:?}", &var_n),
                };
                println!(
                    "var_c = {}, cds_pos = {:#?}",
                    &var_c,
                    &cds_pos.as_ref().unwrap()
                );
                let protein_pos = match &var_p {
                    HgvsVariant::ProtVariant { loc_edit, .. } => match &loc_edit {
                        hgvs::parser::ProtLocEdit::Ordinary { loc, .. } => Some(Pos {
                            ord: loc.inner().start.number,
                            total: Some(prot_len),
                        }),
                        _ => None,
                    },
                    _ => panic!("Not a protein position: {:?}", &var_n),
                };
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
                distance: Some(-1391),
                messages: None,
            }
        );

        Ok(())
    }
}