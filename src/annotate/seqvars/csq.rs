//! Compute molecular consequence of variants.

use std::{collections::HashMap, rc::Rc};

use hgvs::{
    data::interface::{Provider, TxForRegionRecord},
    mapper::assembly::{Config as AssemblyConfig, Mapper as AssemblyMapper},
    static_data::Assembly,
};

use crate::world_flatbuffers::mehari::TranscriptBiotype;

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
            ..Default::default()
        };
        let mapper = AssemblyMapper::new(config, provider);

        ConsequencePredictor {
            mapper,
            chrom_to_acc,
        }
    }

    pub fn predict(&self, var: &VcfVariant) -> Result<Option<Vec<AnnField>>, anyhow::Error> {
        // Obtain accession from chromosome name.
        let acc = self.chrom_to_acc.get(&var.chromosome);
        let acc = if let Some(acc) = acc {
            acc
        } else {
            tracing::debug!(
                "Could not determine chromosome accession for {:?}; giving up on annotation",
                &var
            );
            return Ok(None);
        };

        // Get all affected transcripts.
        let start = var.position - PADDING;
        let end = var.position + var.reference.len() as i32 - 1;
        todo!();
        // let txs = self
        //     .provider
        //     .get_tx_for_region(acc.as_ref(), ALT_ALN_METHOD, start, end)?;

        // // Generate `AnnField` records for each transcript.
        // Ok(Some(
        //     txs.into_iter()
        //         .map(|tx| self.build_ann_field(var, tx, start, end))
        //         .collect::<Result<Vec<_>, _>>()?,
        // ))
    }

    fn build_ann_field(
        &self,
        var: &VcfVariant,
        tx_record: TxForRegionRecord,
        var_start: i32,
        var_end: i32,
    ) -> Result<AnnField, anyhow::Error> {
        todo!()
        // let tx = self.provider.get_tx(&tx_record.tx_ac).unwrap();

        // let alignment = tx.genome_alignments().unwrap().get(0);

        // let mut min_start = std::i32::MAX;
        // let mut max_end = std::i32::MIN;

        // // Find first exon that overlaps with variant or intron that contains the variant.
        // //
        // // Note that exons are stored in genome position order.
        // let mut prev_end = std::i32::MAX;
        // let mut rank = Rank::default();
        // for exon_alignment in alignment.exons().unwrap() {
        //     let exon_start = exon_alignment.alt_start_i();
        //     let exon_end = exon_alignment.alt_end_i();
        //     let intron_start = prev_end;
        //     let intron_end = exon_start;

        //     if var_start <= exon_end && exon_start <= var_end {
        //         // overlaps with exon
        //         rank = Rank {
        //             ord: exon_alignment.ord(),
        //             total: alignment.exons().unwrap().len() as i32,
        //         };
        //         break;
        //     } else if var_start >= intron_start && var_end <= intron_end {
        //         // contained within intron: cannot be in next exon
        //         rank = Rank {
        //             ord: exon_alignment.ord(),
        //             total: alignment.exons().unwrap().len() as i32 - 1,
        //         };
        //         break;
        //     }

        //     min_start = std::cmp::min(min_start, exon_start);
        //     max_end = std::cmp::min(max_end, exon_end);
        //     prev_end = exon_end;
        // }

        // // TODO: fill the BELOW!!!
        // let hgvs_c = Some(String::from("HGVS.c"));
        // let hgvs_p = Some(String::from("HGVS.p"));
        // let mut consequences: Vec<Consequence> = Vec::new();
        // let cdna_pos = Some(Pos {
        //     ord: 1,
        //     total: Some(2),
        // });
        // let cds_pos = Some(Pos {
        //     ord: 1,
        //     total: Some(2),
        // });
        // let protein_pos = Some(Pos {
        //     ord: 1,
        //     total: Some(2),
        // });
        // let distance = Some(1);
        // // TODO: fill the ABOVE!!!

        // // Take a highest-ranking consequence and derive putative impact from it.
        // consequences.sort();
        // let putative_impact = (*consequences.first().unwrap()).into();

        // // Build and return ANN field from the information derived above.
        // Ok(AnnField {
        //     allele: Allele::Alt {
        //         alternative: var.alternative.clone(),
        //     },
        //     consequences,
        //     putative_impact,
        //     gene_id: tx.gene_id().unwrap().to_string(),
        //     feature_type: FeatureType::SoTerm {
        //         term: SoFeature::Transcript,
        //     },
        //     feature_id: tx.id().unwrap().to_string(),
        //     feature_biotype: match tx.biotype() {
        //         TranscriptBiotype::Coding => FeatureBiotype::Coding,
        //         TranscriptBiotype::NonCoding => FeatureBiotype::Noncoding,
        //         _ => panic!("Unexpected flatbuffers biotype: {:?}", tx.biotype()),
        //     },
        //     rank,
        //     hgvs_c,
        //     hgvs_p,
        //     cdna_pos,
        //     cds_pos,
        //     protein_pos,
        //     distance,
        //     messages: None,
        // })
    }
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use crate::annotate::seqvars::load_tx_db;

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
                position: 43_045_682,
                reference: String::from("T"),
                alternative: String::from("C"),
            })?
            .unwrap();

        assert_eq!(format!("{:?}", res), "");

        Ok(())
    }
}
