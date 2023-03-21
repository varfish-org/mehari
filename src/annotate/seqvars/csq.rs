//! Compute molecular consequence of variants.

use crate::world_flatbuffers::mehari::TxSeqDatabase;

use super::{
    ann::{Allele, AnnField, FeatureBiotype, FeatureType, Pos, PutativeImpact, Rank, SoFeature},
    TxIntervalTrees,
};

pub struct ConsequencePredictor<'a> {
    tx_db: TxSeqDatabase<'a>,
    tx_trees: TxIntervalTrees,
}

impl<'a> ConsequencePredictor<'a> {
    pub fn from(tx_db: TxSeqDatabase<'a>, tx_trees: TxIntervalTrees) -> Self {
        ConsequencePredictor { tx_db, tx_trees }
    }

    pub fn predict(&self) -> AnnField {
        AnnField {
            allele: Allele::Alt {
                alternative: String::from("A"),
            },
            consequences: Vec::new(),
            putative_impact: PutativeImpact::Moderate,
            gene_id: String::from("xxx"),
            feature_type: FeatureType::SoTerm {
                term: SoFeature::Transcript,
            },
            feature_id: String::from("xxx"),
            feature_biotype: FeatureBiotype::Coding,
            rank: Rank { ord: 1, total: 2 },
            hgvs_c: Some(String::from("HGVS.c")),
            hgvs_p: Some(String::from("HGVS.p")),
            cdna_pos: Some(Pos {
                ord: 1,
                total: Some(2),
            }),
            cds_pos: Some(Pos {
                ord: 1,
                total: Some(2),
            }),
            protein_pos: Some(Pos {
                ord: 1,
                total: Some(2),
            }),
            distance: Some(1),
            messages: None,
        }
    }
}

#[cfg(test)]
mod test {
    use std::fs::File;

    use memmap2::Mmap;
    use pretty_assertions::assert_eq;

    use crate::{annotate::seqvars::TxIntervalTrees, world_flatbuffers::mehari::TxSeqDatabase};

    use super::*;

    #[test]
    fn annotate_snvs_brca1() -> Result<(), anyhow::Error> {
        let tx_path = "tests/data/annotate/db/seqvars/grch37/txs.bin";
        let tx_file = File::open(tx_path)?;
        let tx_mmap = unsafe { Mmap::map(&tx_file)? };
        let tx_db = flatbuffers::root::<TxSeqDatabase>(&tx_mmap)?;
        let tx_trees = TxIntervalTrees::new(&tx_db);

        let predictor = ConsequencePredictor::from(tx_db, tx_trees);

        let res = predictor.predict();

        Ok(())
    }
}
