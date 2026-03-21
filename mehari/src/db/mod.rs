//! Database construction and introspection tools.

use crate::pbs;
use crate::pbs::txs::TxSeqDatabase;
use biocommons_bioutils::assemblies::Assembly;

pub mod check;
pub mod create;
pub mod dump;
pub mod merge;
pub mod subset;

/// Trait for transcript databases.
pub trait TranscriptDatabase {
    /// Get the assembly of the transcript database.
    fn assembly(&self) -> Assembly;
}

impl TranscriptDatabase for TxSeqDatabase {
    fn assembly(&self) -> Assembly {
        let assembly = self
            .source_version
            .iter()
            .map(|v| pbs::txs::Assembly::try_from(v.assembly).unwrap())
            .next()
            .expect("At least one source_version entry expected");

        match assembly {
            pbs::txs::Assembly::Grch37 => Assembly::Grch37p10, // has MT
            pbs::txs::Assembly::Grch38 => Assembly::Grch38,
            _ => panic!("Unsupported assembly"),
        }
    }
}
