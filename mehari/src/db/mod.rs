//! Database construction and introspection tools.

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
            .map(|v| v.assembly.as_str())
            .next()
            .expect("At least one source_version entry expected");

        match assembly {
            "grch37" | "grch37p10" => Assembly::Grch37p10, // has MT
            "grch38" => Assembly::Grch38,
            _ => panic!("Unsupported assembly"),
        }
    }
}
