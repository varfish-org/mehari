//! Database construction and introspection tools.

use crate::pbs::txs::TxSeqDatabase;

pub mod check;
pub mod create;
pub mod dump;
pub mod merge;
pub mod subset;

/// Trait for transcript databases.
pub trait TranscriptDatabase {
    /// Get the assembly of the transcript database.
    fn assembly(&self) -> String;
}

impl TranscriptDatabase for TxSeqDatabase {
    fn assembly(&self) -> String {
        let assembly = self
            .source_version
            .iter()
            .map(|v| v.assembly.as_str())
            .next()
            .expect("At least one source_version entry expected");

        match assembly {
            "grch37" | "grch37p10" => "grch37".into(),
            "grch38" => "grch38".into(),
            x => x.into(),
        }
    }
}
