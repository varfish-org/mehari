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
        let source_version = self
            .source_version
            .iter()
            .next()
            .expect("At least one source_version entry expected");

        // Prefer the new string field, fall back to deprecated enum field
        let assembly = if !source_version.assembly.trim().is_empty() {
            source_version.assembly.as_str()
        } else {
            // Fall back to deprecated enum field
            #[allow(deprecated)]
            match source_version.assembly_enum() {
                crate::pbs::txs::Assembly::Grch37 => "grch37",
                crate::pbs::txs::Assembly::Grch38 => "grch38",
                _ => "",
            }
        };

        match assembly {
            "grch37" | "grch37p10" => "grch37".into(),
            "grch38" => "grch38".into(),
            x => x.into(),
        }
    }
}
