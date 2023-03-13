//! Transcript database.

use clap::Parser;

/// Command line arguments for `db create txs` sub command.
#[derive(Parser, Debug)]
#[command(about = "Construct mehari transcripts and sequence database", long_about = None)]
pub struct Args {
    /// Path to output flatbuffers file to write to.
    #[arg(long)]
    pub path_out: String,
    /// Paths to the cdot JSON transcripts to import.
    #[arg(long, required = true)]
    pub path_cdot_json: Vec<String>,
    /// Path to the seqrepo instance directory to use.
    #[arg(long)]
    pub path_seqrepo_instance: String,
}

/// Main entry point for `db create txs` sub command.
pub fn run(_common: &crate::common::Args, _args: &Args) -> Result<(), anyhow::Error> {
    todo!()
}
