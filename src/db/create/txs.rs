//! Transcript database.

use clap::Parser;

/// Command line arguments for `db create txs` sub command.
#[derive(Parser, Debug)]
#[command(about = "Construct mehari transcripts database", long_about = None)]
pub struct Args {}

/// Main entry point for `db create txs` sub command.
pub fn run(_common: &crate::common::Args, _args: &Args) -> Result<(), anyhow::Error> {
    todo!()
}
