//! Dump transcript database.

use std::{io::Write, path::PathBuf};

use clap::Parser;

use crate::annotate::seqvars::load_tx_db;

/// Command line arguments for `db dump` sub command.
#[derive(Parser, Debug)]
#[command(about = "Dump transcript database", long_about = None)]
pub struct Args {
    /// Path to database file to dump
    #[arg(long)]
    pub path_db: PathBuf,
}

/// Run with a Write.
pub fn run_with_write<W: Write>(
    _common: &crate::common::Args,
    args: &Args,
    writer: &mut W,
) -> Result<(), anyhow::Error> {
    tracing::info!("Opening transcript database");
    let tx_db = load_tx_db(&format!("{}", args.path_db.display()))?;
    tracing::info!("Dumping ...");
    serde_yaml::to_writer(writer, &tx_db)?;
    tracing::info!("... done");

    Ok(())
}

/// Main entry point for `db create txs` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    run_with_write(common, args, &mut std::io::stdout())
}
