//! Sequence variant frequency annotation.

use clap::Parser;

/// Command line arguments for `db create seqvar-freqs` sub command.
#[derive(Parser, Debug)]
#[command(about = "Construct mehari sequence variant frequencies database", long_about = None)]
pub struct Args {}

/// Main entry point for `db create seqvar_freqs` sub command.
pub(crate) fn run(common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!(
        "Building sequence variant frequencies table\ncommon args: {:#?}\nargs: {:#?}",
        common,
        args
    );

    tracing::info!("Done building sequence variant frequency table");
    Ok(())
}
