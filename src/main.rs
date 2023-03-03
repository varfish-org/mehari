//! Main entry point for the Mehari CLI.

pub mod common;
pub mod db;

use clap::{command, Args, Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "VCF variant effect prediction and annotation"
)]
struct Cli {
    /// Commonly used arguments
    #[command(flatten)]
    common: common::Args,

    /// The sub command to run
    #[command(subcommand)]
    command: Commands,
}

/// Enum supporting the parsing of top-level commands.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Subcommand)]
enum Commands {
    /// Database-related commands.
    Db(Db),
    // /// SV related commands.
    // Sv(Sv),
    // /// Server related commands.
    // Server(Server),
}

/// Parsing of "db *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Db {
    /// The sub command to run
    #[command(subcommand)]
    command: DbCommands,
}

/// Enum supporting the parsing of "db *" sub commands.
#[derive(Debug, Subcommand)]
enum DbCommands {
    Create(DbCreate),
}

/// Parsing of "db create *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct DbCreate {
    /// The sub command to run
    #[command(subcommand)]
    command: DbCreateCommands,
}

/// Enum supporting the parsing of "db create *" sub commands.
#[derive(Debug, Subcommand)]
enum DbCreateCommands {
    Txs(db::create::txs::Args),
    SeqvarFreqs(db::create::seqvar_freqs::Args),
}

fn main() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    // Build a tracing subscriber according to the configuration in `cli.common`.
    let collector = tracing_subscriber::fmt()
        .with_target(false)
        .with_max_level(match cli.common.verbose.log_level() {
            Some(level) => match level {
                log::Level::Error => tracing::Level::ERROR,
                log::Level::Warn => tracing::Level::WARN,
                log::Level::Info => tracing::Level::INFO,
                log::Level::Debug => tracing::Level::DEBUG,
                log::Level::Trace => tracing::Level::TRACE,
            },
            None => tracing::Level::INFO,
        })
        .compact()
        .finish();

    // Install collector and go into sub commands.
    tracing::subscriber::with_default(collector, || {
        tracing::info!("Mehari startup -- letting the camel from the leash...");

        match &cli.command {
            Commands::Db(db) => match &db.command {
                DbCommands::Create(db_create) => match &db_create.command {
                    DbCreateCommands::Txs(args) => db::create::txs::run(&cli.common, args)?,
                    DbCreateCommands::SeqvarFreqs(args) => {
                        db::create::seqvar_freqs::run(&cli.common, args)?
                    }
                },
            },
        }

        tracing::info!("All done. Have a nice day!");

        Ok::<(), anyhow::Error>(())
    })?;

    Ok(())
}
