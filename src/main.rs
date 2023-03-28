//! Mehari is a software package for annotating VCF files with variant effect/consequence.
//!
//! The program uses hgvs-rs for projecting genomic variants to transcripts and proteins and
//! thus has high prediction quality.  The aim is to provide variant effect prediction that
//! is identical with the VariantValidator (and the used `hgvs` Python software package).
//!
//! ## Command Line Usage
//!
//! ## Library Usage
//!
//! You can also use Mehari as a library.  Simply add the `mehari` crate to your Cargo.toml file:
//!
//! ```text
//! $ cargo add mehari
//! ```

pub mod annotate;
pub mod common;
pub mod db;

#[allow(
    non_snake_case,
    unused_imports,
    clippy::extra_unused_lifetimes,
    clippy::missing_safety_doc,
    clippy::derivable_impls,
    clippy::size_of_in_element_count
)]
#[path = "../target/flatbuffers/world_generated.rs"]
pub mod world_flatbuffers;

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
    /// Annotation related commands.
    Annotate(Annotate),
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

/// Parsing of "annotate *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Annotate {
    /// The sub command to run
    #[command(subcommand)]
    command: AnnotateCommands,
}

/// Enum supporting the parsing of "annotate *" sub commands.
#[derive(Debug, Subcommand)]
enum AnnotateCommands {
    Seqvars(annotate::seqvars::Args),
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
        tracing::info!("Mehari startup -- letting the dromedary off the leash...");

        match &cli.command {
            Commands::Db(db) => match &db.command {
                DbCommands::Create(db_create) => match &db_create.command {
                    DbCreateCommands::Txs(args) => db::create::txs::run(&cli.common, args)?,
                    DbCreateCommands::SeqvarFreqs(args) => {
                        db::create::seqvar_freqs::run(&cli.common, args)?
                    }
                },
            },
            Commands::Annotate(annotate) => match &annotate.command {
                AnnotateCommands::Seqvars(args) => annotate::seqvars::run(&cli.common, args)?,
            },
        }

        tracing::info!("All done. Have a nice day!");

        Ok::<(), anyhow::Error>(())
    })?;

    Ok(())
}
