//! Mehari is a software package for annotating VCF files with variant effect/consequence.
//!
//! The program uses hgvs-rs for projecting genomic variants to transcripts and proteins and
//! thus has high prediction quality.  The aim is to provide variant effect prediction that
//! is identical with the VariantValidator (and the used `hgvs` Python software package).
//!
//! **If you are not a Rust programmer than you most likely want to install Mehari using
//! (bio)conda, use the command line interface, and use pre-built databases.**
//!
//! Below, you can find some overview.  For more details, please refer to the
//!
//! - [end user documentation](`self::user_doc`).
//!
//! ## Installation
//!
//! Mehari provides a command line interface (CLI).  You can install the software either using
//! `cargo`:
//!
//! ```text
//! $ cargo install mehari
//! ```
//!
//! Or if you prefer, use the conda package manager from the bioconda channel.
//!
//! ```text
//! $ conda install -c bioconda mehari
//! ```
//!
//! ## Command Line Usage
//!
//! To run Mehari, invoke the `mehari` executable.  The program provides a command line help:
//!
//! ```text
//! $ mehari help
//! VCF variant effect prediction and annotation
//!
//! Usage: mehari [OPTIONS] <COMMAND>
//!
//! Commands:
//!   db        Database-related commands
//!   annotate  Annotation related commands
//!   help      Print this message or the help of the given subcommand(s)
//!
//! Options:
//!   -v, --verbose...  More output per occurrence
//!   -q, --quiet...    Less output per occurrence
//!   -h, --help        Print help
//!   -V, --version     Print version
//! ```
//!
//! The most important commands are ...
//!
//! ... for annotating VCF files:
//!
//! * `mehari annotate seqvars`
//!
//! ... for building the databases:
//!
//! * `mehari db create txs` -- create a database of transcript sequences from
//!   [cdot JSON files](https://github.com/SACGF/cdot/releases)
//! * `mehari db create seqvar-freqs` -- create a database of sequence variant (SNV/indel/MNV)
//!   population frequencies (from [gnomAD](https://gnomad.broadinstitute.org/) and
//!   [HelixMtDb](https://www.helix.com/pages/mitochondrial-variant-database))
//!
//! Full documentation is available in the [user documentation](`self::user_doc`).
//!
//! ## Library Usage
//!
//! You can also use Mehari as a library.  Simply add the `mehari` crate to your Cargo.toml file:
//!
//! ```text
//! $ cargo add mehari
//! ```

use clap::{command, Args, Parser, Subcommand};

use mehari::{annotate, common, db, server, verify};

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
    /// Verification related commands.
    Verify(Verify),
    /// Server related commands.
    RunServer(server::Args),
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
    Strucvars(annotate::strucvars::Args),
}

/// Parsing of "verify *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Verify {
    /// The sub command to run
    #[command(subcommand)]
    command: VerifyCommands,
}

/// Enum supporting the parsing of "annotate *" sub commands.
#[derive(Debug, Subcommand)]
enum VerifyCommands {
    Seqvars(verify::seqvars::Args),
}

#[tokio::main]
async fn main() -> Result<(), anyhow::Error> {
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
    tracing::subscriber::set_global_default(collector)?;

    // Install collector and go into sub commands.
    tracing::info!("Mehari startup -- letting the dromedary off the leash...");

    match &cli.command {
        Commands::Db(db) => match &db.command {
            DbCommands::Create(db_create) => match &db_create.command {
                DbCreateCommands::Txs(args) => db::create::txs::run(&cli.common, args)?,
            },
        },
        Commands::Annotate(annotate) => match &annotate.command {
            AnnotateCommands::Seqvars(args) => annotate::seqvars::run(&cli.common, args)?,
            AnnotateCommands::Strucvars(args) => {
                annotate::strucvars::run(&cli.common, args).await?
            }
        },
        Commands::Verify(verify) => match &verify.command {
            VerifyCommands::Seqvars(args) => verify::seqvars::run(&cli.common, args)?,
        },
        Commands::RunServer(args) => server::run(&cli.common, args).await?,
    }

    tracing::info!("... the dromedary is back in the stable.");
    tracing::info!("All done. Have a nice day!");

    Ok(())
}
