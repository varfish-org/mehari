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
//! * `mehari db transcripts create` -- create a database of transcript sequences from
//!   [cdot JSON files](https://github.com/SACGF/cdot/releases)
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

#[cfg(all(
    feature = "jemalloc",
    not(target_env = "msvc"),
    not(feature = "dhat-heap")
))]
use tikv_jemallocator::Jemalloc;

#[cfg(all(
    feature = "jemalloc",
    not(target_env = "msvc"),
    not(feature = "dhat-heap")
))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

#[cfg(all(feature = "dhat-heap", not(feature = "jemalloc")))]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use clap::{Args, Parser, Subcommand};
use clap_verbosity_flag::{InfoLevel, Verbosity};
use mehari::db::transcripts::create;
#[cfg(feature = "server")]
use mehari::server;
use mehari::{annotate, common, db, verify};
use tracing::level_filters::LevelFilter;
use tracing_subscriber::filter::Targets;
use tracing_subscriber::layer::SubscriberExt;

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

    #[cfg(feature = "server")]
    /// Server related commands.
    Server(Server),

    /// Verification related commands.
    Verify(Verify),

    /// Post-processing utilities for Mehari-annotated files.
    Postprocess(Postprocess),
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
    /// Commands related to transcript database
    Transcripts(TranscriptsArgs),
}

/// Subcommands under "db transcripts"
#[derive(Debug, Parser)]
#[command(args_conflicts_with_subcommands = true)]
pub struct TranscriptsArgs {
    /// The sub command to run
    #[command(subcommand)]
    pub command: TranscriptsCommands,
}

#[derive(Debug, Subcommand)]
pub enum TranscriptsCommands {
    /// Construct transcripts database
    Create(create::cli::Args),
    /// Introspect / check transcripts database
    Check(db::transcripts::check::Args),
    /// Dump transcripts database
    Dump(db::transcripts::dump::Args),
    /// Merge transcript databases
    Merge(db::transcripts::merge::Args),
    /// Subset transcripts database
    Subset(db::transcripts::subset::Args),
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
#[allow(clippy::large_enum_variant)]
enum AnnotateCommands {
    Seqvars(annotate::seqvars::Args),
    Strucvars(annotate::strucvars::Args),
}

#[cfg(feature = "server")]
/// Parsing of "server *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Server {
    /// The sub command to run
    #[command(subcommand)]
    command: ServerCommands,
}

#[cfg(feature = "server")]
/// Enum supporting the parsing of "server *" sub commands.
#[derive(Debug, Subcommand)]
#[allow(clippy::large_enum_variant)]
enum ServerCommands {
    Run(server::run::Args),
    Schema(server::schema::Args),
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

#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Postprocess {
    #[command(subcommand)]
    command: PostprocessCommands,
}

/// Enum supporting the parsing of "postprocess *" sub commands.
#[derive(Debug, Subcommand)]
enum PostprocessCommands {
    /// Experimental: Extract the proteome from a mehari-annotated (seqvars, with aa sequences) VCF file.
    Proteome(mehari::postprocess::proteome::Args),

    /// Convert a mehari-annotated (seqvars) VCF file to a VarFish-compatible TSV file.
    VarfishSeqvars(mehari::postprocess::varfish::Args),
}

#[tokio::main]
async fn main() -> Result<(), anyhow::Error> {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let cli = Cli::parse();

    // Build a tracing subscriber according to the configuration in `cli.common`.
    let level = max_level(&cli.common.verbose);
    // hgvs logs details of single variants at INFO, as biocommons/hgvs does.
    let hgvs_level = if level == LevelFilter::INFO {
        LevelFilter::WARN
    } else {
        level
    };
    // A server runs for long, so it logs every warning.
    #[cfg(feature = "server")]
    let is_server = matches!(cli.command, Commands::Server(_));
    #[cfg(not(feature = "server"))]
    let is_server = false;
    let log_limit = (!is_server).then(common::log_limit::LogLimit::default);
    let collector = tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .with_target(false)
        .with_max_level(level)
        .compact()
        .finish()
        .with(
            Targets::new()
                .with_default(level)
                .with_target("hgvs", hgvs_level),
        )
        .with(log_limit.clone());
    tracing::subscriber::set_global_default(collector)?;

    // Install collector and go into sub commands.
    tracing::info!("Mehari startup -- letting the dromedary off the leash...");

    common::log_limit::summarize(log_limit.as_ref(), run(&cli)).await?;
    tracing::info!("... the dromedary is back in the stable.");
    tracing::info!("All done. Have a nice day!");

    Ok(())
}

/// Run the command of `cli`.
async fn run(cli: &Cli) -> Result<(), anyhow::Error> {
    match &cli.command {
        Commands::Db(db) => match &db.command {
            DbCommands::Transcripts(transcripts) => match &transcripts.command {
                TranscriptsCommands::Create(args) => {
                    db::transcripts::create::run(&cli.common, args)?
                }
                TranscriptsCommands::Check(args) => db::transcripts::check::run(&cli.common, args)?,
                TranscriptsCommands::Dump(args) => db::transcripts::dump::run(&cli.common, args)?,
                TranscriptsCommands::Subset(args) => {
                    db::transcripts::subset::run(&cli.common, args)?
                }
                TranscriptsCommands::Merge(args) => db::transcripts::merge::run(&cli.common, args)?,
            },
        },
        Commands::Annotate(annotate) => match &annotate.command {
            AnnotateCommands::Seqvars(args) => annotate::seqvars::run(&cli.common, args).await?,
            AnnotateCommands::Strucvars(args) => {
                annotate::strucvars::run(&cli.common, args).await?
            }
        },
        #[cfg(feature = "server")]
        Commands::Server(server) => match &server.command {
            ServerCommands::Run(args) => server::run::run(&cli.common, args).await?,
            ServerCommands::Schema(args) => server::schema::run(&cli.common, args)?,
        },
        Commands::Verify(verify) => match &verify.command {
            VerifyCommands::Seqvars(args) => verify::seqvars::run(&cli.common, args)?,
        },
        Commands::Postprocess(postprocess) => match &postprocess.command {
            PostprocessCommands::Proteome(args) => {
                mehari::postprocess::proteome::run(&cli.common, args).await?
            }
            PostprocessCommands::VarfishSeqvars(args) => {
                mehari::postprocess::varfish::run(&cli.common, args).await?
            }
        },
    }

    Ok(())
}

/// The most verbose level that `verbose` logs.
fn max_level(verbose: &Verbosity<InfoLevel>) -> LevelFilter {
    match verbose.log_level() {
        Some(level) => match level {
            log::Level::Error => LevelFilter::ERROR,
            log::Level::Warn => LevelFilter::WARN,
            log::Level::Info => LevelFilter::INFO,
            log::Level::Debug => LevelFilter::DEBUG,
            log::Level::Trace => LevelFilter::TRACE,
        },
        None => LevelFilter::OFF,
    }
}

#[cfg(test)]
mod test {
    use clap_verbosity_flag::Verbosity;
    use tracing::level_filters::LevelFilter;

    /// mehari logs from INFO on; each `-v` logs one level more, each `-q` one less, down to
    /// nothing.
    #[rstest::rstest]
    #[case::trace(2, 0, LevelFilter::TRACE)]
    #[case::debug(1, 0, LevelFilter::DEBUG)]
    #[case::info(0, 0, LevelFilter::INFO)]
    #[case::warn(0, 1, LevelFilter::WARN)]
    #[case::error(0, 2, LevelFilter::ERROR)]
    #[case::off(0, 3, LevelFilter::OFF)]
    fn max_level(#[case] verbose: u8, #[case] quiet: u8, #[case] expected: LevelFilter) {
        assert_eq!(super::max_level(&Verbosity::new(verbose, quiet)), expected);
    }
}
