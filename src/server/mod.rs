use crate::{
    annotate::{
        seqvars::{
            csq::ConsequencePredictor as SeqvarConsequencePredictor, load_tx_db, path_component,
            provider::Provider as MehariProvider,
        },
        strucvars::csq::ConsequencePredictor as StrucvarConsequencePredictor,
    },
    common::GenomeRelease,
};
use clap::ValueEnum;
use std::path::PathBuf;
use std::sync::Arc;

/// Implementation of Actix server.
pub mod actix_server;

/// Command line arguments for "run-server` command.
#[derive(clap::Parser, Debug)]
#[command(about = "Run Mehari REST API server", long_about = None)]
pub struct Args {
    /// Paths to genome fasta files.
    #[arg(long)]
    pub reference: Vec<PathBuf>,

    /// Paths to transcript databases.
    #[arg(long)]
    pub transcripts: Vec<PathBuf>,

    /// Whether to suppress printing hints.
    #[arg(long, default_value_t = false)]
    pub suppress_hints: bool,

    /// IP to listen on.
    #[arg(long, default_value = "127.0.0.1")]
    pub listen_host: String,
    /// Port to listen on.
    #[arg(long, default_value_t = 8080)]
    pub listen_port: u16,
}

/// Print some hints via `tracing::info!`.
pub fn print_hints(args: &Args) {
    tracing::info!(
        "Launching server main on http://{}:{} ...",
        args.listen_host.as_str(),
        args.listen_port
    );

    // Short-circuit if no hints are to be
    if args.suppress_hints {
        return;
    }

    // The endpoint `/genes/txs` provides transcript information.
    tracing::info!(
        "  try: http://{}:{}/genes/txs?hgncId=HGNC:1100&\
        genomeBuild=GENOME_BUILD_GRCH37",
        args.listen_host.as_str(),
        args.listen_port
    );
    // The endpoint `/tx/csq` to comput ethe consequence of a variant; without and with filtering
    // for HGNC gene ID.
    tracing::info!(
        "  try: http://{}:{}/seqvars/csq?genome_release=grch37\
        &chromosome=17&position=48275363&reference=C&alternative=A",
        args.listen_host.as_str(),
        args.listen_port
    );
    tracing::info!(
        "  try: http://{}:{}/seqvars/csq?genome_release=grch37\
        &chromosome=17&position=48275363&reference=C&alternative=A&hgnc_id=HGNC:2197",
        args.listen_host.as_str(),
        args.listen_port
    );
    // The endpoint `/strucvars/csq` computes the consequence of an SV.
    tracing::info!(
        "  try: http://{}:{}/strucvars/csq?genome_release=grch37\
        &chromosome=17&start=48275360&&stop=48275370&sv_type=DEL",
        args.listen_host.as_str(),
        args.listen_port
    );
    // The endpoint `/structvars/csq` computes the consequence of an SV.
    tracing::info!(
        "  try: http://{}:{}/strucvars/csq?genome_release=grch37\
        &chromosome=17&start=48275360&&stop=48275370&sv_type=DEL",
        args.listen_host.as_str(),
        args.listen_port
    );
}

/// Main entry point for `run-server` sub command.
///
/// # Errors
///
/// In the case that there is an error running the server.
pub async fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("args_common = {:?}", &args_common);
    tracing::info!("args = {:?}", &args);

    if let Some(level) = args_common.verbose.log_level() {
        match level {
            log::Level::Trace | log::Level::Debug => {
                std::env::set_var("RUST_LOG", "debug");
                env_logger::init_from_env(env_logger::Env::new().default_filter_or("info"));
            }
            _ => (),
        }
    }

    // Load data that we need for running the server.
    tracing::info!("Loading database...");
    let before_loading = std::time::Instant::now();
    let mut data = actix_server::WebServerData::default();
    for (reference_path, txdb_path) in args.reference.iter().zip(args.transcripts.iter()) {
        let reference = crate::annotate::seqvars::reference::genome_reference(reference_path)?;
        let assembly = reference.guess_assembly().unwrap_or_else(|| {
            panic!("Could not guess assembly for {}", &reference_path.display());
        });
        let genome_release = GenomeRelease::from(assembly);

        if !reference_path.exists() {
            tracing::warn!("No reference fasta found at {}", reference_path.display());
            continue;
        }
        if !txdb_path.exists() {
            tracing::warn!("No transcript database found at {}", txdb_path.display());
            continue;
        }

        tracing::info!("  - loading {}", txdb_path.display());
        let tx_db = load_tx_db(&txdb_path)?;
        assert_eq!(
            GenomeRelease::from_str(&tx_db.genome_release.as_ref().unwrap(), true).unwrap(),
            assembly.into(),
            "Genome release mismatch between reference fasta and txdb"
        );

        tracing::info!("  - building interval trees");

        let provider = Arc::new(MehariProvider::new(
            tx_db,
            reference.clone(),
            Default::default(),
        ));
        data.provider.insert(genome_release, provider.clone());

        tracing::info!("  - building seqvars predictors");
        data.seqvars_predictors.insert(
            genome_release,
            SeqvarConsequencePredictor::new(provider.clone(), Default::default()),
        );

        tracing::info!("  - building strucvars predictors");
        data.strucvars_predictors.insert(
            genome_release,
            StrucvarConsequencePredictor::new(provider.clone()),
        );
    }
    let data = actix_web::web::Data::new(data);
    tracing::info!("...done loading data {:?}", before_loading.elapsed());

    // Print the server URL and some hints (the latter: unless suppressed).
    print_hints(args);
    // Launch the Actix web server.
    actix_server::main(args, data).await?;

    tracing::info!("All done. Have a nice day!");
    Ok(())
}
