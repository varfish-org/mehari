//! Run the server.

use std::sync::Arc;

use hgvs::static_data::Assembly;

use crate::annotate::seqvars::{
    csq::ConsequencePredictor, load_tx_db, path_component, provider::MehariProvider,
};

/// Implementation of Actix server.
pub mod actix_server {
    use hgvs::static_data::Assembly;

    use crate::annotate::seqvars::csq::ConsequencePredictor;

    /// Data structure for the web server data.
    #[derive(Debug, Default)]
    pub struct WebServerData {
        /// The consequence predictors for each assembly.
        pub predictors: std::collections::HashMap<Assembly, ConsequencePredictor>,
    }

    /// Main entry point for running the REST server.
    #[allow(clippy::unused_async)]
    #[actix_web::main]
    pub async fn main(
        args: &super::Args,
        data: actix_web::web::Data<WebServerData>,
    ) -> std::io::Result<()> {
        actix_web::HttpServer::new(move || {
            actix_web::App::new()
                .app_data(data.clone())
                // .service(hpo_genes::handle)
                .wrap(actix_web::middleware::Logger::default())
        })
        .bind((args.listen_host.as_str(), args.listen_port))?
        .run()
        .await
    }
}

/// Command line arguments for "run-server` command.
#[derive(clap::Parser, Debug)]
#[command(about = "Run Mehari REST API server", long_about = None)]
pub struct Args {
    /// Path to the mehari database folder.
    #[arg(long)]
    pub path_db: String,

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
    // // The endpoint `/hpo/sim/term-gene` allows to compute the same for a list of `terms` and
    // // `gene_symbols`.
    // tracing::info!(
    //     "  try: http://{}:{}/hpo/sim/term-gene?terms=HP:0001166,HP:0000098&gene_symbols=FBN1,TGDS,TTN",
    //     args.listen_host.as_str(),
    //     args.listen_port
    // );
}

/// Main entry point for `run-server` sub command.
///
/// # Errors
///
/// In the case that there is an error running the server.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
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
    for assembly in [Assembly::Grch37, Assembly::Grch38] {
        let path = format!("{}/{}/txs.bin.zst", &args.path_db, path_component(assembly));
        tracing::info!("  - loading {}", &path);
        let tx_db = load_tx_db(&path)?;
        tracing::info!("  - building interval trees");
        let provider = Arc::new(MehariProvider::new(tx_db, assembly));
        let predictor = ConsequencePredictor::new(provider, assembly);
        data.predictors.insert(assembly, predictor);
    }
    let data = actix_web::web::Data::new(data);
    tracing::info!("...done loading data {:?}", before_loading.elapsed());

    // Print the server URL and some hints (the latter: unless suppressed).
    print_hints(args);
    // Launch the Actix web server.
    actix_server::main(args, data)?;

    tracing::info!("All done. Have a nice day!");
    Ok(())
}
