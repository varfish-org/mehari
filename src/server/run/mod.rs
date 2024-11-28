use crate::annotate::cli::{Sources, TranscriptSettings};
use crate::annotate::seqvars::{setup_seqvars_annotator, ClinvarAnnotator, FrequencyAnnotator};
use crate::db::merge::merge_transcript_databases;
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
use anyhow::Error;
use std::sync::Arc;

/// Implementation of Actix server.
pub mod actix_server;

/// Module with OpenAPI documentation.
pub mod openapi {
    use crate::annotate::seqvars::ann::{
        Consequence, FeatureBiotype, FeatureType, Message, Pos, PutativeImpact, Rank, SoFeature,
    };
    use crate::annotate::strucvars::csq::interface::StrucvarsSvType;
    use crate::annotate::strucvars::csq::{
        StrucvarsGeneTranscriptEffects, StrucvarsTranscriptEffect,
    };
    use crate::common::GenomeRelease;
    use crate::server::run::actix_server::gene_txs::{
        ExonAlignment, GenesTranscriptsListQuery, GenesTranscriptsListResponse, GenomeAlignment,
        Strand, Transcript, TranscriptBiotype, TranscriptTag,
    };
    use crate::server::run::actix_server::seqvars_csq::{
        SeqvarsCsqQuery, SeqvarsCsqResponse, SeqvarsCsqResultEntry,
    };
    use crate::server::run::actix_server::strucvars_csq::{
        StrucvarsCsqQuery, StrucvarsCsqResponse,
    };
    use crate::server::run::actix_server::versions::{
        Assembly, DataVersionEntry, SoftwareVersions, VersionsInfoResponse,
    };

    use super::actix_server::{gene_txs, seqvars_csq, strucvars_csq, versions, CustomError};

    /// Utoipa-based `OpenAPI` generation helper.
    #[derive(utoipa::OpenApi)]
    #[openapi(
        paths(
            versions::handle,
            gene_txs::handle_with_openapi,
            seqvars_csq::handle_with_openapi,
            strucvars_csq::handle_with_openapi,
            strucvars_csq::handle_with_openapi
        ),
        components(schemas(
            Assembly,
            CustomError,
            VersionsInfoResponse,
            SoftwareVersions,
            DataVersionEntry,
            StrucvarsCsqResponse,
            StrucvarsCsqQuery,
            StrucvarsGeneTranscriptEffects,
            StrucvarsSvType,
            GenomeRelease,
            StrucvarsTranscriptEffect,
            SeqvarsCsqQuery,
            SeqvarsCsqResponse,
            SeqvarsCsqResultEntry,
            Consequence,
            PutativeImpact,
            FeatureType,
            FeatureBiotype,
            Rank,
            Pos,
            Message,
            SoFeature,
            ExonAlignment,
            GenesTranscriptsListQuery,
            GenesTranscriptsListResponse,
            GenomeAlignment,
            Strand,
            Transcript,
            TranscriptBiotype,
            TranscriptTag,
        ))
    )]
    pub struct ApiDoc;
}

/// Command line arguments for "server run` command.
#[derive(clap::Parser, Debug)]
#[command(about = "Run Mehari REST API server", long_about = None)]
pub struct Args {
    /// What to annotate and which source to use.
    #[command(flatten)]
    pub sources: Sources,

    /// Transcript related settings.
    #[command(flatten)]
    pub transcript_settings: TranscriptSettings,

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

/// Main entry point for `server run` sub command.
///
/// # Errors
///
/// In the case that there is an error running the server.
pub async fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("args_common = {:?}", &args_common);
    tracing::info!("args = {:?}", &args);

    if let Some(log::Level::Trace | log::Level::Debug) = args_common.verbose.log_level() {
        std::env::set_var("RUST_LOG", "debug");
        env_logger::init_from_env(env_logger::Env::new().default_filter_or("info"));
    }

    // Load data that we need for running the server.
    tracing::info!("Loading database...");
    let before_loading = std::time::Instant::now();
    let mut data = actix_server::WebServerData::default();
    for genome_release in [GenomeRelease::Grch37, GenomeRelease::Grch38] {
        let assembly = genome_release.into();

        let annotator = setup_seqvars_annotator(&args.sources, &args.transcript_settings)?;
        let seqvars_csq_predictor = &annotator.seqvars().unwrap().predictor;
        let provider = seqvars_csq_predictor.provider.clone();
        data.provider.insert(genome_release, provider.clone());
        tracing::info!("  - building seqvars predictors");
        data.seqvars_predictors.insert(
            genome_release,
            SeqvarConsequencePredictor::new(provider.clone(), Default::default()),
        );
        tracing::info!("  - building strucvars predictors");
        data.strucvars_predictors.insert(
            genome_release,
            StrucvarConsequencePredictor::new(provider.clone(), assembly),
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
