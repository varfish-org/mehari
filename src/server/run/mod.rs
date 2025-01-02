use crate::annotate::cli::{Sources, TranscriptSettings};
use crate::annotate::seqvars::csq::ConfigBuilder;
use crate::annotate::seqvars::{
    load_transcript_dbs_for_assembly, ConsequenceAnnotator, FrequencyAnnotator,
};
use crate::db::merge::merge_transcript_databases;
use crate::{
    annotate::{
        seqvars::csq::ConsequencePredictor as SeqvarConsequencePredictor,
        strucvars::csq::ConsequencePredictor as StrucvarConsequencePredictor,
    },
    common::GenomeRelease,
};
use clap::ValueEnum;

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
    use crate::server::run::actix_server::frequencies::{
        AutosomalResultEntry, FrequencyQuery, FrequencyResponse, FrequencyResultEntry,
        GonosomalResultEntry, MitochondrialResultEntry,
    };
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

    use super::actix_server::{
        frequencies, gene_txs, seqvars_csq, strucvars_csq, versions, CustomError,
    };

    /// Utoipa-based `OpenAPI` generation helper.
    #[derive(utoipa::OpenApi)]
    #[openapi(
        paths(
            versions::handle,
            gene_txs::handle_with_openapi,
            seqvars_csq::handle_with_openapi,
            strucvars_csq::handle_with_openapi,
            strucvars_csq::handle_with_openapi,
            frequencies::handle_with_openapi,
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
            FrequencyQuery,
            FrequencyResponse,
            FrequencyResultEntry,
            AutosomalResultEntry,
            GonosomalResultEntry,
            MitochondrialResultEntry,
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
        "  try: http://{}:{}/api/v1/genes/transcripts?hgnc_id=HGNC:1100&genome_build=grch37&\
        genomeBuild=GENOME_BUILD_GRCH37",
        args.listen_host.as_str(),
        args.listen_port
    );
    // The endpoint `/tx/csq` to compute the consequence of a variant; without and with filtering
    // for HGNC gene ID.
    tracing::info!(
        "  try: http://{}:{}/api/v1/seqvars/csq?genome_release=grch37\
        &chromosome=17&position=48275363&reference=C&alternative=A",
        args.listen_host.as_str(),
        args.listen_port
    );
    tracing::info!(
        "  try: http://{}:{}/api/v1/seqvars/csq?genome_release=grch37\
        &chromosome=17&position=48275363&reference=C&alternative=A&hgnc_id=HGNC:2197",
        args.listen_host.as_str(),
        args.listen_port
    );
    // The endpoint `/strucvars/csq` computes the consequence of an SV.
    tracing::info!(
        "  try: http://{}:{}/api/v1/strucvars/csq?genome_release=grch37\
        &chromosome=17&start=48275360&&stop=48275370&sv_type=DEL",
        args.listen_host.as_str(),
        args.listen_port
    );
    // The endpoint `/structvars/csq` computes the consequence of an SV.
    tracing::info!(
        "  try: http://{}:{}/api/v1/strucvars/csq?genome_release=grch37\
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
    for genome_release in GenomeRelease::value_variants().iter().copied() {
        tracing::info!("  - loading genome release {:?}", genome_release);
        let assembly = genome_release.into();

        if let Some(tx_db_paths) = args.sources.transcripts.as_ref() {
            tracing::info!("  - building seqvars predictors");
            let tx_dbs = load_transcript_dbs_for_assembly(tx_db_paths, Some(assembly))?;
            if tx_dbs.is_empty() {
                tracing::warn!(
                    "No transcript databases loaded, respective endpoint will be unavailable."
                );
            } else {
                let tx_db = merge_transcript_databases(tx_dbs)?;
                let annotator =
                    ConsequenceAnnotator::from_db_and_settings(tx_db, &args.transcript_settings)?;
                let config = ConfigBuilder::default()
                    .report_most_severe_consequence_by(
                        args.transcript_settings.report_most_severe_consequence_by,
                    )
                    .transcript_source(args.transcript_settings.transcript_source)
                    .build()?;

                let provider = annotator.predictor.provider.clone();
                data.provider.insert(genome_release, provider.clone());
                data.seqvars_predictors.insert(
                    genome_release,
                    SeqvarConsequencePredictor::new(provider.clone(), config),
                );

                tracing::info!("  - building strucvars predictors");
                data.strucvars_predictors.insert(
                    genome_release,
                    StrucvarConsequencePredictor::new(provider.clone(), assembly),
                );
            }
        } else {
            tracing::warn!(
                "No predictors for genome release {:?}, respective endpoint will be unavailable.",
                genome_release
            );
        }

        if let Some(path) = args.sources.frequencies.as_ref() {
            if path
                .to_ascii_lowercase()
                .contains(&genome_release.name().to_ascii_lowercase())
            {
                if let Ok(freq) = FrequencyAnnotator::from_path(path) {
                    data.frequency_annotators.insert(genome_release, freq);
                } else {
                    tracing::warn!(
                        "Failed to load frequencies predictor for genome release {:?}, respective endpoint will be unavailable.",
                        genome_release
                    );
                }
            } else {
                tracing::warn!(
                    "No frequencies predictor for genome release {:?}, respective endpoint will be unavailable.",
                    genome_release
                );
            }
        } else {
            tracing::warn!(
                "No frequencies predictor for genome release {:?}, respective endpoint will be unavailable.",
                genome_release
            );
        }
    }
    let data = actix_web::web::Data::new(data);
    tracing::info!("... done loading data {:?}", before_loading.elapsed());

    // Print the server URL and some hints (the latter: unless suppressed).
    print_hints(args);
    // Launch the Actix web server.
    actix_server::main(args, data).await?;

    tracing::info!("All done. Have a nice day!");
    Ok(())
}
