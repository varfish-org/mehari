use crate::annotate::cli::{Sources, TranscriptSettings};
use crate::annotate::seqvars::csq::ConfigBuilder;
use crate::annotate::seqvars::{
    initialize_clinvar_annotators_for_assembly, initialize_frequency_annotators_for_assembly,
    load_transcript_dbs_for_assembly, ConsequenceAnnotator,
};
use crate::db::merge::merge_transcript_databases;
use crate::{
    annotate::{
        seqvars::csq::ConsequencePredictor as SeqvarConsequencePredictor,
        strucvars::csq::ConsequencePredictor as StrucvarConsequencePredictor,
    },
    common::GenomeRelease,
};
use anyhow::anyhow;
use clap::ValueEnum;
use std::collections::HashMap;
use std::path::PathBuf;
use strum::EnumString;

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
    use crate::server::run::actix_server::seqvars_clinvar::{
        ClinvarQuery, ClinvarResponse, ClinvarResultEntry,
    };
    use crate::server::run::actix_server::seqvars_csq::{
        SeqvarsCsqQuery, SeqvarsCsqResponse, SeqvarsCsqResultEntry,
    };
    use crate::server::run::actix_server::seqvars_frequencies::{
        AutosomalResultEntry, FrequencyQuery, FrequencyResponse, FrequencyResultEntry,
        GonosomalResultEntry, MitochondrialResultEntry,
    };
    use crate::server::run::actix_server::strucvars_csq::{
        StrucvarsCsqQuery, StrucvarsCsqResponse,
    };
    use crate::server::run::actix_server::versions::{
        Assembly, DataVersionEntry, SoftwareVersions, VersionsInfoResponse,
    };

    use super::actix_server::{
        gene_txs, seqvars_clinvar, seqvars_csq, seqvars_frequencies, strucvars_csq, versions,
        CustomError,
    };

    /// Utoipa-based `OpenAPI` generation helper.
    #[derive(utoipa::OpenApi)]
    #[openapi(
        paths(
            versions::handle,
            gene_txs::handle_with_openapi,
            seqvars_csq::handle_with_openapi,
            strucvars_csq::handle_with_openapi,
            seqvars_frequencies::handle_with_openapi,
            seqvars_clinvar::handle_with_openapi,
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
            ClinvarQuery,
            ClinvarResponse,
            ClinvarResultEntry,
        ))
    )]
    pub struct ApiDoc;
}

/// Command line arguments for "server run` command.
#[derive(clap::Parser, Debug)]
#[command(about = "Run Mehari REST API server", long_about = None)]
pub struct Args {
    /// Path to the reference genome(s), with accompanying index.
    ///
    /// Note that (at the moment) the reference path must contain "GRCh37" or "GRCh38"
    /// to be able to match the reference to the transcript db correctly.
    #[arg(long)]
    pub reference: Vec<PathBuf>,

    /// Read the reference genome into memory.
    #[arg(long, requires = "reference")]
    pub in_memory_reference: bool,

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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, EnumString)]
enum Endpoint {
    Transcripts,
    Frequency,
    Clinvar,
}

/// Print some hints via `tracing::info!`.
fn print_hints(args: &Args, enabled_sources: &[(GenomeRelease, Endpoint)]) {
    tracing::info!(
        "Launching server main on http://{}:{} ...",
        args.listen_host.as_str(),
        args.listen_port
    );

    // Short-circuit if no hints are to be
    if args.suppress_hints {
        return;
    }

    use Endpoint::*;
    use GenomeRelease::*;

    let prefix = format!(
        "try: http://{host}:{port}/api/v1/",
        host = args.listen_host,
        port = args.listen_port
    );
    let examples: HashMap<(GenomeRelease, Endpoint), Vec<&str>> = HashMap::from([
        (
            (Grch37, Transcripts),
            vec![
                r#"genes/transcripts?hgnc_id=HGNC:1100&genome_build=grch37"#,
                r#"seqvars/csq?genome_release=grch37&chromosome=17&position=48275363&reference=C&alternative=A"#,
                r#"seqvars/csq?genome_release=grch37&chromosome=17&position=48275363&reference=C&alternative=A&hgnc_id=HGNC:2197"#,
                r#"strucvars/csq?genome_release=grch37&chromosome=17&start=48275360&&stop=48275370&sv_type=DEL""#,
            ],
        ),
        (
            (Grch37, Frequency),
            vec![
                r#"seqvars/frequency?genome_release=grch37&chromosome=17&position=48275363&reference=C&alternative=A"#,
            ],
        ),
        (
            (Grch37, Clinvar),
            vec![
                r#"seqvars/clinvar?genome_release=grch37&chromosome=17&position=48275363&reference=C&alternative=A"#,
            ],
        ),
        (
            (Grch38, Transcripts),
            vec![
                r#"genes/transcripts?hgnc_id=HGNC:1100&genome_build=grch38"#,
                r#"seqvars/csq?genome_release=grch38&chromosome=2&position=26364839&reference=C&alternative=T"#,
            ],
        ),
        (
            (Grch38, Frequency),
            vec![
                r#"seqvars/frequency?genome_release=grch38&chromosome=2&position=26364839&reference=C&alternative=T"#,
            ],
        ),
        (
            (Grch38, Clinvar),
            vec![
                r#"seqvars/clinvar?genome_release=grch38&chromosome=2&position=26364839&reference=C&alternative=T"#,
            ],
        ),
    ]);
    for (genome_release, endpoint) in enabled_sources {
        if let Some(examples) = examples.get(&(*genome_release, *endpoint)) {
            for example in examples {
                tracing::info!("{}{}", prefix, example);
            }
        }
    }
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
        // SAFETY: This environment variable is set during server initialization,
        // before any worker threads are spawned. At this point, only the main thread
        // is running, making this operation thread-safe.
        unsafe { std::env::set_var("RUST_LOG", "debug") };
        env_logger::init_from_env(env_logger::Env::new().default_filter_or("info"));
    }

    // Load data that we need for running the server.
    tracing::info!("Loading data...");
    let before_loading = std::time::Instant::now();
    let mut data = actix_server::WebServerData::default();

    let mut enabled_sources = vec![];
    use Endpoint::*;

    if !(0..=GenomeRelease::value_variants().len()).contains(&args.reference.len()) {
        return Err(anyhow!("Supplied number of reference genomes does not match the number of supported genome releases ({:?}).",
        GenomeRelease::value_variants()));
    }

    for genome_release in GenomeRelease::value_variants().iter().copied() {
        tracing::info!("Loading genome release {:?}", genome_release);
        let assembly = genome_release.into();
        let reference = {
            args.reference
                .iter()
                .filter_map(|reference_path| {
                    if reference_path
                        .file_name()
                        .map(|f| {
                            f.to_string_lossy()
                                .to_lowercase()
                                .contains(&genome_release.name().to_lowercase())
                        })
                        .unwrap_or(false)
                    {
                        tracing::info!(
                            "Found reference genome for {:?}: {}",
                            genome_release,
                            reference_path.display()
                        );
                        Some(reference_path.clone())
                    } else {
                        None
                    }
                })
                .next()
        };

        if !args.reference.is_empty() && reference.is_none() {
            tracing::warn!("No reference genome found for {:?}.", genome_release);
        }

        if let Some(tx_db_paths) = args.sources.transcripts.as_ref() {
            tracing::info!("Building seqvars predictors...");
            let tx_dbs = load_transcript_dbs_for_assembly(tx_db_paths, Some(assembly))?;
            if tx_dbs.is_empty() {
                tracing::warn!(
                    "No transcript databases loaded, respective endpoint will be unavailable."
                );
            } else {
                let tx_db = merge_transcript_databases(tx_dbs)?;
                let annotator = ConsequenceAnnotator::from_db_and_settings(
                    tx_db,
                    reference.as_ref(),
                    args.in_memory_reference,
                    &args.transcript_settings,
                )?;
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
                tracing::info!("Finished building seqvars predictors.");

                tracing::info!("Building strucvars predictors...");
                data.strucvars_predictors.insert(
                    genome_release,
                    StrucvarConsequencePredictor::new(provider.clone(), assembly),
                );
                tracing::info!("Finished building strucvars predictors.");
                enabled_sources.push((genome_release, Transcripts));
            }
        } else {
            tracing::warn!(
                "No predictors for genome release {:?}, respective endpoint will be unavailable.",
                genome_release
            );
        }

        if let Some(paths) = args.sources.frequencies.as_ref() {
            let annotators = initialize_frequency_annotators_for_assembly(paths, Some(assembly))?;

            match annotators.len() {
                0 => tracing::warn!(
                    "No frequency databases loaded, respective endpoint will be unavailable."
                ),
                1 => {
                    let frequency_db = annotators.into_iter().next().unwrap();
                    data.frequency_annotators
                        .insert(genome_release, frequency_db);
                    enabled_sources.push((genome_release, Frequency));
                }
                _ => tracing::warn!(
                    "Multiple frequency databases loaded. This is not supported. The respective endpoint will be unavailable."
                ),
            }
        }

        if let Some(paths) = args.sources.clinvar.as_ref() {
            let annotators = initialize_clinvar_annotators_for_assembly(paths, Some(assembly))?;

            match annotators.len() {
                0 => tracing::warn!(
                    "No clinvar databases loaded, respective endpoint will be unavailable."
                ),
                1 => {
                    let annotator = annotators.into_iter().next().unwrap();
                    data.clinvar_annotators.insert(genome_release, annotator);
                    enabled_sources.push((genome_release, Clinvar));
                }
                _ => tracing::warn!(
                    "Multiple clinvar databases specified. This is not supported. The respective endpoint will be unavailable."
                ),
            }
        }
    }
    let data = actix_web::web::Data::new(data);
    tracing::info!("... done loading data {:?}", before_loading.elapsed());

    tracing::info!("Loaded the following sources: {:?}", &enabled_sources);

    // Print the server URL and some hints (the latter: unless suppressed).
    print_hints(args, &enabled_sources);
    // Launch the Actix web server.
    actix_server::main(args, data).await?;

    tracing::info!("All done. Have a nice day!");
    Ok(())
}
