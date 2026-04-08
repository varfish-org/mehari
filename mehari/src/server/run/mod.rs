use crate::annotate::cli::{PredictorSettings, Sources};
use crate::annotate::seqvars::csq::ConfigBuilder;
use crate::annotate::seqvars::{
    ConsequenceAnnotator, initialize_clinvar_annotators_for_assembly,
    initialize_frequency_annotators_for_assembly, load_transcript_dbs_for_assembly,
};
use crate::annotate::{
    seqvars::csq::ConsequencePredictor as SeqvarConsequencePredictor,
    strucvars::csq::ConsequencePredictor as StrucvarConsequencePredictor,
};
use crate::common::contig::ContigManager;
use crate::db::merge::merge_transcript_databases;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::Arc;
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
        CustomError, gene_txs, seqvars_clinvar, seqvars_csq, seqvars_frequencies, strucvars_csq,
        versions,
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

#[derive(Debug, Clone)]
pub struct AssemblyReference {
    pub assembly: String,
    pub path: PathBuf,
}

impl FromStr for AssemblyReference {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.split_once('=') {
            Some((assembly, path)) if !assembly.is_empty() && !path.is_empty() => {
                Ok(AssemblyReference {
                    assembly: assembly.to_lowercase(),
                    path: PathBuf::from(path),
                })
            }
            _ => Err(format!(
                "Invalid reference format `{}`. Expected `ASSEMBLY=PATH`",
                s
            )),
        }
    }
}

/// Command line arguments for "server run` command.
#[derive(clap::Parser, Debug)]
#[command(about = "Run Mehari REST API server", long_about = None)]
pub struct Args {
    /// Path to the reference genome(s), with accompanying index.
    ///
    /// Provide one reference per assembly you want to serve (e.g., GRCh37 and/or GRCh38).
    /// Example: --reference grch37=/path/to/grch37.fa
    #[arg(long, value_name = "ASSEMBLY=PATH")]
    pub reference: Vec<AssemblyReference>,

    /// Read the reference genome into memory.
    #[arg(long, requires = "reference")]
    pub in_memory_reference: bool,

    /// What to annotate and which source to use.
    #[command(flatten)]
    pub sources: Sources,

    /// Transcript related settings.
    #[command(flatten)]
    pub predictor_settings: PredictorSettings,

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
fn print_hints(args: &Args, enabled_sources: &[(String, Endpoint)]) {
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

    let prefix = format!(
        "try: http://{host}:{port}/api/v1/",
        host = args.listen_host,
        port = args.listen_port
    );
    let examples: HashMap<(String, Endpoint), Vec<&str>> = HashMap::from([
        (
            ("grch37".into(), Transcripts),
            vec![
                r#"genes/transcripts?hgnc_id=HGNC:1100&genome_build=grch37"#,
                r#"seqvars/csq?genome_release=grch37&chromosome=17&position=48275363&reference=C&alternative=A"#,
                r#"seqvars/csq?genome_release=grch37&chromosome=17&position=48275363&reference=C&alternative=A&hgnc_id=HGNC:2197"#,
                r#"strucvars/csq?genome_release=grch37&chromosome=17&start=48275360&&stop=48275370&sv_type=DEL""#,
            ],
        ),
        (
            ("grch37".into(), Frequency),
            vec![
                r#"seqvars/frequency?genome_release=grch37&chromosome=17&position=48275363&reference=C&alternative=A"#,
            ],
        ),
        (
            ("grch37".into(), Clinvar),
            vec![
                r#"seqvars/clinvar?genome_release=grch37&chromosome=17&position=48275363&reference=C&alternative=A"#,
            ],
        ),
        (
            ("grch38".into(), Transcripts),
            vec![
                r#"genes/transcripts?hgnc_id=HGNC:1100&genome_build=grch38"#,
                r#"seqvars/csq?genome_release=grch38&chromosome=2&position=26364839&reference=C&alternative=T"#,
            ],
        ),
        (
            ("grch38".into(), Frequency),
            vec![
                r#"seqvars/frequency?genome_release=grch38&chromosome=2&position=26364839&reference=C&alternative=T"#,
            ],
        ),
        (
            ("grch38".into(), Clinvar),
            vec![
                r#"seqvars/clinvar?genome_release=grch38&chromosome=2&position=26364839&reference=C&alternative=T"#,
            ],
        ),
    ]);
    for (genome_release, endpoint) in enabled_sources {
        if let Some(examples) = examples.get(&(genome_release.clone(), *endpoint)) {
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

    let reference_paths: HashMap<String, PathBuf> = args
        .reference
        .iter()
        .map(|r| (r.assembly.clone(), r.path.clone()))
        .collect();

    let associate_db_path = |path: &str| -> Option<String> {
        let path_str = path.to_lowercase();
        if path_str.contains("grch37") {
            Some("grch37".into())
        } else if path_str.contains("grch38") {
            Some("grch38".into())
        } else {
            None
        }
    };

    let mut transcript_paths: HashMap<String, Vec<String>> = HashMap::new();
    if let Some(paths) = args.sources.transcripts.as_ref() {
        for path in paths {
            if let Some(release) = associate_db_path(path) {
                transcript_paths
                    .entry(release)
                    .or_default()
                    .push(path.clone());
            } else {
                tracing::warn!(
                    "Could not determine assembly for transcript db: {}. Skipping.",
                    path
                );
            }
        }
    }

    let mut frequency_paths = HashMap::new();
    if let Some(paths) = args.sources.frequencies.as_ref() {
        for path in paths {
            if let Some(assembly) = associate_db_path(path) {
                if let Some(prev) = frequency_paths.insert(assembly.clone(), path.clone()) {
                    tracing::warn!(
                        "Duplicate frequency DB for {:?}: {} overwritten by {}",
                        assembly,
                        prev,
                        path
                    );
                }
            } else {
                tracing::warn!(
                    "Could not determine assembly for frequency db: {}. Skipping.",
                    path
                );
            }
        }
    }

    let mut clinvar_paths = HashMap::new();
    if let Some(paths) = args.sources.clinvar.as_ref() {
        for path in paths {
            if let Some(assembly) = associate_db_path(path) {
                if let Some(prev) = clinvar_paths.insert(assembly.clone(), path.clone()) {
                    tracing::warn!(
                        "Duplicate ClinVar DB for {:?}: {} overwritten by {}",
                        assembly,
                        prev,
                        path
                    );
                }
            } else {
                tracing::warn!(
                    "Could not determine assembly for clinvar db: {}. Skipping.",
                    path
                );
            }
        }
    }

    let mut all_releases: HashSet<String> = HashSet::new();
    all_releases.extend(transcript_paths.keys().cloned());
    all_releases.extend(frequency_paths.keys().cloned());
    all_releases.extend(clinvar_paths.keys().cloned());

    let mut enabled_sources = vec![];
    for assembly in all_releases.iter() {
        let contig_manager = Arc::new(ContigManager::new(&assembly));
        let reference_path = reference_paths.get(&assembly.clone());

        tracing::info!("Loading data for assembly {:?}...", &assembly);

        // Load transcript data if available for this release.
        // At the moment, this is the only type of data that can benefit from a reference FASTA.
        if let Some(tx_db_paths) = transcript_paths.get(&assembly.clone()) {
            if reference_path.is_none() {
                tracing::warn!(
                    "No reference FASTA provided for {}. Some features like HGVSg normalization will be unavailable.",
                    &assembly
                );
            }

            tracing::info!("Building seqvars predictors for {}...", &assembly);
            let tx_dbs = load_transcript_dbs_for_assembly(tx_db_paths, &assembly)?;
            if tx_dbs.is_empty() {
                tracing::warn!(
                    "No transcript databases loaded for {}, respective endpoints will be unavailable.",
                    &assembly
                );
            } else {
                let tx_db = merge_transcript_databases(tx_dbs)?;
                let annotator = ConsequenceAnnotator::from_db_and_settings(
                    tx_db,
                    reference_path.cloned(),
                    args.in_memory_reference,
                    &args.predictor_settings,
                )?;
                let config = ConfigBuilder::default()
                    .report_most_severe_consequence_by(
                        args.predictor_settings
                            .transcript_settings
                            .report_most_severe_consequence_by,
                    )
                    .transcript_source(
                        args.predictor_settings
                            .transcript_settings
                            .transcript_source,
                    )
                    .build()?;

                let provider = annotator.predictor.provider.clone();
                data.provider.insert(assembly.clone(), provider.clone());
                data.seqvars_predictors.insert(
                    assembly.clone(),
                    SeqvarConsequencePredictor::new(provider.clone(), config),
                );
                tracing::info!("Finished building seqvars predictors.");

                tracing::info!("Building strucvars predictors...");
                data.strucvars_predictors.insert(
                    assembly.clone(),
                    StrucvarConsequencePredictor::new(provider.clone()),
                );
                enabled_sources.push((assembly.clone(), Endpoint::Transcripts));
                tracing::info!(
                    "Finished building seqvars predictors for {:?}.",
                    assembly.clone()
                );
            }
        }

        if let Some(freq_path) = frequency_paths.get(&assembly.clone()) {
            tracing::info!("Loading frequency data for {:?}...", assembly);
            let annotators = initialize_frequency_annotators_for_assembly(
                std::slice::from_ref(freq_path),
                &assembly,
                contig_manager.clone(),
            )?;
            if let Some(frequency_db) = annotators.into_iter().next() {
                data.frequency_annotators
                    .insert(assembly.clone(), frequency_db);
                enabled_sources.push((assembly.clone(), Endpoint::Frequency));
                tracing::info!(
                    "... done loading frequency data for {:?}.",
                    assembly.clone()
                );
            } else {
                tracing::warn!("Could not load frequency data for {:?}.", assembly.clone());
            }
        }

        if let Some(clinvar_path) = clinvar_paths.get(&assembly.clone()) {
            tracing::info!("Loading ClinVar data for {:?}...", assembly.clone());
            let annotators = initialize_clinvar_annotators_for_assembly(
                std::slice::from_ref(clinvar_path),
                &assembly,
                contig_manager.clone(),
            )?;
            if let Some(annotator) = annotators.into_iter().next() {
                data.clinvar_annotators.insert(assembly.clone(), annotator);
                enabled_sources.push((assembly.clone(), Endpoint::Clinvar));
                tracing::info!("... done loading ClinVar data for {:?}.", assembly.clone());
            } else {
                tracing::warn!("Could not load ClinVar data for {:?}.", assembly.clone());
            }
        }
    }
    let data = actix_web::web::Data::new(data);
    tracing::info!("... done loading data {:?}", before_loading.elapsed());

    if enabled_sources.is_empty() {
        tracing::error!("No data sources loaded. Exiting.");
        return Err(anyhow::anyhow!("No data sources loaded."));
    }

    tracing::info!("Loaded the following sources: {:?}", &enabled_sources);

    // Print the server URL and some hints (the latter: unless suppressed).
    print_hints(args, &enabled_sources);
    // Launch the Actix web server.
    actix_server::main(args, data).await?;

    tracing::info!("All done. Have a nice day!");
    Ok(())
}
