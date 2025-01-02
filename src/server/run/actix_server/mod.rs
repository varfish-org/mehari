//! Run the server.

use std::sync::Arc;

use actix_web::ResponseError;
use utoipa::OpenApi as _;

use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::annotate::seqvars::{ClinvarAnnotator, FrequencyAnnotator};
use crate::annotate::strucvars::csq::ConsequencePredictor as StrucvarConsequencePredictor;
use crate::{annotate::seqvars::csq::ConsequencePredictor, common::GenomeRelease};

pub mod gene_txs;
pub mod seqvars_clinvar;
pub mod seqvars_csq;
pub mod seqvars_frequencies;
pub mod strucvars_csq;
pub mod versions;

#[derive(Debug, serde::Serialize, utoipa::ToSchema)]
pub struct CustomError {
    err: String,
}

impl std::fmt::Display for CustomError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.err)
    }
}

impl CustomError {
    fn new(err: anyhow::Error) -> Self {
        CustomError {
            err: err.to_string(),
        }
    }
}

impl ResponseError for CustomError {}

/// Data structure for the web server data.
#[derive(Default, derivative::Derivative)]
#[derivative(Debug)]
pub struct WebServerData {
    /// `MehariProvider` to provide the transcript info.
    #[derivative(Debug = "ignore")]
    pub provider: std::collections::HashMap<GenomeRelease, Arc<MehariProvider>>,

    /// The sequence variant consequence predictors for each assembly.
    pub seqvars_predictors: std::collections::HashMap<GenomeRelease, ConsequencePredictor>,

    /// The structural variant consequence predictors for each assembly.
    pub strucvars_predictors:
        std::collections::HashMap<GenomeRelease, StrucvarConsequencePredictor>,

    /// The frequency annotators for each assembly.
    pub frequency_annotators: std::collections::HashMap<GenomeRelease, FrequencyAnnotator>,

    /// The clinvar annotators for each assembly.
    pub clinvar_annotators: std::collections::HashMap<GenomeRelease, ClinvarAnnotator>,
}

/// Main entry point for running the REST server.
#[allow(clippy::unused_async)]
pub async fn main(
    args: &super::Args,
    data: actix_web::web::Data<WebServerData>,
) -> std::io::Result<()> {
    actix_web::HttpServer::new(move || {
        actix_web::App::new()
            .app_data(data.clone())
            .service(gene_txs::handle)
            .service(gene_txs::handle_with_openapi)
            .service(seqvars_csq::handle)
            .service(seqvars_csq::handle_with_openapi)
            .service(strucvars_csq::handle)
            .service(strucvars_csq::handle_with_openapi)
            .service(seqvars_frequencies::handle)
            .service(seqvars_frequencies::handle_with_openapi)
            .service(seqvars_clinvar::handle)
            .service(seqvars_clinvar::handle_with_openapi)
            .service(versions::handle)
            .service(
                utoipa_swagger_ui::SwaggerUi::new("/swagger-ui/{_:.*}")
                    .url("/api-docs/openapi.json", super::openapi::ApiDoc::openapi()),
            )
            .wrap(actix_web::middleware::Logger::default())
    })
    .bind((args.listen_host.as_str(), args.listen_port))?
    .run()
    .await
}
