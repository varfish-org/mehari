//! Run the server.

use std::sync::Arc;

use actix_web::ResponseError;
use utoipa::OpenApi as _;

use crate::annotate::seqvars::clinvar::ClinvarAnnotator;
use crate::annotate::seqvars::consequence::logic::ConsequencePredictor;
use crate::annotate::seqvars::frequency::FrequencyAnnotator;
use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::annotate::strucvars::csq::ConsequencePredictor as StrucvarConsequencePredictor;

pub mod gene_txs;
pub mod seqvars_cadd;
pub mod seqvars_clinvar;
pub mod seqvars_csq;
pub mod seqvars_dbsnp;
pub mod seqvars_frequencies;
pub mod seqvars_spliceai;
pub mod strucvars_csq;
pub mod versions;

#[derive(Debug, serde::Serialize, utoipa::ToSchema)]
pub struct CustomError {
    err: String,
    #[serde(skip)]
    status: actix_web::http::StatusCode,
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
            status: actix_web::http::StatusCode::INTERNAL_SERVER_ERROR,
        }
    }

    /// An error in the query, with status 400.
    fn bad_request(err: anyhow::Error) -> Self {
        CustomError {
            status: actix_web::http::StatusCode::BAD_REQUEST,
            ..Self::new(err)
        }
    }
}

/// A query that does not deserialize, as a bad request with a JSON body like every other error.
fn query_error(
    err: actix_web::error::QueryPayloadError,
    _req: &actix_web::HttpRequest,
) -> actix_web::Error {
    CustomError::bad_request(anyhow::anyhow!("{err}")).into()
}

impl ResponseError for CustomError {
    fn status_code(&self) -> actix_web::http::StatusCode {
        self.status
    }

    /// The error as JSON, as the OpenAPI schema says.
    fn error_response(&self) -> actix_web::HttpResponse {
        actix_web::HttpResponse::build(self.status).json(self)
    }
}

/// Data structure for the web server data.
#[derive(Default, derivative::Derivative)]
#[derivative(Debug)]
pub struct WebServerData {
    /// `MehariProvider` to provide the transcript info.
    #[derivative(Debug = "ignore")]
    pub provider: std::collections::HashMap<String, Arc<MehariProvider>>,

    /// The sequence variant consequence predictors for each assembly.
    pub seqvars_predictors: std::collections::HashMap<String, ConsequencePredictor>,

    /// The structural variant consequence predictors for each assembly.
    pub strucvars_predictors: std::collections::HashMap<String, StrucvarConsequencePredictor>,

    /// The frequency annotators for each assembly.
    pub frequency_annotators: std::collections::HashMap<String, FrequencyAnnotator>,

    /// The clinvar annotators for each assembly.
    pub clinvar_annotators: std::collections::HashMap<String, ClinvarAnnotator>,

    /// The CADD annotators for each assembly.
    pub cadd_annotators:
        std::collections::HashMap<String, crate::annotate::seqvars::cadd::CaddAnnotator>,

    /// The SpliceAI annotators for each assembly.
    pub spliceai_annotators:
        std::collections::HashMap<String, crate::annotate::seqvars::spliceai::SpliceAiAnnotator>,

    /// The dbSNP annotators for each assembly.
    pub dbsnp_annotators:
        std::collections::HashMap<String, crate::annotate::seqvars::dbsnp::DbsnpAnnotator>,
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
            .app_data(actix_web::web::QueryConfig::default().error_handler(query_error))
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
            .service(seqvars_cadd::handle)
            .service(seqvars_cadd::handle_with_openapi)
            .service(seqvars_spliceai::handle)
            .service(seqvars_spliceai::handle_with_openapi)
            .service(seqvars_dbsnp::handle)
            .service(seqvars_dbsnp::handle_with_openapi)
            .service(versions::handle)
            .service(
                utoipa_swagger_ui::SwaggerUi::new("/swagger-ui/{_:.*}")
                    .url("/api-docs/openapi.json", super::openapi::ApiDoc::openapi()),
            )
            .wrap(actix_web::middleware::Logger::default())
    })
    .bind((args.listen_host.as_str(), args.listen_port))?
    .workers(args.threads.get())
    .run()
    .await
}

#[cfg(test)]
mod test {
    /// A query that does not deserialize gives status 400 with a JSON body, as the schema says.
    #[actix_web::test]
    async fn query_error_is_json() -> Result<(), anyhow::Error> {
        let err = actix_web::web::Query::<super::seqvars_csq::SeqvarsCsqQuery>::from_query(
            "assembly=grch37&chromosome=17&position=x&reference=N&alternative=A",
        )
        .err()
        .ok_or_else(|| anyhow::anyhow!("no error"))?;
        let request = actix_web::test::TestRequest::default().to_http_request();
        let error = super::query_error(err, &request);
        let response = error.error_response();
        assert_eq!(response.status(), actix_web::http::StatusCode::BAD_REQUEST);
        assert_eq!(
            response.headers().get("content-type").map(|v| v.as_bytes()),
            Some(b"application/json".as_slice())
        );
        Ok(())
    }
}
