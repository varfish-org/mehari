//! Run the server.

use actix_web::ResponseError;

use crate::{annotate::seqvars::csq::ConsequencePredictor, common::GenomeRelease};

pub mod tx_csq;

#[derive(Debug)]
struct CustomError {
    err: anyhow::Error,
}

impl std::fmt::Display for CustomError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.err)
    }
}

impl CustomError {
    fn new(err: anyhow::Error) -> Self {
        CustomError { err }
    }
}

impl ResponseError for CustomError {}

/// Data structure for the web server data.
#[derive(Debug, Default)]
pub struct WebServerData {
    /// The consequence predictors for each assembly.
    pub predictors: std::collections::HashMap<GenomeRelease, ConsequencePredictor>,
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
            .service(tx_csq::handle)
            .wrap(actix_web::middleware::Logger::default())
    })
    .bind((args.listen_host.as_str(), args.listen_port))?
    .run()
    .await
}
