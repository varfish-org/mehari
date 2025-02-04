//! Implementation of endpoint `/api/v1/seqvars/clinvar`.
//!
//! Also includes the implementation of the `/seqvars/clinvar` endpoint (deprecated).

use actix_web::{
    get,
    web::{self, Data, Json, Path},
};

use crate::{annotate::seqvars::csq::VcfVariant, common::GenomeRelease};

use super::{versions::VersionsInfoResponse, CustomError};

/// Query parameters of the `/api/v1/seqvars/clinvar` endpoint.
#[derive(
    Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::IntoParams, utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
#[serde_with::skip_serializing_none]
pub(crate) struct ClinvarQuery {
    /// The assembly.
    pub genome_release: GenomeRelease,
    /// SPDI sequence.
    pub chromosome: String,
    /// SPDI position.
    pub position: u32,
    /// SPDI deletion.
    pub reference: String,
    /// SPDI insertion.
    pub alternative: String,
}

/// One entry in `ClinvarResponse`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct ClinvarResultEntry {
    pub clinvar_vcv: Vec<String>,
    pub clinvar_germline_classification: Vec<String>,
}

/// Response of the `/api/v1/seqvars/clinvar` endpoint.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct ClinvarResponse {
    /// Version information.
    pub version: VersionsInfoResponse,

    /// The original query records.
    pub query: ClinvarQuery,

    /// The resulting records for the scored genes.
    pub result: Vec<ClinvarResultEntry>,
}

/// Implementation of endpoints.
async fn handle_impl(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<ClinvarQuery>,
) -> actix_web::Result<Json<ClinvarResponse>, super::CustomError> {
    let ClinvarQuery {
        genome_release,
        chromosome,
        position,
        reference,
        alternative,
    } = query.clone().into_inner();

    let annotator = data
        .clinvar_annotators
        .get(&genome_release)
        .ok_or_else(|| {
            super::CustomError::new(anyhow::anyhow!(
                "genome release not supported: {:?}",
                &query.genome_release
            ))
        })?;

    let mut result = Vec::new();
    let g_var = VcfVariant {
        chromosome,
        position: position as i32,
        reference: reference.into(),
        alternative: alternative.into(),
    };
    let annotations = annotator
        .annotate_variant(&g_var)
        .map_err(|e| super::CustomError::new(anyhow::anyhow!("annotation failed: {}", &e)))?;
    if let Some(annotations) = annotations {
        result.push(annotations);
    }

    let result = ClinvarResponse {
        version: VersionsInfoResponse::from_web_server_data(data.into_inner().as_ref())
            .map_err(|e| CustomError::new(anyhow::anyhow!("Problem determining version: {}", e)))?,
        query: query.into_inner(),
        result,
    };

    Ok(Json(result))
}

/// Query for ClinVar information of a variant.
#[allow(clippy::unused_async)]
#[get("/seqvars/clinvar")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<ClinvarQuery>,
) -> actix_web::Result<Json<ClinvarResponse>, super::CustomError> {
    handle_impl(data, _path, query).await
}

/// Query for ClinVar information of a variant.
#[allow(clippy::unused_async)]
#[utoipa::path(
    get,
    operation_id = "seqvarsClinvar",
    params(
        ClinvarQuery
    ),
    responses(
        (status = 200, description = "Clinvar information.", body = ClinvarResponse),
        (status = 500, description = "Internal server error.", body = CustomError)
    )
)]
#[get("/api/v1/seqvars/clinvar")]
async fn handle_with_openapi(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<ClinvarQuery>,
) -> actix_web::Result<Json<ClinvarResponse>, super::CustomError> {
    handle_impl(data, _path, query).await
}
