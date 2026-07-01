//! Implementation of endpoint `/api/v1/seqvars/cadd`.

use super::{CustomError, versions::VersionsInfoResponse};
use crate::db::keys::Var;
use actix_web::{
    get,
    web::{self, Data, Json, Path},
};

/// Query parameters of the `/api/v1/seqvars/cadd` endpoint.
#[derive(
    Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::IntoParams, utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
pub(crate) struct CaddQuery {
    /// The assembly.
    pub assembly: String,
    /// Chromosome name.
    pub chromosome: String,
    /// 1-based position.
    pub position: u32,
    /// Reference allele.
    pub reference: String,
    /// Alternative allele.
    pub alternative: String,
}

/// One entry in `CaddResponse`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct CaddResultEntry {
    pub raw_score: f32,
    pub phred: f32,
}

/// Response of the `/api/v1/seqvars/cadd` endpoint.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct CaddResponse {
    /// Version information.
    pub version: VersionsInfoResponse,

    /// The original query parameters.
    pub query: CaddQuery,

    /// The resulting CADD score if found.
    pub result: Option<CaddResultEntry>,
}

/// Implementation of endpoints.
async fn handle_impl(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<CaddQuery>,
) -> actix_web::Result<Json<CaddResponse>, CustomError> {
    let CaddQuery {
        assembly,
        chromosome,
        position,
        reference,
        alternative,
    } = query.clone().into_inner();

    let annotator = data.cadd_annotators.get(&assembly).ok_or_else(|| {
        CustomError::new(anyhow::anyhow!(
            "genome assembly not supported: {:?}",
            &query.assembly
        ))
    })?;

    let var = Var {
        chrom: chromosome,
        pos: position as i32,
        reference,
        alternative,
    };

    let annotation = annotator
        .annotate(&var)
        .map_err(|e| super::CustomError::new(anyhow::anyhow!("annotation failed: {}", &e)))?;

    let result = annotation.map(|res| CaddResultEntry {
        raw_score: res.raw_score,
        phred: res.phred,
    });

    let response = CaddResponse {
        version: VersionsInfoResponse::from_web_server_data(data.into_inner().as_ref())
            .map_err(|e| CustomError::new(anyhow::anyhow!("Problem determining version: {}", e)))?,
        query: query.into_inner(),
        result,
    };

    Ok(Json(response))
}

#[allow(clippy::unused_async)]
#[get("/seqvars/cadd")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<CaddQuery>,
) -> actix_web::Result<Json<CaddResponse>, CustomError> {
    handle_impl(data, _path, query).await
}

#[allow(clippy::unused_async)]
#[utoipa::path(
    get,
    operation_id = "seqvarsCadd",
    params(
        CaddQuery
    ),
    responses(
        (status = 200, description = "CADD scores.", body = CaddResponse),
        (status = 500, description = "Internal server error.", body = CustomError)
    )
)]
#[get("/api/v1/seqvars/cadd")]
async fn handle_with_openapi(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<CaddQuery>,
) -> actix_web::Result<Json<CaddResponse>, super::CustomError> {
    handle_impl(data, _path, query).await
}
