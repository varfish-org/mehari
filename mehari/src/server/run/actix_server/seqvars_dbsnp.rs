//! Implementation of endpoint `/api/v1/seqvars/dbsnp`.

use super::{CustomError, versions::VersionsInfoResponse};
use crate::db::keys::Var;
use actix_web::{
    get,
    web::{self, Data, Json, Path},
};

/// Query parameters of the `/api/v1/seqvars/dbsnp` endpoint.
#[derive(
    Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::IntoParams, utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
pub(crate) struct DbsnpQuery {
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

/// One entry in `DbsnpResponse`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct DbsnpResultEntry {
    pub allele: String,
    pub rs_id: String,
}

/// Response of the `/api/v1/seqvars/dbsnp` endpoint.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct DbsnpResponse {
    /// Version information.
    pub version: VersionsInfoResponse,

    /// The original query parameters.
    pub query: DbsnpQuery,

    /// The resulting dbSNP record if found.
    pub result: Option<DbsnpResultEntry>,
}

/// Implementation of endpoints.
async fn handle_impl(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<DbsnpQuery>,
) -> actix_web::Result<Json<DbsnpResponse>, CustomError> {
    let DbsnpQuery {
        assembly,
        chromosome,
        position,
        reference,
        alternative,
    } = query.clone().into_inner();

    let annotator = data.dbsnp_annotators.get(&assembly).ok_or_else(|| {
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

    let result = annotation.map(|res| DbsnpResultEntry {
        allele: res.allele,
        rs_id: res.rs_id,
    });

    let response = DbsnpResponse {
        version: VersionsInfoResponse::from_web_server_data(data.into_inner().as_ref())
            .map_err(|e| CustomError::new(anyhow::anyhow!("Problem determining version: {}", e)))?,
        query: query.into_inner(),
        result,
    };

    Ok(Json(response))
}

#[allow(clippy::unused_async)]
#[get("/seqvars/dbsnp")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<DbsnpQuery>,
) -> actix_web::Result<Json<DbsnpResponse>, CustomError> {
    handle_impl(data, _path, query).await
}

#[allow(clippy::unused_async)]
#[utoipa::path(
    get,
    operation_id = "seqvarsDbsnp",
    params(
        DbsnpQuery
    ),
    responses(
        (status = 200, description = "dbSNP records.", body = DbsnpResponse),
        (status = 500, description = "Internal server error.", body = CustomError)
    )
)]
#[get("/api/v1/seqvars/dbsnp")]
async fn handle_with_openapi(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<DbsnpQuery>,
) -> actix_web::Result<Json<DbsnpResponse>, super::CustomError> {
    handle_impl(data, _path, query).await
}
