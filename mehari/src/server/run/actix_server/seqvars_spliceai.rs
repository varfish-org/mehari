//! Implementation of endpoint `/api/v1/seqvars/spliceai`.

use super::{CustomError, versions::VersionsInfoResponse};
use crate::db::keys::Var;
use actix_web::{
    get,
    web::{self, Data, Json, Path},
};

/// Query parameters of the `/api/v1/seqvars/spliceai` endpoint.
#[derive(
    Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::IntoParams, utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
pub(crate) struct SpliceAiQuery {
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

/// One prediction entry in `SpliceAiResponse`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct SpliceAiPredictionEntry {
    pub allele: String,
    pub symbol: String,
    pub ds_ag: f32,
    pub ds_al: f32,
    pub ds_dg: f32,
    pub ds_dl: f32,
    pub dp_ag: i32,
    pub dp_al: i32,
    pub dp_dg: i32,
    pub dp_dl: i32,
}

/// Response of the `/api/v1/seqvars/spliceai` endpoint.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct SpliceAiResponse {
    /// Version information.
    pub version: VersionsInfoResponse,

    /// The original query parameters.
    pub query: SpliceAiQuery,

    /// The resulting predictions.
    pub result: Vec<SpliceAiPredictionEntry>,
}

/// Implementation of endpoints.
async fn handle_impl(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<SpliceAiQuery>,
) -> actix_web::Result<Json<SpliceAiResponse>, CustomError> {
    let SpliceAiQuery {
        assembly,
        chromosome,
        position,
        reference,
        alternative,
    } = query.clone().into_inner();

    let annotator = data.spliceai_annotators.get(&assembly).ok_or_else(|| {
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

    let result = match annotation {
        Some(res) => res
            .predictions
            .into_iter()
            .map(|p| SpliceAiPredictionEntry {
                allele: p.allele,
                symbol: p.symbol,
                ds_ag: p.ds_ag,
                ds_al: p.ds_al,
                ds_dg: p.ds_dg,
                ds_dl: p.ds_dl,
                dp_ag: p.dp_ag,
                dp_al: p.dp_al,
                dp_dg: p.dp_dg,
                dp_dl: p.dp_dl,
            })
            .collect(),
        None => Vec::new(),
    };

    let response = SpliceAiResponse {
        version: VersionsInfoResponse::from_web_server_data(data.into_inner().as_ref())
            .map_err(|e| CustomError::new(anyhow::anyhow!("Problem determining version: {}", e)))?,
        query: query.into_inner(),
        result,
    };

    Ok(Json(response))
}

#[allow(clippy::unused_async)]
#[get("/seqvars/spliceai")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<SpliceAiQuery>,
) -> actix_web::Result<Json<SpliceAiResponse>, CustomError> {
    handle_impl(data, _path, query).await
}

#[allow(clippy::unused_async)]
#[utoipa::path(
    get,
    operation_id = "seqvarsSpliceai",
    params(
        SpliceAiQuery
    ),
    responses(
        (status = 200, description = "SpliceAI predictions.", body = SpliceAiResponse),
        (status = 500, description = "Internal server error.", body = CustomError)
    )
)]
#[get("/api/v1/seqvars/spliceai")]
async fn handle_with_openapi(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<SpliceAiQuery>,
) -> actix_web::Result<Json<SpliceAiResponse>, super::CustomError> {
    handle_impl(data, _path, query).await
}
