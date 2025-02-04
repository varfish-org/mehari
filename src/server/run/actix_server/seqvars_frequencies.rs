//! Implementation of endpoint `/api/v1/seqvars/frequency`.
//!
//! Also includes the implementation of the `/seqvars/frequency` endpoint (deprecated).

use actix_web::{
    get,
    web::{self, Data, Json, Path},
};

use crate::{annotate::seqvars::csq::VcfVariant, common::GenomeRelease};

use super::{versions::VersionsInfoResponse, CustomError};

/// Query parameters of the `/api/v1/seqvars/frequency` endpoint.
#[derive(
    Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::IntoParams, utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
#[serde_with::skip_serializing_none]
pub(crate) struct FrequencyQuery {
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

/// One entry in `FrequencyResponse`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) enum FrequencyResultEntry {
    Autosomal(AutosomalResultEntry),
    Gonosomal(GonosomalResultEntry),
    Mitochondrial(MitochondrialResultEntry),
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct AutosomalResultEntry {
    pub gnomad_exomes_an: u32,

    pub gnomad_exomes_hom: u32,

    pub gnomad_exomes_het: u32,

    pub gnomad_genomes_an: u32,

    pub gnomad_genomes_hom: u32,

    pub gnomad_genomes_het: u32,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct GonosomalResultEntry {
    pub gnomad_exomes_an: u32,

    pub gnomad_exomes_hom: u32,

    pub gnomad_exomes_het: u32,

    pub gnomad_exomes_hemi: u32,

    pub gnomad_genomes_an: u32,

    pub gnomad_genomes_hom: u32,

    pub gnomad_genomes_het: u32,

    pub gnomad_genomes_hemi: u32,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct MitochondrialResultEntry {
    pub helix_an: u32,

    pub helix_hom: u32,

    pub helix_het: u32,

    pub gnomad_genomes_an: u32,

    pub gnomad_genomes_hom: u32,

    pub gnomad_genomes_het: u32,
}

/// Response of the `/api/v1/seqvars/frequency` endpoint.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct FrequencyResponse {
    /// Version information.
    pub version: VersionsInfoResponse,
    /// The original query records.
    pub query: FrequencyQuery,
    /// The resulting records for the scored genes.
    pub result: Vec<FrequencyResultEntry>,
}

/// Implementation of endpoints.
async fn handle_impl(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<FrequencyQuery>,
) -> actix_web::Result<Json<FrequencyResponse>, super::CustomError> {
    let FrequencyQuery {
        genome_release,
        chromosome,
        position,
        reference,
        alternative,
    } = query.clone().into_inner();

    let annotator = data
        .frequency_annotators
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
    let frequencies = annotator
        .annotate_variant(&g_var)
        .map_err(|e| super::CustomError::new(anyhow::anyhow!("annotation failed: {}", &e)))?;
    if let Some(frequencies) = frequencies {
        result.push(frequencies);
    }

    let result = FrequencyResponse {
        version: VersionsInfoResponse::from_web_server_data(data.into_inner().as_ref())
            .map_err(|e| CustomError::new(anyhow::anyhow!("Problem determining version: {}", e)))?,
        query: query.into_inner(),
        result,
    };

    Ok(Json(result))
}

/// Query for gnomAD frequencies of a variant.
#[allow(clippy::unused_async)]
#[get("/seqvars/frequency")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<FrequencyQuery>,
) -> actix_web::Result<Json<FrequencyResponse>, super::CustomError> {
    handle_impl(data, _path, query).await
}

/// Query for gnomAD frequencies of a variant.
#[allow(clippy::unused_async)]
#[utoipa::path(
    get,
    operation_id = "seqvarsFrequency",
    params(
        FrequencyQuery
    ),
    responses(
        (status = 200, description = "Frequency information.", body = FrequencyResponse),
        (status = 500, description = "Internal server error.", body = CustomError)
    )
)]
#[get("/api/v1/seqvars/frequency")]
async fn handle_with_openapi(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<FrequencyQuery>,
) -> actix_web::Result<Json<FrequencyResponse>, super::CustomError> {
    handle_impl(data, _path, query).await
}
