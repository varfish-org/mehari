//! Implementation of endpoint `/api/v1/strucvars/csq`.
//!
//! Also includes the implementation of the `/strucvars/csq` endpoint (deprecated).

use actix_web::{
    get,
    web::{self, Data, Json, Path},
    Responder,
};

use crate::{
    annotate::strucvars::csq::{
        interface::{self, StrucvarsSvType},
        StrucvarsGeneTranscriptEffects,
    },
    common::GenomeRelease,
    server::run::actix_server::CustomError,
};

use super::versions::VersionsInfoResponse;

/// Parameters for `/strucvars/csq`.
///
#[derive(serde::Serialize, serde::Deserialize, Debug, Clone)]
#[serde(rename_all = "snake_case")]
#[serde_with::skip_serializing_none]
struct Query {
    /// The assembly.
    pub genome_release: GenomeRelease,
    /// Chromosome.
    pub chromosome: String,
    /// 1-based start position.
    pub start: i32,
    /// 1-based stop position, ignored for INS.
    pub stop: Option<i32>,
    /// The variant type to use for annotation.
    pub sv_type: String,
}

impl interface::StrucVar for Query {
    fn chrom(&self) -> String {
        self.chromosome.clone()
    }

    fn chrom2(&self) -> String {
        self.chromosome.clone()
    }

    fn start(&self) -> i32 {
        self.start
    }

    fn stop(&self) -> i32 {
        if self.sv_type() == interface::StrucvarsSvType::Ins {
            self.start
        } else {
            self.stop.unwrap_or(self.start)
        }
    }

    fn sv_type(&self) -> interface::StrucvarsSvType {
        match self.sv_type.to_uppercase().as_ref() {
            "DEL" => interface::StrucvarsSvType::Del,
            "DUP" => interface::StrucvarsSvType::Dup,
            "INS" => interface::StrucvarsSvType::Ins,
            "BND" => interface::StrucvarsSvType::Bnd,
            "INV" => interface::StrucvarsSvType::Inv,
            _ => interface::StrucvarsSvType::Del,
        }
    }

    fn strand_orientation(&self) -> interface::StrandOrientation {
        interface::StrandOrientation::NotApplicable
    }
}

/// Container for the result.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
struct Container {
    /// Version information.
    pub version: crate::common::Version,
    /// The original query records.
    pub query: Query,
    /// The resulting records for the scored genes.
    pub result: Vec<StrucvarsGeneTranscriptEffects>,
}

/// Query for consequence of a variant.
#[allow(clippy::unused_async)]
#[get("/strucvars/csq")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<Query>,
) -> actix_web::Result<impl Responder, super::CustomError> {
    let predictor = data
        .strucvars_predictors
        .get(&query.genome_release)
        .ok_or_else(|| {
            super::CustomError::new(anyhow::anyhow!(
                "genome release not supported: {:?}",
                &query.genome_release
            ))
        })?;

    let result = predictor.compute_tx_effects(&query.clone().into_inner());

    let result = Container {
        version: crate::common::Version::new(predictor.data_version()),
        query: query.into_inner(),
        result,
    };

    Ok(Json(result))
}

/// Query parameters of the `/api/v1/strucvars/csq` endpoint.
#[derive(
    Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::IntoParams, utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
#[serde_with::skip_serializing_none]
pub(crate) struct StrucvarsCsqQuery {
    /// The assembly.
    pub genome_release: GenomeRelease,
    /// Chromosome.
    pub chromosome: String,
    /// 1-based start position.
    pub start: i32,
    /// 1-based stop position, ignored for INS.
    pub stop: Option<i32>,
    /// The variant type to use for annotation.
    pub sv_type: StrucvarsSvType,
}

impl interface::StrucVar for StrucvarsCsqQuery {
    fn chrom(&self) -> String {
        self.chromosome.clone()
    }

    fn chrom2(&self) -> String {
        self.chromosome.clone()
    }

    fn start(&self) -> i32 {
        self.start
    }

    fn stop(&self) -> i32 {
        if self.sv_type() == interface::StrucvarsSvType::Ins {
            self.start
        } else {
            self.stop.unwrap_or(self.start)
        }
    }

    fn sv_type(&self) -> interface::StrucvarsSvType {
        self.sv_type
    }

    fn strand_orientation(&self) -> interface::StrandOrientation {
        interface::StrandOrientation::NotApplicable
    }
}

/// Response of the `/api/v1/strucvars/csq` endpoint.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct StrucvarsCsqResponse {
    /// Version information.
    pub version: VersionsInfoResponse,
    /// The original query record.
    pub query: StrucvarsCsqQuery,
    /// The resulting records for the affected genes.
    pub result: Vec<StrucvarsGeneTranscriptEffects>,
}

/// Query for consequence of a variant.
#[allow(clippy::unused_async)]
#[utoipa::path(
    get,
    operation_id = "strucvarsCsq",
    params(
        StrucvarsCsqQuery
    ),
    responses(
        (status = 200, description = "Strucvars consequence information.", body = StrucvarsCsqResponse),
        (status = 500, description = "Internal server error.", body = CustomError)
    )
)]
#[get("/api/v1/strucvars/csq")]
async fn handle_with_openapi(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<StrucvarsCsqQuery>,
) -> actix_web::Result<Json<StrucvarsCsqResponse>, CustomError> {
    let predictor = data
        .strucvars_predictors
        .get(&query.genome_release)
        .ok_or_else(|| {
            super::CustomError::new(anyhow::anyhow!(
                "genome release not supported: {:?}",
                &query.genome_release
            ))
        })?;

    let result = predictor.compute_tx_effects(&query.clone().into_inner());

    let result = StrucvarsCsqResponse {
        version: VersionsInfoResponse::from_web_server_data(data.into_inner().as_ref())
            .map_err(|e| CustomError::new(anyhow::anyhow!("Problem determining version: {}", e)))?,
        query: query.into_inner(),
        result,
    };

    Ok(Json(result))
}
