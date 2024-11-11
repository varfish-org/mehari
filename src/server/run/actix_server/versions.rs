use actix_web::{
    get,
    web::{self, Data, Json, Path},
};

use crate::{annotate::seqvars::provider::Provider, pbs};

use super::CustomError;

/// Assembly to be passed on the command line.
///
/// Copy from annonars with extension to derive `utoipa::ToSchema`.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    clap::ValueEnum,
    serde::Deserialize,
    serde::Serialize,
    utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Assembly {
    /// GRCh37
    Grch37,
    /// GRCh38
    Grch38,
}

impl From<biocommons_bioutils::assemblies::Assembly> for Assembly {
    fn from(assembly: biocommons_bioutils::assemblies::Assembly) -> Self {
        match assembly {
            biocommons_bioutils::assemblies::Assembly::Grch37
            | biocommons_bioutils::assemblies::Assembly::Grch37p10 => Assembly::Grch37,
            biocommons_bioutils::assemblies::Assembly::Grch38 => Assembly::Grch38,
        }
    }
}

impl From<Assembly> for pbs::txs::GenomeBuild {
    fn from(val: Assembly) -> Self {
        match val {
            Assembly::Grch37 => pbs::txs::GenomeBuild::Grch37,
            Assembly::Grch38 => pbs::txs::GenomeBuild::Grch38,
        }
    }
}

impl TryFrom<pbs::txs::GenomeBuild> for Assembly {
    type Error = anyhow::Error;

    fn try_from(value: pbs::txs::GenomeBuild) -> Result<Self, Self::Error> {
        match value {
            pbs::txs::GenomeBuild::Grch37 => Ok(Self::Grch37),
            pbs::txs::GenomeBuild::Grch38 => Ok(Self::Grch38),
            _ => Err(anyhow::anyhow!("Unsupported assembly")),
        }
    }
}

/// Software version specification.
#[derive(Clone, Debug, PartialEq, Eq, serde::Deserialize, serde::Serialize, utoipa::ToSchema)]
pub struct SoftwareVersions {
    /// Version of `mehari`.
    pub mehari: String,
    /// Version of the `hgvs` crate.
    pub hgvs_rs: String,
}

impl SoftwareVersions {
    /// Create a new `SoftwareVersions` instance.
    pub fn new() -> Result<Self, anyhow::Error> {
        let mehari = crate::built_info::PKG_VERSION.to_string();
        let hgvs_rs = crate::built_info::DEPENDENCIES
            .iter()
            .find(|(name, _)| name == &"hgvs")
            .map(|(_, version)| version.to_string())
            .ok_or_else(|| anyhow::anyhow!("Failed to find hgvs version"))?;

        Ok(Self { mehari, hgvs_rs })
    }
}

/// Specification of data version for a given genome build.
#[derive(Clone, Debug, PartialEq, Eq, serde::Deserialize, serde::Serialize, utoipa::ToSchema)]
pub struct DataVersionEntry {
    /// Assembly for which the data version is specified.
    pub genome_build: Assembly,
    /// Version of the RefSeq database, if any.
    pub version_refseq: Option<String>,
    /// Version of the Ensembl database, if any.
    pub version_ensembl: Option<String>,
}

impl DataVersionEntry {
    /// Create a new `DataVersionEntry` instance from `Provider`.`
    pub fn from_provider(provider: &Provider) -> Self {
        let genome_build = Assembly::from(provider.assembly());
        Self {
            genome_build,
            version_refseq: None,
            version_ensembl: None,
        }
    }
}

/// Response of the `/api/v1/version` endpoint.
#[derive(Clone, Debug, PartialEq, Eq, serde::Deserialize, serde::Serialize, utoipa::ToSchema)]
pub struct VersionsInfoResponse {
    /// Software versions specification.
    pub software: SoftwareVersions,
    /// Data versions specification.
    pub data: Vec<DataVersionEntry>,
}

impl VersionsInfoResponse {
    /// Create a new `VersionsInfoResponse` instance from the given `WebServerData``.
    pub fn from_web_server_data(data: &super::WebServerData) -> Result<Self, anyhow::Error> {
        let software = SoftwareVersions::new()?;
        let data = data
            .provider
            .values()
            .map(|provider| DataVersionEntry::from_provider(provider.as_ref()))
            .collect();

        Ok(Self { software, data })
    }
}

/// Query for consequence of a variant.
#[allow(clippy::unused_async)]
#[utoipa::path(
    get,
    operation_id = "versionsInfo",
    responses(
        (status = 200, description = "Version information.", body = VersionsInfoResponse),
        (status = 500, description = "Internal server error.", body = CustomError)
    )
)]
#[get("/api/v1/versionsInfo")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    _query: web::Query<()>,
) -> actix_web::Result<Json<VersionsInfoResponse>, CustomError> {
    Ok(Json(
        VersionsInfoResponse::from_web_server_data(data.into_inner().as_ref())
            .map_err(|e| CustomError::new(anyhow::anyhow!("Problem determining version: {}", e)))?,
    ))
}
