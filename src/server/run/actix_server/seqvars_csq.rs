//! Implementation of endpoint `/api/v1/seqvars/csq`.
//!
//! Also includes the implementation of the `/seqvars/csq` endpoint (deprecated).

use actix_web::{
    get,
    web::{self, Data, Json, Path},
};

use crate::{
    annotate::seqvars::{
        ann::{
            AnnField, Consequence, FeatureBiotype, FeatureType, Message, Pos, PutativeImpact, Rank,
        },
        csq::VcfVariant,
    },
    common::GenomeRelease,
};

use super::{versions::VersionsInfoResponse, CustomError};

/// Query parameters of the `/api/v1/seqvars/csq` endpoint.
#[derive(
    Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::IntoParams, utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
#[serde_with::skip_serializing_none]
pub(crate) struct SeqvarsCsqQuery {
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
    /// Optionally, the HGNC ID of the gene to limit to.
    pub hgnc_id: Option<String>,
}

/// One entry in `SeqvarsCsqResponse`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct SeqvarsCsqResultEntry {
    /// The consequences of the allele.
    pub consequences: Vec<Consequence>,
    /// The putative impact.
    pub putative_impact: PutativeImpact,
    /// The gene symbol.
    pub gene_symbol: String,
    /// The gene identifier.
    pub gene_id: String,
    /// The feature type.
    pub feature_type: FeatureType,
    /// The feature identifier.
    pub feature_id: String,
    /// The feature biotype.
    pub feature_biotype: FeatureBiotype,
    /// The feature tags.
    pub feature_tag: Vec<FeatureBiotype>,
    /// The exon / intron rank.
    pub rank: Option<Rank>,
    /// HGVS c. notation.
    pub hgvs_t: Option<String>,
    /// HGVS p. notation.
    pub hgvs_p: Option<String>,
    /// cDNA position.
    pub tx_pos: Option<Pos>,
    /// CDS position.
    pub cds_pos: Option<Pos>,
    /// Protein position.
    pub protein_pos: Option<Pos>,
    /// Distance to feature.
    pub distance: Option<i32>,
    /// Strand of the alignment
    pub strand: i32,
    /// Optional list of warnings and error messages.
    pub messages: Option<Vec<Message>>,
}

/// Response of the `/api/v1/seqvars/csq` endpoint.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::ToSchema)]
pub(crate) struct SeqvarsCsqResponse {
    /// Version information.
    pub version: VersionsInfoResponse,
    /// The original query records.
    pub query: SeqvarsCsqQuery,
    /// The resulting records for the scored genes.
    pub result: Vec<SeqvarsCsqResultEntry>,
}

/// Implementation of endpoints.
async fn handle_impl(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<SeqvarsCsqQuery>,
) -> actix_web::Result<Json<SeqvarsCsqResponse>, super::CustomError> {
    let SeqvarsCsqQuery {
        genome_release,
        chromosome,
        position,
        reference,
        alternative,
        hgnc_id,
    } = query.clone().into_inner();

    let predictor = data
        .seqvars_predictors
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
    let ann_fields = predictor
        .predict(&g_var)
        .map_err(|e| super::CustomError::new(anyhow::anyhow!("prediction failed: {}", &e)))?;
    if let Some(ann_fields) = ann_fields {
        for ann_field in ann_fields {
            if let Some(hgnc_id) = &hgnc_id {
                // Skip if HGNC gene ID is given but does not match.
                if ann_field.gene_id != *hgnc_id {
                    continue;
                }
            }
            let AnnField {
                consequences,
                putative_impact,
                gene_symbol,
                gene_id,
                feature_type,
                feature_id,
                feature_biotype,
                rank,
                hgvs_t,
                hgvs_p,
                tx_pos,
                cds_pos,
                protein_pos,
                distance,
                strand,
                messages,
                ..
            } = ann_field;
            let entry = SeqvarsCsqResultEntry {
                consequences,
                putative_impact,
                gene_symbol,
                gene_id,
                feature_type,
                feature_id,
                feature_biotype: if feature_biotype.contains(&FeatureBiotype::Coding) {
                    FeatureBiotype::Coding
                } else {
                    FeatureBiotype::Noncoding
                },
                feature_tag: feature_biotype
                    .iter()
                    .cloned()
                    .filter(|b| *b != FeatureBiotype::Coding && *b != FeatureBiotype::Noncoding)
                    .collect(),
                rank,
                hgvs_t,
                hgvs_p,
                tx_pos,
                cds_pos,
                protein_pos,
                distance,
                strand,
                messages,
            };
            result.push(entry);
        }
    }

    result.sort_by_key(|e| e.putative_impact);

    let result = SeqvarsCsqResponse {
        version: VersionsInfoResponse::from_web_server_data(data.into_inner().as_ref())
            .map_err(|e| CustomError::new(anyhow::anyhow!("Problem determining version: {}", e)))?,
        query: query.into_inner(),
        result,
    };

    Ok(Json(result))
}

/// Query for consequence of a variant.
#[allow(clippy::unused_async)]
#[get("/seqvars/csq")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<SeqvarsCsqQuery>,
) -> actix_web::Result<Json<SeqvarsCsqResponse>, super::CustomError> {
    handle_impl(data, _path, query).await
}

/// Query for consequence of a variant.
#[allow(clippy::unused_async)]
#[utoipa::path(
    get,
    operation_id = "seqvarsCsq",
    params(
        SeqvarsCsqQuery
    ),
    responses(
        (status = 200, description = "Seqvars consequence information.", body = SeqvarsCsqResponse),
        (status = 500, description = "Internal server error.", body = CustomError)
    )
)]
#[get("/api/v1/seqvars/csq")]
async fn handle_with_openapi(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<SeqvarsCsqQuery>,
) -> actix_web::Result<Json<SeqvarsCsqResponse>, super::CustomError> {
    handle_impl(data, _path, query).await
}
