//! Implementation of `/seqvars/csq` endpoint.

use actix_web::{
    get,
    web::{self, Data, Json, Path},
    Responder,
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

/// Parameters for `/seqvars/csq`.
///
#[derive(serde::Serialize, serde::Deserialize, Debug, Clone)]
#[serde(rename_all = "snake_case")]
#[serde_with::skip_serializing_none]
struct Query {
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

/// Result entry for the API.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
struct ResultEntry {
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
    /// Optional list of warnings and error messages.
    pub messages: Option<Vec<Message>>,
}

/// Container for the result.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
struct Container {
    /// Version information.
    pub version: crate::common::Version,
    /// The original query records.
    pub query: Query,
    /// The resulting records for the scored genes.
    pub result: Vec<ResultEntry>,
}

/// Query for consequence of a variant.
#[allow(clippy::unused_async)]
#[get("/seqvars/csq")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<Query>,
) -> actix_web::Result<impl Responder, super::CustomError> {
    let Query {
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
        reference,
        alternative,
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
                messages,
                ..
            } = ann_field;
            let entry = ResultEntry {
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
                messages,
            };
            result.push(entry);
        }
    }

    result.sort_by_key(|e| e.putative_impact);

    let result = Container {
        version: crate::common::Version::new(predictor.data_version()),
        query: query.into_inner(),
        result,
    };

    Ok(Json(result))
}
