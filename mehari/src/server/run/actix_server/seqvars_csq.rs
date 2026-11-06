//! Implementation of endpoint `/api/v1/seqvars/csq`.
//!
//! Also includes the implementation of the `/seqvars/csq` endpoint (deprecated).

use std::collections::BTreeMap;

use actix_web::{
    get,
    web::{self, Data, Json, Path},
};

use super::{CustomError, versions::VersionsInfoResponse};
use crate::annotate::seqvars::consequence::terms::FeatureTag;
use crate::annotate::seqvars::consequence::terms::{
    AnnField, Consequence, FeatureBiotype, FeatureType, Message, Pos, PutativeImpact, Rank,
};
use crate::annotate::seqvars::consequence::{
    SvClaim, SvInfo, Variant, VcfVariant, round_copy_number,
};
use crate::annotate::seqvars::provider::Provider as MehariProvider;
use crate::errors::SeqvarsError;

/// Query parameters of the `/api/v1/seqvars/csq` endpoint.
#[derive(
    Debug, Clone, serde::Serialize, serde::Deserialize, utoipa::IntoParams, utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
#[serde_with::skip_serializing_none]
pub(crate) struct SeqvarsCsqQuery {
    /// The assembly.
    pub assembly: String,
    /// SPDI sequence.
    pub chromosome: String,
    /// SPDI position.
    pub position: u32,
    /// SPDI deletion.
    pub reference: String,
    /// SPDI insertion; also a symbolic allele such as `<DEL>`, or a breakend such as
    /// `G]13:123456]`.
    pub alternative: String,
    /// Optionally, the HGNC ID of the gene to limit to.
    pub hgnc_id: Option<String>,
    /// For a symbolic allele, the position of the last changed base (`END`).
    pub end: Option<i32>,
    /// For a symbolic allele, its length (`SVLEN`); for `<INS>`, the number of inserted bases.
    pub sv_length: Option<i32>,
    /// For a symbolic allele, what the call claims (`SVCLAIM`): `D`, `J` or `DJ`.
    pub sv_claim: Option<SvClaim>,
    /// For a symbolic allele or a breakend, whether the breakpoints are imprecise (`IMPRECISE`).
    pub imprecise: Option<bool>,
    /// For a custom symbolic allele such as `<FOO>`, the type it counts as (`SVTYPE`): `DEL`,
    /// `INS`, `DUP`, `INV`, `CNV` or `BND`.
    pub sv_type: Option<String>,
    /// For a symbolic allele, the copy number of the allele (`CN` as of VCF 4.4), rounded: for a
    /// `<CNV>`, 0 is a deletion, 1 no variant and 2 or more a duplication.
    pub cn: Option<f64>,
    /// For a symbolic allele or a breakend, the confidence interval around `position` (`CIPOS`)
    /// as `min,max`, e.g. `-500,500`.
    pub cipos: Option<String>,
    /// For a symbolic allele, the confidence interval around `end` (`CIEND`) as `min,max`;
    /// `cipos` if unset.
    pub ciend: Option<String>,
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
    /// The gene ID of the annotation, e.g., `"672"` for BRCA1 in RefSeq.
    pub gene_id: String,
    /// The HGNC ID of the gene, e.g., `"HGNC:1100"` for BRCA1; unset for a gene without one.
    pub hgnc_id: Option<String>,
    /// The feature type.
    pub feature_type: FeatureType,
    /// The feature identifier.
    pub feature_id: String,
    /// The feature biotype.
    pub feature_biotype: FeatureBiotype,
    /// The feature tags.
    pub feature_tags: Vec<FeatureTag>,
    /// The exon / intron rank.
    pub rank: Option<Rank>,
    /// HGVS g. notation.
    pub hgvs_g: Option<String>,
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
    /// Values of the custom `ANN` columns, e.g. from plugins, by column, as written into `ANN`.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub custom_fields: BTreeMap<String, String>,
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
    /// The gene fusions of a structural variant or a breakend as `GENE1::GENE2`, as in the
    /// `FUSION` INFO field of `annotate seqvars`.
    pub fusions: Vec<String>,
}

/// Implementation of endpoints.
async fn handle_impl(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<SeqvarsCsqQuery>,
) -> actix_web::Result<Json<SeqvarsCsqResponse>, super::CustomError> {
    let SeqvarsCsqQuery {
        assembly: genome_release,
        chromosome,
        position,
        reference,
        alternative,
        hgnc_id,
        end,
        sv_length,
        sv_claim,
        imprecise,
        sv_type,
        cn,
        cipos,
        ciend,
    } = query.clone().into_inner();

    let predictor = data
        .seqvars_predictors
        .get(&genome_release)
        .ok_or_else(|| {
            super::CustomError::new(anyhow::anyhow!(
                "genome release not supported: {:?}",
                query.assembly
            ))
        })?;

    let mut result = Vec::new();
    let position = i32::try_from(position).map_err(|_| {
        CustomError::bad_request(anyhow::anyhow!("position {position} is too large"))
    })?;
    let mut fusions = Vec::new();
    let g_var = VcfVariant {
        chromosome,
        position,
        reference,
        alternative,
    };
    let info = SvInfo {
        end,
        sv_length,
        sv_type: sv_type.and_then(|sv_type| sv_type.parse().ok()),
        claim: sv_claim,
        copy_number: (cn.map(round_copy_number).transpose())
            .map_err(|e| CustomError::bad_request(e.into()))?,
        cipos: interval("cipos", cipos.as_deref())?,
        ciend: interval("ciend", ciend.as_deref())?,
        imprecise: imprecise.unwrap_or_default(),
    };
    let annotations = Variant::from_info(g_var, &info)
        .and_then(|variant| variant.map_or(Ok(None), |variant| predictor.predict_variant(&variant)))
        .map_err(|e| match e {
            SeqvarsError::InvalidVariant(_)
            | SeqvarsError::InvalidCoordinates(_, _)
            | SeqvarsError::UnknownChromosomeAccession => super::CustomError::bad_request(e.into()),
            _ => super::CustomError::new(anyhow::anyhow!("prediction failed: {}", e)),
        })?;
    if let Some(annotations) = annotations {
        fusions = annotations.fusions;
        for ann_field in annotations.anns {
            if let Some(hgnc_id) = &hgnc_id {
                // Skip if HGNC gene ID is given but does not match.
                if !has_hgnc_id(&predictor.provider, &ann_field.feature_id, hgnc_id) {
                    continue;
                }
            }
            let AnnField {
                consequences,
                putative_impact,
                gene_symbol,
                gene_id,
                hgnc_id: gene_hgnc_id,
                feature_type,
                feature_id,
                feature_biotype,
                feature_tags,
                rank,
                hgvs_g,
                hgvs_c: hgvs_t,
                hgvs_p,
                cdna_pos: tx_pos,
                cds_pos,
                protein_pos,
                distance,
                strand,
                messages,
                custom_fields,
                ..
            } = ann_field;
            let entry = SeqvarsCsqResultEntry {
                consequences,
                putative_impact,
                gene_symbol,
                gene_id,
                hgnc_id: gene_hgnc_id,
                feature_type,
                feature_id,
                feature_biotype: if feature_biotype.contains(&FeatureBiotype::Coding) {
                    FeatureBiotype::Coding
                } else {
                    FeatureBiotype::Noncoding
                },
                feature_tags,
                rank,
                hgvs_g,
                hgvs_t,
                hgvs_p,
                tx_pos,
                cds_pos,
                protein_pos,
                distance,
                strand,
                messages,
                custom_fields: custom_fields
                    .into_iter()
                    .filter_map(|(column, value)| Some((column, value?)))
                    .collect(),
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
        fusions,
    };

    Ok(Json(result))
}

/// The confidence interval `value` of the query parameter `name` as `min,max`, e.g. `-500,500`.
fn interval(name: &str, value: Option<&str>) -> Result<Option<(i32, i32)>, CustomError> {
    value
        .map(|value| {
            (value.split_once(','))
                .and_then(|(min, max)| Some((min.parse().ok()?, max.parse().ok()?)))
                .ok_or_else(|| {
                    CustomError::bad_request(anyhow::anyhow!(
                        "{name} {value:?} is not two integers min,max"
                    ))
                })
        })
        .transpose()
}

/// Whether the gene of the transcript `tx_ac` has the HGNC ID `hgnc_id`.
fn has_hgnc_id(provider: &MehariProvider, tx_ac: &str, hgnc_id: &str) -> bool {
    provider.get_tx(tx_ac).and_then(|tx| tx.hgnc_id.as_deref()) == Some(hgnc_id)
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
        (status = 400, description = "Invalid variant, e.g. a symbolic allele without its end.", body = CustomError),
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

#[cfg(test)]
mod test {
    /// The filter by HGNC ID uses `hgnc_id`, not `gene_id`.
    #[rstest::rstest]
    #[case::hgnc_id("HGNC:1", true)]
    #[case::gene_id("11", false)]
    fn has_hgnc_id(#[case] hgnc_id: &str, #[case] expected: bool) {
        use crate::annotate::seqvars::provider::{Provider, test_db};

        let provider = Provider::new(test_db::one_gene_db(), None::<&str>, true);

        assert_eq!(
            super::has_hgnc_id(&provider, "NM_000001.1", hgnc_id),
            expected
        );
    }

    /// The response for `alternative` at `position` on chr17 of GRCh37, with `reference` and the
    /// further query parameters `params`, e.g. `end=41197760&cipos=-500,500`.
    async fn csq(
        position: u32,
        reference: &str,
        alternative: &str,
        params: &str,
    ) -> Result<Result<super::SeqvarsCsqResponse, super::CustomError>, anyhow::Error> {
        use crate::annotate::seqvars::consequence::logic::ConsequencePredictor;
        use crate::annotate::seqvars::consequence::{ConfigBuilder, load_tx_db};
        use crate::annotate::seqvars::provider::Provider;
        use crate::server::run::actix_server::WebServerData;
        use actix_web::web::{Data, Path, Query};
        use std::sync::Arc;

        let tx_db = load_tx_db("tests/data/annotate/db/grch37/txs.bin.zst")?;
        let provider = Arc::new(Provider::new(tx_db, None::<&str>, true));
        let predictor =
            ConsequencePredictor::new(provider.clone(), ConfigBuilder::default().build()?);
        let data = WebServerData {
            provider: [("grch37".to_string(), provider)].into(),
            seqvars_predictors: [("grch37".to_string(), predictor)].into(),
            ..Default::default()
        };
        let query = Query::<super::SeqvarsCsqQuery>::from_query(&format!(
            "assembly=grch37&chromosome=17&position={position}&reference={reference}\
             &alternative={alternative}&{params}"
        ))?;
        let response = super::handle_impl(Data::new(data), Path::from(()), query).await;
        Ok(response.map(|response| response.into_inner()))
    }

    /// A symbolic allele with its end or length, and a breakend, are annotated as in `annotate
    /// seqvars`; all lie in BRCA1. A custom allele counts as its `sv_type`, which other alleles
    /// ignore. A `<CNV>` with 0 copies (`cn`, rounded) and a call with `sv_claim=D` are imprecise
    /// deletions, and `imprecise` makes the breakpoints imprecise, within `cipos` and `ciend`
    /// where given, also of a breakend.
    #[rstest::rstest]
    #[case::deletion(
        41197710,
        "N",
        "<DEL>",
        "end=41197760",
        "g.41197711_41197760del",
        false
    )]
    #[case::custom_allele(
        41197710,
        "N",
        "<FOO>",
        "end=41197760&sv_type=DEL",
        "g.41197711_41197760del",
        false
    )]
    #[case::copy_number_loss(
        41197710,
        "N",
        "<CNV>",
        "end=41197760&cn=0.4",
        "g.(?_41197711)_(41197760_?)del",
        true
    )]
    #[case::imprecise_deletion(
        41197710,
        "N",
        "<DEL>",
        "end=41197760&imprecise=true&cipos=-5,5",
        "g.(41197706_41197716)_(41197755_41197765)del",
        true
    )]
    #[case::imprecise_without_intervals(
        41197710,
        "N",
        "<DEL>",
        "end=41197760&imprecise=true",
        "g.(?_41197711)_(41197760_?)del",
        true
    )]
    #[case::end_interval(
        41197710,
        "N",
        "<DEL>",
        "end=41197760&imprecise=true&cipos=-5,5&ciend=-10,10",
        "g.(41197706_41197716)_(41197750_41197770)del",
        true
    )]
    #[case::length(
        41197710,
        "N",
        "<DEL>",
        "sv_length=50",
        "g.41197711_41197760del",
        false
    )]
    #[case::copy_number_call(
        41197710,
        "N",
        "<DEL>",
        "end=41197760&sv_claim=D",
        "g.(?_41197711)_(41197760_?)del",
        true
    )]
    #[case::type_of_standard_allele(
        41197710,
        "N",
        "<DEL>",
        "end=41197760&sv_type=CPX",
        "g.41197711_41197760del",
        false
    )]
    #[case::breakend(
        41215930,
        "G",
        "G]13:32900000]",
        "",
        "g.41215931_qterdelins[NC_000013.10:g.pter_32900000inv]",
        false
    )]
    #[case::imprecise_breakend(
        41215930,
        "G",
        "G]13:32900000]",
        "imprecise=true&cipos=-500,500",
        "g.(41215431_41216431)_qterdelins[NC_000013.10:g.pter_32900000inv]",
        true
    )]
    #[case::type_of_breakend(
        41215930,
        "G",
        "G]13:32900000]",
        "sv_type=TRA",
        "g.41215931_qterdelins[NC_000013.10:g.pter_32900000inv]",
        false
    )]
    #[tokio::test]
    async fn structural_variant(
        #[case] position: u32,
        #[case] reference: &str,
        #[case] alternative: &str,
        #[case] params: &str,
        #[case] hgvs_g: &str,
        #[case] imprecise: bool,
    ) -> Result<(), anyhow::Error> {
        use crate::annotate::seqvars::consequence::terms::Message;

        let response = csq(position, reference, alternative, params)
            .await?
            .map_err(|e| anyhow::anyhow!("{e}"))?;
        let brca1 = (response.result.iter())
            .find(|entry| entry.feature_id == "NM_007294.4")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_007294.4"))?;
        assert_eq!(brca1.hgvs_g.as_deref(), Some(hgvs_g));
        let messages = brca1.messages.as_deref().unwrap_or_default();
        assert_eq!(
            messages.contains(&Message::InfoImpreciseBreakpoints),
            imprecise
        );
        Ok(())
    }

    /// A `<CNV>` that keeps the copy number of the reference, `cn=1`, is no variant, as in
    /// `annotate seqvars`.
    #[tokio::test]
    async fn copy_number_kept() -> Result<(), anyhow::Error> {
        let response = csq(41197710, "N", "<CNV>", "end=41197760&cn=1")
            .await?
            .map_err(|e| anyhow::anyhow!("{e}"))?;
        assert!(response.result.is_empty());
        Ok(())
    }

    /// The response lists the gene fusions, here of the breakend that joins the 3' part of BRCA1
    /// on chr17 to the reversed 5' part of OPA1 on chr3, which reads into it; the rows of BRCA1
    /// get the fusion terms.
    #[tokio::test]
    async fn fusions() -> Result<(), anyhow::Error> {
        use crate::annotate::seqvars::consequence::terms::Consequence;

        let response = csq(41215930, "G", "G]3:193360000]", "")
            .await?
            .map_err(|e| anyhow::anyhow!("{e}"))?;
        assert_eq!(response.fusions, ["OPA1::BRCA1"]);
        let brca1 = (response.result.iter())
            .find(|entry| entry.feature_id == "NM_007294.4")
            .ok_or_else(|| anyhow::anyhow!("no annotation for NM_007294.4"))?;
        assert!(brca1.consequences.contains(&Consequence::GeneFusion));
        Ok(())
    }

    /// `*`, `<*>`, `<NON_REF>` and `ALT` `.` are no variants of their own, so they have no
    /// annotations, as in `annotate seqvars`.
    #[rstest::rstest]
    #[case::spanning_deletion("*")]
    #[case::any_allele("<*>")]
    #[case::non_ref("<NON_REF>")]
    #[case::no_allele(".")]
    #[tokio::test]
    async fn no_variant_of_its_own(#[case] alternative: &str) -> Result<(), anyhow::Error> {
        let response = csq(41197701, "G", alternative, "")
            .await?
            .map_err(|e| anyhow::anyhow!("{e}"))?;
        assert!(response.result.is_empty());
        Ok(())
    }

    /// A symbolic allele that cannot be annotated, a malformed breakend, an unknown type, a
    /// negative copy number, a malformed confidence interval and a position past `i32::MAX` are a
    /// bad request, with the reason as JSON.
    #[rstest::rstest]
    #[case::deletion_without_end(41197710, "<DEL>", "", "has no end past its position")]
    #[case::end_before_position(41197710, "<DEL>", "end=41197700", "has no end past its position")]
    #[case::unknown_symbolic_allele(41197710, "<FOO>", "end=41197760", "is not supported")]
    #[case::malformed_breakend(41197710, "G]13:x]", "", "invalid mate position")]
    #[case::unknown_type(41197710, "<FOO>", "end=41197760&sv_type=FOO", "SVTYPE")]
    #[case::negative_copy_number(41197710, "<CNV>", "end=41197760&cn=-1", "CN")]
    #[case::malformed_interval(41197710, "<DEL>", "end=41197760&cipos=-500", "cipos")]
    #[case::reversed_interval(41197710, "<DEL>", "end=41197760&ciend=500,-500", "CIEND")]
    #[case::position_too_large(1 << 31, "C", "", "position")]
    #[tokio::test]
    async fn bad_request(
        #[case] position: u32,
        #[case] alternative: &str,
        #[case] params: &str,
        #[case] message: &str,
    ) -> Result<(), anyhow::Error> {
        use actix_web::ResponseError as _;

        let result = csq(position, "N", alternative, params).await?;

        let error = result.err().ok_or_else(|| anyhow::anyhow!("no error"))?;
        let response = error.error_response();
        assert_eq!(response.status(), actix_web::http::StatusCode::BAD_REQUEST);
        assert_eq!(
            response
                .headers()
                .get(actix_web::http::header::CONTENT_TYPE),
            Some(&actix_web::http::header::HeaderValue::from_static(
                "application/json"
            ))
        );
        let body = actix_web::body::to_bytes(response.into_body())
            .await
            .map_err(|e| anyhow::anyhow!("{e}"))?;
        let body: serde_json::Value = serde_json::from_slice(&body)?;
        assert!(
            body["err"]
                .as_str()
                .is_some_and(|err| err.contains(message)),
            "{body}"
        );
        Ok(())
    }

    /// An `sv_claim` other than `D`, `J` and `DJ` fails to deserialize, a bad request.
    #[test]
    fn unknown_claim() -> Result<(), anyhow::Error> {
        use actix_web::ResponseError as _;

        let error = actix_web::web::Query::<super::SeqvarsCsqQuery>::from_query(
            "assembly=grch37&chromosome=17&position=41197710&reference=N&alternative=<DEL>\
             &end=41197760&sv_claim=X",
        )
        .err()
        .ok_or_else(|| anyhow::anyhow!("no error"))?;
        assert_eq!(
            error.status_code(),
            actix_web::http::StatusCode::BAD_REQUEST
        );
        assert!(error.to_string().contains("unknown variant `X`"), "{error}");
        Ok(())
    }
}
