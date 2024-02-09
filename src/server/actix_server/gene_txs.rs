//! Implementation of `/seqvars/csq` endpoint.

use crate::common::GenomeRelease;
use crate::pbs::server::{GeneTranscriptsQuery, GeneTranscriptsResponse};
use crate::pbs::txs::GenomeBuild;
use actix_web::{
    get,
    web::{self, Data, Json, Path},
    Responder,
};
use hgvs::data::interface::Provider as _;

/// Maximal page size.
static PAGE_SIZE_MAX: i32 = 1000;
/// Default page size.
static PAGE_SIZE_DEFAULT: i32 = 100;

/// Query for consequence of a variant.
#[allow(clippy::unused_async)]
#[get("/genes/txs")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<GeneTranscriptsQuery>,
) -> actix_web::Result<impl Responder, super::CustomError> {
    let GeneTranscriptsQuery {
        genome_build,
        hgnc_id,
        page_size,
        next_page_token,
    } = query.clone().into_inner();
    let genome_build = GenomeBuild::try_from(genome_build.unwrap_or(GenomeBuild::Grch37 as i32))
        .map_err(|e| super::CustomError::new(anyhow::anyhow!("Invalid genome build: {}", e)))?;
    let genome_release = GenomeRelease::try_from(genome_build)
        .map_err(|e| super::CustomError::new(anyhow::anyhow!("Invalid genome build: {}", e)))?;
    let hgnc_id = hgnc_id
        .as_ref()
        .ok_or_else(|| super::CustomError::new(anyhow::anyhow!("No HGNC ID provided.")))?;
    let page_size = page_size
        .unwrap_or(PAGE_SIZE_DEFAULT)
        .min(PAGE_SIZE_MAX)
        .max(1);

    let provider = data
        .provider
        .get(&genome_release)
        .ok_or_else(|| super::CustomError::new(anyhow::anyhow!("No provider available.")))?;
    let tx_acs = provider
        .get_tx_for_gene(hgnc_id)
        .map_err(|e| super::CustomError::new(anyhow::anyhow!("No transcripts found: {}", e)))?
        .into_iter()
        .map(|tx| tx.tx_ac)
        .collect::<Vec<_>>();
    let first = next_page_token
        .as_ref()
        .and_then(|next_page_token| tx_acs.iter().position(|tx_ac| tx_ac == next_page_token))
        .unwrap_or(0);
    let last = (first + page_size as usize).min(tx_acs.len());

    Ok(Json(GeneTranscriptsResponse {
        transcripts: tx_acs[first..last]
            .iter()
            .filter_map(|tx_ac| provider.get_tx(tx_ac))
            .collect::<Vec<_>>(),
        next_page_token: if last < tx_acs.len() {
            Some(tx_acs[last].clone())
        } else {
            None
        },
    }))
}
