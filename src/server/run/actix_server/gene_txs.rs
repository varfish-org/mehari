//! Implementation of `/seqvars/csq` endpoint.

use crate::common::GenomeRelease;
use crate::pbs;
use crate::pbs::server::{GeneTranscriptsQuery, GeneTranscriptsResponse};
use crate::pbs::txs::GenomeBuild;
use crate::server::run::actix_server::CustomError;
use actix_web::{
    get,
    web::{self, Data, Json, Path},
};
use hgvs::data::interface::Provider as _;

use super::versions::Assembly;

/// Maximal page size.
static PAGE_SIZE_MAX: i32 = 1000;
/// Default page size.
static PAGE_SIZE_DEFAULT: i32 = 100;

/// Core implementation of the `/genes/txs` and `/api/v1/genes/transcripts` endpoints.
///
/// For now takes the `GeneTranscriptsQuery` as the argument and returns
/// the `GeneTranscriptsResponse` as the result.
fn genes_tx_impl(
    data: Data<super::WebServerData>,
    query: GeneTranscriptsQuery,
) -> Result<GeneTranscriptsResponse, CustomError> {
    let GeneTranscriptsQuery {
        genome_build,
        hgnc_id,
        page_size,
        next_page_token,
    } = query;
    let genome_build = GenomeBuild::try_from(genome_build.unwrap_or(GenomeBuild::Grch37 as i32))
        .map_err(|e| CustomError::new(anyhow::anyhow!("Invalid genome build: {}", e)))?;
    let genome_release = GenomeRelease::try_from(genome_build)
        .map_err(|e| CustomError::new(anyhow::anyhow!("Invalid genome build: {}", e)))?;
    let hgnc_id = hgnc_id
        .as_ref()
        .ok_or_else(|| CustomError::new(anyhow::anyhow!("No HGNC ID provided.")))?;
    let page_size = page_size
        .unwrap_or(PAGE_SIZE_DEFAULT)
        .min(PAGE_SIZE_MAX)
        .max(1);

    let provider = data
        .provider
        .get(&genome_release)
        .ok_or_else(|| CustomError::new(anyhow::anyhow!("No provider available.")))?;
    let tx_acs = provider
        .get_tx_for_gene(hgnc_id)
        .map_err(|e| CustomError::new(anyhow::anyhow!("No transcripts found: {}", e)))?
        .into_iter()
        .map(|tx| tx.tx_ac)
        .collect::<Vec<_>>();
    let first = next_page_token
        .as_ref()
        .and_then(|next_page_token| tx_acs.iter().position(|tx_ac| tx_ac == next_page_token))
        .unwrap_or(0);
    let last = (first + page_size as usize).min(tx_acs.len());

    Ok(GeneTranscriptsResponse {
        transcripts: tx_acs[first..last]
            .iter()
            .filter_map(|tx_ac| provider.get_tx(tx_ac))
            .collect::<Vec<_>>(),
        next_page_token: if last < tx_acs.len() {
            Some(tx_acs[last].clone())
        } else {
            None
        },
    })
}

/// Implementation of the `/genes/txs` endpoint.
#[allow(clippy::unused_async)]
#[get("/genes/txs")]
async fn handle(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<GeneTranscriptsQuery>,
) -> actix_web::Result<Json<GeneTranscriptsResponse>, CustomError> {
    Ok(Json(genes_tx_impl(data, query.into_inner())?))
}

/// Query arguments for the `/api/v1/genes/transcripts` endpoint.
#[derive(Debug, Clone, serde::Deserialize, utoipa::ToSchema)]
struct GenesTranscriptsListQuery {
    /// HGNC gene ID.
    pub hgnc_id: String,
    /// Genome build.
    pub genome_build: Assembly,
    /// Page size.
    pub page_size: Option<i32>,
    /// Next page token.
    pub next_page_token: Option<String>,
}

impl From<GenesTranscriptsListQuery> for pbs::server::GeneTranscriptsQuery {
    fn from(val: GenesTranscriptsListQuery) -> Self {
        pbs::server::GeneTranscriptsQuery {
            genome_build: Some(Into::<pbs::txs::GenomeBuild>::into(val.genome_build) as i32),
            hgnc_id: Some(val.hgnc_id),
            page_size: val.page_size,
            next_page_token: val.next_page_token,
        }
    }
}

/// Enumeration for `Transcript::biotype`.
#[derive(Debug, Clone, serde::Serialize, utoipa::ToSchema)]
#[serde(rename_all = "snake_case")]
enum TranscriptBiotype {
    /// Coding transcript.
    Coding,
    /// Non-coding transcript.
    NonCoding,
}

impl TryFrom<pbs::txs::TranscriptBiotype> for TranscriptBiotype {
    type Error = anyhow::Error;

    fn try_from(value: pbs::txs::TranscriptBiotype) -> Result<Self, Self::Error> {
        match value {
            pbs::txs::TranscriptBiotype::Coding => Ok(TranscriptBiotype::Coding),
            pbs::txs::TranscriptBiotype::NonCoding => Ok(TranscriptBiotype::NonCoding),
            _ => Err(anyhow::anyhow!("Invalid biotype: {:?}", value)),
        }
    }
}

// Bit values for the transcript tags.
#[derive(Debug, Clone, serde::Serialize, utoipa::ToSchema)]
#[serde(rename_all = "snake_case")]
enum TranscriptTag {
    /// Member of Ensembl basic.
    Basic,
    /// Member of Ensembl canonical.
    EnsemblCanonical,
    /// Member of MANE Select.
    ManeSelect,
    /// Member of MANE Plus Clinical.
    ManePlusClinical,
    /// Member of RefSeq Select.
    RefSeqSelect,
    /// Flagged as being a selenoprotein (UGA => selenon).
    Selenoprotein,
    /// Member of GENCODE Primary
    GencodePrimary,
    /// Catchall for other tags.
    Other,
}

impl TryFrom<pbs::txs::TranscriptTag> for TranscriptTag {
    type Error = anyhow::Error;

    fn try_from(value: pbs::txs::TranscriptTag) -> Result<Self, Self::Error> {
        match value {
            pbs::txs::TranscriptTag::Basic => Ok(TranscriptTag::Basic),
            pbs::txs::TranscriptTag::EnsemblCanonical => Ok(TranscriptTag::EnsemblCanonical),
            pbs::txs::TranscriptTag::ManeSelect => Ok(TranscriptTag::ManeSelect),
            pbs::txs::TranscriptTag::ManePlusClinical => Ok(TranscriptTag::ManePlusClinical),
            pbs::txs::TranscriptTag::RefSeqSelect => Ok(TranscriptTag::RefSeqSelect),
            pbs::txs::TranscriptTag::Selenoprotein => Ok(TranscriptTag::Selenoprotein),
            pbs::txs::TranscriptTag::GencodePrimary => Ok(TranscriptTag::GencodePrimary),
            pbs::txs::TranscriptTag::Other => Ok(TranscriptTag::Other),
            _ => Err(anyhow::anyhow!("Invalid transcript tag: {:?}", value)),
        }
    }
}

/// Enumeration for the two strands of the genome.
#[derive(Debug, Clone, serde::Serialize, utoipa::ToSchema)]
#[serde(rename_all = "snake_case")]
enum Strand {
    /// unknown
    Unknown,
    /// Forward / plus
    Plus,
    /// Reverse / minus
    Minus,
}

impl From<pbs::txs::Strand> for Strand {
    fn from(value: pbs::txs::Strand) -> Self {
        match value {
            pbs::txs::Strand::Unknown => Strand::Unknown,
            pbs::txs::Strand::Plus => Strand::Plus,
            pbs::txs::Strand::Minus => Strand::Minus,
        }
    }
}

/// Store the alignment of one exon to the reference.
#[derive(Debug, Clone, serde::Serialize, utoipa::ToSchema)]
struct ExonAlignment {
    /// Start position on reference.
    pub alt_start_i: i32,
    /// End position on reference.
    pub alt_end_i: i32,
    /// Exon number.
    pub ord: i32,
    /// CDS start coordinate.
    pub alt_cds_start_i: Option<i32>,
    /// CDS end coordinate.
    pub alt_cds_end_i: Option<i32>,
    /// CIGAR string of alignment, empty indicates full matches.
    pub cigar: String,
}

impl From<pbs::txs::ExonAlignment> for ExonAlignment {
    fn from(value: pbs::txs::ExonAlignment) -> Self {
        ExonAlignment {
            alt_start_i: value.alt_start_i,
            alt_end_i: value.alt_end_i,
            ord: value.ord,
            alt_cds_start_i: value.alt_cds_start_i,
            alt_cds_end_i: value.alt_cds_end_i,
            cigar: value.cigar.clone(),
        }
    }
}

/// Store information about a transcript aligning to a genome.
#[derive(Debug, Clone, serde::Serialize, utoipa::ToSchema)]
struct GenomeAlignment {
    /// The genome build identifier.
    pub genome_build: Assembly,
    /// Accession of the contig sequence.
    pub contig: String,
    /// CDS end position, `-1` to indicate `None`.
    pub cds_start: Option<i32>,
    /// CDS end position, `-1` to indicate `None`.
    pub cds_end: Option<i32>,
    /// The strand.
    pub strand: Strand,
    /// Exons of the alignment.
    pub exons: Vec<ExonAlignment>,
}

impl TryFrom<pbs::txs::GenomeAlignment> for GenomeAlignment {
    type Error = anyhow::Error;

    fn try_from(value: pbs::txs::GenomeAlignment) -> Result<Self, Self::Error> {
        Ok(GenomeAlignment {
            genome_build: Assembly::try_from(pbs::txs::GenomeBuild::try_from(value.genome_build)?)?,
            contig: value.contig.clone(),
            cds_start: value.cds_start,
            cds_end: value.cds_end,
            strand: Strand::from(pbs::txs::Strand::try_from(value.strand)?),
            exons: value.exons.into_iter().map(Into::into).collect(),
        })
    }
}

/// Transcript information.
#[derive(Debug, Clone, serde::Serialize, utoipa::ToSchema)]
struct Transcript {
    /// Transcript accession with version, e.g., `"NM_007294.3"` or `"ENST00000461574.1"` for BRCA1.
    pub id: String,
    /// HGNC symbol, e.g., `"BRCA1"`
    pub gene_symbol: String,
    /// HGNC gene identifier, e.g., `"1100"` for BRCA1.
    pub gene_id: String,
    /// Transcript biotype.
    pub biotype: TranscriptBiotype,
    /// Transcript flags.
    pub tags: Vec<TranscriptTag>,
    /// Identifier of the corresponding protein.
    pub protein: Option<String>,
    /// CDS start codon.
    pub start_codon: Option<i32>,
    /// CDS stop codon.
    pub stop_codon: Option<i32>,
    /// Alignments on the different genome builds.
    pub genome_alignments: Vec<GenomeAlignment>,
    /// Whether this transcript has an issue (e.g. MissingStopCodon), cf. `mehari::db::create::mod::Reason`.
    pub filtered: Option<bool>,
    /// Reason for filtering.
    pub filter_reason: ::core::option::Option<u32>,
}

impl TryFrom<pbs::txs::Transcript> for Transcript {
    type Error = anyhow::Error;

    fn try_from(value: pbs::txs::Transcript) -> Result<Self, Self::Error> {
        Ok(Transcript {
            id: value.id.clone(),
            gene_symbol: value.gene_symbol.clone(),
            gene_id: value.gene_id.clone(),
            biotype: TranscriptBiotype::try_from(pbs::txs::TranscriptBiotype::try_from(
                value.biotype,
            )?)?,
            tags: value
                .tags
                .into_iter()
                .map(|i32_tag| -> Result<_, anyhow::Error> {
                    TranscriptTag::try_from(pbs::txs::TranscriptTag::try_from(i32_tag)?)
                })
                .collect::<Result<Vec<_>, _>>()?,
            protein: value.protein.clone(),
            start_codon: value.start_codon,
            stop_codon: value.stop_codon,
            genome_alignments: value
                .genome_alignments
                .into_iter()
                .map(TryInto::try_into)
                .collect::<Result<_, _>>()?,
            filtered: value.filtered,
            filter_reason: value.filter_reason,
        })
    }
}

/// Response of the `/api/v1/genes/transcripts` endpoint.
#[derive(Debug, Clone, serde::Serialize, utoipa::ToSchema)]
struct GenesTranscriptsListResponse {
    /// The transcripts for the gene.
    pub transcripts: Vec<Transcript>,
    /// The token to continue from a previous query.
    pub next_page_token: Option<String>,
}

impl TryFrom<pbs::server::GeneTranscriptsResponse> for GenesTranscriptsListResponse {
    type Error = anyhow::Error;

    fn try_from(value: pbs::server::GeneTranscriptsResponse) -> Result<Self, Self::Error> {
        Ok(GenesTranscriptsListResponse {
            transcripts: value
                .transcripts
                .into_iter()
                .map(TryInto::try_into)
                .collect::<Result<_, _>>()?,
            next_page_token: value.next_page_token,
        })
    }
}

/// Query for consequence of a variant.
#[allow(clippy::unused_async)]
#[utoipa::path(
    get,
    operation_id = "genesTranscriptsList",
    responses(
        (status = 200, description = "Transcripts for the selected gene.", body = GenesTranscriptsListResponse),
        (status = 500, description = "Internal server error.", body = CustomError)
    )
)]
#[get("/api/v1/genes/transcripts")]
async fn handle_with_openapi(
    data: Data<super::WebServerData>,
    _path: Path<()>,
    query: web::Query<GenesTranscriptsListQuery>,
) -> actix_web::Result<Json<GenesTranscriptsListResponse>, CustomError> {
    let result = genes_tx_impl(data, query.into_inner().into())?;
    Ok(Json(result.try_into().map_err(|e| {
        CustomError::new(anyhow::anyhow!("Conversion error: {}", e))
    })?))
}
