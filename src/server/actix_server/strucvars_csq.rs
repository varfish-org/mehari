

use actix_web::{
    get,
    web::{self, Data, Json, Path},
    Responder,
};

use crate::{
    annotate::{
        strucvars::csq::{interface, GeneTranscriptEffects},
    },
    common::GenomeRelease,
};

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
    pub variant_type: String,
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
        if self.sv_type() == interface::StrucVarType::Ins {
            self.start
        } else {
            self.stop.unwrap_or(self.start)
        }
    }

    fn sv_type(&self) -> interface::StrucVarType {
        match self.variant_type.to_uppercase().as_ref() {
            "DEL" => interface::StrucVarType::Del,
            "DUP" => interface::StrucVarType::Dup,
            "INS" => interface::StrucVarType::Ins,
            "BND" => interface::StrucVarType::Bnd,
            "INV" => interface::StrucVarType::Inv,
            _ => interface::StrucVarType::Del,
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
    pub result: Vec<GeneTranscriptEffects>,
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
        .strucvar_predictors
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
