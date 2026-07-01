use crate::annotate::seqvars::consequence::VcfVariant;
use crate::common::contig::ContigManager;
use crate::db::keys;
use anyhow::Error;
use prost::Message;
use rocksdb::DBWithThreadMode;
use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::io::Cursor;
use std::path::Path;
use std::sync::Arc;

#[derive(Debug)]
pub struct ClinvarAnnotator {
    db: DBWithThreadMode<rocksdb::MultiThreaded>,
    contig_manager: Arc<ContigManager>,
}

impl ClinvarAnnotator {
    pub fn new(
        db: DBWithThreadMode<rocksdb::MultiThreaded>,
        contig_manager: Arc<ContigManager>,
    ) -> Self {
        Self { db, contig_manager }
    }

    pub(crate) fn from_path(
        path: impl AsRef<Path> + Display,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        tracing::info!("Opening ClinVar database");
        tracing::debug!("RocksDB path = {}", &path);
        let options = rocksdb::Options::default();
        let db_clinvar =
            rocksdb::DB::open_cf_for_read_only(&options, &path, ["meta", "clinvar"], false)?;
        Ok(Self::new(db_clinvar, contig_manager))
    }

    /// Annotate record with ClinVar information
    pub fn annotate_record_clinvar(&self, key: &[u8]) -> Result<Option<ClinvarResult>, Error> {
        if let Some(raw_value) = self
            .db
            .get_cf(self.db.cf_handle("clinvar").as_ref().unwrap(), key)?
        {
            let record_list = annonars::pbs::clinvar::minimal::ExtractedVcvRecordList::decode(
                &mut Cursor::new(&raw_value),
            )?;

            let mut vcv = Vec::new();
            let mut germline_classification = Vec::new();
            for r in record_list.records.iter() {
                let accession = r.accession.as_ref().expect("must have VCV");
                vcv.push(format!("{}.{}", accession.accession, accession.version));
                if let Some(gc) = &r
                    .classifications
                    .as_ref()
                    .expect("has cls")
                    .germline_classification
                {
                    germline_classification
                        .push(gc.description.as_ref().expect("has desc").to_string());
                }
            }
            Ok(Some(ClinvarResult {
                vcv,
                germline_classification,
            }))
        } else {
            Ok(None)
        }
    }

    pub(crate) fn annotate(&self, vcf_var: &keys::Var) -> anyhow::Result<Option<ClinvarResult>> {
        // Only attempt lookups into RocksDB for canonical contigs.
        if self
            .contig_manager
            .is_homo_sapiens_canonical_alias(vcf_var.chrom.as_str())
        {
            // Build key for RocksDB database from `vcf_var`.
            let annonars_var: annonars::common::keys::Var = vcf_var.into();
            let key: Vec<u8> = annonars_var.into();
            return self.annotate_record_clinvar(&key);
        }
        Ok(None)
    }

    #[cfg(feature = "server")]
    pub(crate) fn annotate_variant(
        &self,
        vcf_var: &VcfVariant,
    ) -> anyhow::Result<Option<crate::server::run::actix_server::seqvars_clinvar::ClinvarResultEntry>>
    {
        // Build key for RocksDB database
        let vcf_var = annonars::common::keys::Var::from(
            &vcf_var.chromosome,
            vcf_var.position,
            &vcf_var.reference,
            &vcf_var.alternative,
        );
        let key: Vec<u8> = vcf_var.clone().into();

        match self
            .db
            .get_cf(self.db.cf_handle("clinvar").as_ref().unwrap(), key)?
        {
            Some(raw_value) => {
                let record_list = annonars::pbs::clinvar::minimal::ExtractedVcvRecordList::decode(
                    &mut Cursor::new(&raw_value),
                )?;

                let mut clinvar_vcvs = Vec::new();
                let mut clinvar_germline_classifications = Vec::new();
                for clinvar_record in record_list.records.iter() {
                    let accession = clinvar_record.accession.as_ref().expect("must have VCV");
                    let vcv = format!("{}.{}", accession.accession, accession.version);
                    let classifications = clinvar_record
                        .classifications
                        .as_ref()
                        .expect("must have classifications");
                    if let Some(germline_classification) = &classifications.germline_classification
                    {
                        let description = germline_classification
                            .description
                            .as_ref()
                            .expect("description missing")
                            .to_string();
                        clinvar_vcvs.push(vcv);
                        clinvar_germline_classifications.push(description);
                    }
                }

                Ok(Some(
                    crate::server::run::actix_server::seqvars_clinvar::ClinvarResultEntry {
                        clinvar_vcv: clinvar_vcvs,
                        clinvar_germline_classification: clinvar_germline_classifications,
                    },
                ))
            }
            _ => Ok(None),
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ClinvarResult {
    pub vcv: Vec<String>,
    pub germline_classification: Vec<String>,
}
