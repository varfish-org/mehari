use crate::annotate::seqvars::consequence::VcfVariant;
use crate::common::contig::ContigManager;
use crate::db::keys;
use annonars::freqs::serialized::{auto, mt, xy};
use anyhow::Error;
use rocksdb::{DBWithThreadMode, MultiThreaded};
use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::path::Path;
use std::sync::Arc;

#[derive(Debug)]
pub struct FrequencyAnnotator {
    db: DBWithThreadMode<MultiThreaded>,
    contig_manager: Arc<ContigManager>,
}

impl FrequencyAnnotator {
    pub fn new(db: DBWithThreadMode<MultiThreaded>, contig_manager: Arc<ContigManager>) -> Self {
        Self { db, contig_manager }
    }

    pub(crate) fn from_path(
        path: impl AsRef<Path> + Display,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        // Open the frequency RocksDB database in read only mode.
        tracing::info!("Opening frequency database");
        tracing::debug!("RocksDB path = {}", &path);
        let options = rocksdb::Options::default();
        let db_freq = rocksdb::DB::open_cf_for_read_only(
            &options,
            &path,
            ["meta", "autosomal", "gonosomal", "mitochondrial"],
            false,
        )?;
        Ok(Self::new(db_freq, contig_manager))
    }

    /// Annotate record on autosomal chromosome with gnomAD exomes/genomes.
    pub fn annotate_record_auto(&self, key: &[u8]) -> Result<Option<FreqResult>, Error> {
        if let Some(freq) = self
            .db
            .get_cf(self.db.cf_handle("autosomal").as_ref().unwrap(), key)?
        {
            let auto_record = auto::Record::from_buf(&freq);
            Ok(Some(FreqResult {
                gnomad_exomes_an: auto_record.gnomad_exomes.an as i32,
                gnomad_exomes_hom: auto_record.gnomad_exomes.ac_hom as i32,
                gnomad_exomes_het: auto_record.gnomad_exomes.ac_het as i32,
                gnomad_genomes_an: auto_record.gnomad_genomes.an as i32,
                gnomad_genomes_hom: auto_record.gnomad_genomes.ac_hom as i32,
                gnomad_genomes_het: auto_record.gnomad_genomes.ac_het as i32,
                ..Default::default()
            }))
        } else {
            Ok(None)
        }
    }

    /// Annotate record on gonosomal chromosome with gnomAD exomes/genomes.
    pub fn annotate_record_xy(&self, key: &[u8]) -> Result<Option<FreqResult>, Error> {
        if let Some(freq) = self
            .db
            .get_cf(self.db.cf_handle("gonosomal").as_ref().unwrap(), key)?
        {
            let xy_record = xy::Record::from_buf(&freq);
            Ok(Some(FreqResult {
                gnomad_exomes_an: xy_record.gnomad_exomes.an as i32,
                gnomad_exomes_hom: xy_record.gnomad_exomes.ac_hom as i32,
                gnomad_exomes_het: xy_record.gnomad_exomes.ac_het as i32,
                gnomad_exomes_hemi: Some(xy_record.gnomad_exomes.ac_hemi as i32),
                gnomad_genomes_an: xy_record.gnomad_genomes.an as i32,
                gnomad_genomes_hom: xy_record.gnomad_genomes.ac_hom as i32,
                gnomad_genomes_het: xy_record.gnomad_genomes.ac_het as i32,
                gnomad_genomes_hemi: Some(xy_record.gnomad_genomes.ac_hemi as i32),
                ..Default::default()
            }))
        } else {
            Ok(None)
        }
    }

    /// Annotate record on mitochondrial genome with gnomAD mtDNA and HelixMtDb.
    pub fn annotate_record_mt(&self, key: &[u8]) -> Result<Option<FreqResult>, Error> {
        if let Some(freq) = self
            .db
            .get_cf(self.db.cf_handle("mitochondrial").as_ref().unwrap(), key)?
        {
            let mt_record = mt::Record::from_buf(&freq);
            Ok(Some(FreqResult {
                helix_an: Some(mt_record.helixmtdb.an as i32),
                helix_hom: Some(mt_record.helixmtdb.ac_hom as i32),
                helix_het: Some(mt_record.helixmtdb.ac_het as i32),
                gnomad_genomes_an: mt_record.gnomad_mtdna.an as i32,
                gnomad_genomes_hom: mt_record.gnomad_mtdna.ac_hom as i32,
                gnomad_genomes_het: mt_record.gnomad_mtdna.ac_het as i32,
                ..Default::default()
            }))
        } else {
            Ok(None)
        }
    }

    pub(crate) fn annotate(&self, vcf_var: &keys::Var) -> anyhow::Result<Option<FreqResult>> {
        let contig_manager = &self.contig_manager;
        // Only attempt lookups into RocksDB for canonical contigs.
        if contig_manager.is_homo_sapiens_canonical_alias(vcf_var.chrom.as_str()) {
            // Build key for RocksDB database from `vcf_var`.
            let annonars_var: annonars::common::keys::Var = vcf_var.into();
            let key: Vec<u8> = annonars_var.into();

            // Annotate with frequency.
            if contig_manager.is_autosomal_alias(&vcf_var.chrom) {
                return self.annotate_record_auto(&key);
            } else if contig_manager.is_gonosomal_alias(&vcf_var.chrom) {
                return self.annotate_record_xy(&key);
            } else if contig_manager.is_mitochondrial_alias(&vcf_var.chrom) {
                return self.annotate_record_mt(&key);
            }
        }
        Ok(None)
    }

    #[cfg(feature = "server")]
    pub(crate) fn annotate_variant(
        &self,
        vcf_var: &VcfVariant,
    ) -> anyhow::Result<
        Option<crate::server::run::actix_server::seqvars_frequencies::FrequencyResultEntry>,
    > {
        let contig_manager = &self.contig_manager;
        // Only attempt lookups into RocksDB for canonical contigs.
        if !contig_manager.is_homo_sapiens_canonical_alias(&vcf_var.chromosome) {
            return Ok(None);
        }

        // Build key for RocksDB database
        let vcf_var = annonars::common::keys::Var::from(
            &vcf_var.chromosome.clone(),
            vcf_var.position,
            &vcf_var.reference.clone(),
            &vcf_var.alternative.clone(),
        );
        let key: Vec<u8> = vcf_var.clone().into();
        use crate::server::run::actix_server::seqvars_frequencies::*;
        use anyhow::anyhow;
        // Annotate with frequency.
        if contig_manager.is_autosomal_alias(vcf_var.chrom.as_str()) {
            match self
                .db
                .get_cf(self.db.cf_handle("autosomal").as_ref().unwrap(), key)?
            {
                Some(freq) => {
                    let val = auto::Record::from_buf(&freq);
                    Ok(Some(FrequencyResultEntry::Autosomal(
                        AutosomalResultEntry {
                            gnomad_exomes_an: val.gnomad_exomes.an,
                            gnomad_exomes_hom: val.gnomad_exomes.ac_hom,
                            gnomad_exomes_het: val.gnomad_exomes.ac_het,
                            gnomad_genomes_an: val.gnomad_genomes.an,
                            gnomad_genomes_hom: val.gnomad_genomes.ac_hom,
                            gnomad_genomes_het: val.gnomad_genomes.ac_het,
                        },
                    )))
                }
                _ => Err(anyhow!("No frequency data found for variant {:?}", vcf_var)),
            }
        } else if contig_manager.is_gonosomal_alias(vcf_var.chrom.as_str()) {
            match self
                .db
                .get_cf(self.db.cf_handle("gonosomal").as_ref().unwrap(), key)?
            {
                Some(freq) => {
                    let val = xy::Record::from_buf(&freq);
                    Ok(Some(FrequencyResultEntry::Gonosomal(
                        GonosomalResultEntry {
                            gnomad_exomes_an: val.gnomad_exomes.an,
                            gnomad_exomes_hom: val.gnomad_exomes.ac_hom,
                            gnomad_exomes_het: val.gnomad_exomes.ac_het,
                            gnomad_exomes_hemi: val.gnomad_exomes.ac_hemi,
                            gnomad_genomes_an: val.gnomad_genomes.an,
                            gnomad_genomes_hom: val.gnomad_genomes.ac_hom,
                            gnomad_genomes_het: val.gnomad_genomes.ac_het,
                            gnomad_genomes_hemi: val.gnomad_genomes.ac_hemi,
                        },
                    )))
                }
                _ => Err(anyhow!("No frequency data found for variant {:?}", vcf_var)),
            }
        } else if contig_manager.is_mitochondrial_alias(vcf_var.chrom.as_str()) {
            match self
                .db
                .get_cf(self.db.cf_handle("mitochondrial").as_ref().unwrap(), key)?
            {
                Some(freq) => {
                    let val = mt::Record::from_buf(&freq);
                    Ok(Some(FrequencyResultEntry::Mitochondrial(
                        MitochondrialResultEntry {
                            helix_an: val.helixmtdb.an,
                            helix_hom: val.helixmtdb.ac_hom,
                            helix_het: val.helixmtdb.ac_het,
                            gnomad_genomes_an: val.gnomad_mtdna.an,
                            gnomad_genomes_hom: val.gnomad_mtdna.ac_hom,
                            gnomad_genomes_het: val.gnomad_mtdna.ac_het,
                        },
                    )))
                }
                _ => Err(anyhow!("No frequency data found for variant {:?}", vcf_var)),
            }
        } else {
            tracing::trace!(
                "Record @{:?} on non-canonical chromosome, skipping.",
                &vcf_var
            );
            Ok(None)
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct FreqResult {
    pub gnomad_exomes_an: i32,
    pub gnomad_exomes_hom: i32,
    pub gnomad_exomes_het: i32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomad_exomes_hemi: Option<i32>,

    pub gnomad_genomes_an: i32,
    pub gnomad_genomes_hom: i32,
    pub gnomad_genomes_het: i32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomad_genomes_hemi: Option<i32>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub helix_an: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub helix_hom: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub helix_het: Option<i32>,
}
