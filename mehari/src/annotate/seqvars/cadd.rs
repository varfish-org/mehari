use std::path::Path;
use std::sync::Arc;

use crate::common::contig::ContigManager;
use crate::db::keys;
use annonars::common::contig::ContigDict;
use anyhow::{Context, Error};
use prost::Message;
use rocksdb::{DBWithThreadMode, MultiThreaded};

/// How the database stores the CADD scores.
#[derive(Debug)]
enum Layout {
    /// `annonars cadd import`: contig-ID keys and `CaddRecord` values.
    Cadd(ContigDict),
    /// `annonars tsv import`: chromosome-name keys and raw TSV lines, with the
    /// indices of the score columns.
    Tsv { raw_score: usize, phred: usize },
}

#[derive(Debug)]
pub struct CaddAnnotator {
    db: DBWithThreadMode<MultiThreaded>,
    contig_manager: Arc<ContigManager>,
    layout: Layout,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CaddResult {
    pub raw_score: f32,
    pub phred: f32,
}

impl CaddAnnotator {
    pub fn new(
        db: DBWithThreadMode<MultiThreaded>,
        contig_manager: Arc<ContigManager>,
        contig_dict: ContigDict,
    ) -> Self {
        Self {
            db,
            contig_manager,
            layout: Layout::Cadd(contig_dict),
        }
    }

    /// Open a CADD database built with `annonars cadd import` or `annonars tsv import`.
    pub(crate) fn from_path(
        path: impl AsRef<Path>,
        contig_manager: Arc<ContigManager>,
    ) -> anyhow::Result<Self> {
        tracing::info!("Opening CADD database at {}", path.as_ref().display());
        let options = rocksdb::Options::default();
        let cfs = rocksdb::DB::list_cf(&options, &path)?;
        if cfs.iter().any(|cf| cf == "cadd") {
            let db_cadd =
                rocksdb::DB::open_cf_for_read_only(&options, &path, ["meta", "cadd"], false)?;
            let contig_dict = annonars::seqvars::read_contig_dict(&db_cadd)?;
            Ok(Self::new(db_cadd, contig_manager, contig_dict))
        } else if cfs.iter().any(|cf| cf == "tsv_data") {
            let db_cadd =
                rocksdb::DB::open_cf_for_read_only(&options, &path, ["meta", "tsv_data"], false)?;
            let layout = read_tsv_layout(&db_cadd)?;
            Ok(Self {
                db: db_cadd,
                contig_manager,
                layout,
            })
        } else {
            anyhow::bail!(
                "no cadd or tsv_data column family in {}",
                path.as_ref().display()
            )
        }
    }

    pub fn annotate_record_cadd(&self, key: &[u8]) -> Result<Option<CaddResult>, Error> {
        if let Some(raw_value) = self
            .db
            .get_cf(self.db.cf_handle("cadd").as_ref().unwrap(), key)?
        {
            let record = annonars::pbs::seqvars::base::CaddRecord::decode(&raw_value[..])?;
            Ok(Some(CaddResult {
                raw_score: record.raw_score,
                phred: record.phred,
            }))
        } else {
            Ok(None)
        }
    }

    pub fn annotate(&self, vcf_var: &keys::Var) -> anyhow::Result<Option<CaddResult>> {
        let contig_dict = match &self.layout {
            Layout::Cadd(contig_dict) => contig_dict,
            Layout::Tsv { raw_score, phred } => {
                return self.annotate_tsv(vcf_var, *raw_score, *phred);
            }
        };

        // Normalize chrom to primary name to build standard key
        let chrom_std = self
            .contig_manager
            .get_primary_name(&vcf_var.chrom)
            .cloned()
            .unwrap_or_else(|| vcf_var.chrom.clone());

        if let Some(chrom_id) = contig_dict.id_of(&chrom_std) {
            let annonars_var: annonars::common::keys::Var = vcf_var.into();
            let key = annonars_var.encode_with_id(chrom_id);
            self.annotate_record_cadd(&key)
        } else {
            Ok(None)
        }
    }

    fn annotate_tsv(
        &self,
        vcf_var: &keys::Var,
        raw_score: usize,
        phred: usize,
    ) -> anyhow::Result<Option<CaddResult>> {
        let Some(chrom) = self.contig_manager.get_annonars_chrom(&vcf_var.chrom) else {
            return Ok(None);
        };
        let key: Vec<u8> = annonars::common::keys::Var::from(
            chrom,
            vcf_var.pos,
            &vcf_var.reference,
            &vcf_var.alternative,
        )
        .into();
        let cf = self
            .db
            .cf_handle("tsv_data")
            .context("missing tsv_data column family")?;
        let Some(raw_value) = self.db.get_cf(&cf, key)? else {
            return Ok(None);
        };

        // `annonars tsv import` has no delimiter option, it always uses tabs.
        let line = std::str::from_utf8(&raw_value)?;
        let fields: Vec<&str> = line.split('\t').collect();
        let score = |i: usize| -> anyhow::Result<f32> {
            let field = fields
                .get(i)
                .with_context(|| format!("CADD line has no column {i}"))?;
            field
                .parse()
                .with_context(|| format!("invalid CADD score {field:?}"))
        };
        Ok(Some(CaddResult {
            raw_score: score(raw_score)?,
            phred: score(phred)?,
        }))
    }
}

/// Locate the score columns in the schema that `annonars tsv import` stores.
fn read_tsv_layout(db: &rocksdb::DB) -> anyhow::Result<Layout> {
    let cf_meta = db.cf_handle("meta").context("missing meta column family")?;
    let schema = db
        .get_cf(&cf_meta, "db-schema")?
        .context("missing db-schema in meta column family")?;
    let schema: annonars::tsv::schema::FileSchema = serde_json::from_slice(&schema)?;
    let column = |name: &str| {
        schema
            .columns
            .iter()
            .position(|column| column.name == name)
            .with_context(|| format!("no {name} column in db-schema"))
    };
    Ok(Layout::Tsv {
        raw_score: column("RawScore")?,
        phred: column("PHRED")?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use temp_testdir::TempDir;

    #[test]
    fn test_cadd_import_and_annotate() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let fai_path = temp.join("ref.fa.fai");
        let input_path = temp.join("cadd_test.tsv");
        let output_path = temp.join("cadd_db");

        std::fs::write(&fai_path, "chr1\t248956422\t6\t60\t61\n")?;
        std::fs::write(
            &input_path,
            "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n\
             1\t10000\tC\tT\t1.5\t15.0\n\
             chr1\t10005\tA\tG\t2.3\t20.5\n",
        )?;

        let common_args = annonars::common::cli::Args {
            verbose: Default::default(),
        };
        let import_args = annonars::cadd::cli::import::Args {
            path_in_tsv: input_path.to_string_lossy().into(),
            path_reference_fai: fai_path.to_string_lossy().into(),
            assembly: "GRCh38".to_string(),
            path_out_rocksdb: output_path.to_string_lossy().into(),
            cf_name: "cadd".to_string(),
            path_wal_dir: None,
        };
        annonars::cadd::cli::import::run(&common_args, &import_args)?;

        let contig_manager = Arc::new(ContigManager::new("GRCh38"));
        let annotator = CaddAnnotator::from_path(&output_path, contig_manager)?;

        let var1 = keys::Var {
            chrom: "chr1".to_string(),
            pos: 10000,
            reference: "C".to_string(),
            alternative: "T".to_string(),
        };
        let res1 = annotator.annotate(&var1)?.unwrap();
        assert_eq!(res1.raw_score, 1.5);
        assert_eq!(res1.phred, 15.0);

        let var2 = keys::Var {
            chrom: "1".to_string(),
            pos: 10005,
            reference: "A".to_string(),
            alternative: "G".to_string(),
        };
        let res2 = annotator.annotate(&var2)?.unwrap();
        assert_eq!(res2.raw_score, 2.3);
        assert_eq!(res2.phred, 20.5);

        let var3 = keys::Var {
            chrom: "1".to_string(),
            pos: 99999,
            reference: "C".to_string(),
            alternative: "G".to_string(),
        };
        assert!(annotator.annotate(&var3)?.is_none());

        Ok(())
    }

    /// CADD databases built with the generic `annonars tsv import`.
    #[test]
    fn test_cadd_tsv_import_and_annotate() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let input_path = temp.join("cadd_test.tsv");
        let output_path = temp.join("cadd_db");

        std::fs::write(
            &input_path,
            "## CADD GRCh38-v1.6\n\
             #Chrom\tPos\tRef\tAlt\tType\tRawScore\tPHRED\n\
             1\t10000\tC\tT\tSNV\t1.5\t15.0\n\
             X\t10005\tA\tG\tSNV\t2.3\t20.5\n",
        )?;

        let common_args = annonars::common::cli::Args {
            verbose: Default::default(),
        };
        let import_args = annonars::tsv::cli::import::Args {
            genome_release: annonars::common::cli::GenomeRelease::Grch38,
            path_in_tsv: vec![input_path.to_string_lossy().into()],
            path_out_rocksdb: output_path.to_string_lossy().into(),
            path_schema_json: None,
            db_name: "CADD".to_string(),
            db_version: "1.6".to_string(),
            inference_row_count: 100,
            skip_row_count: 1,
            tbi_window_size: 100000,
            cf_name: "tsv_data".to_string(),
            path_wal_dir: None,
            col_chrom: "Chrom".to_string(),
            col_start: "Pos".to_string(),
            col_ref: "Ref".to_string(),
            col_alt: "Alt".to_string(),
            null_values: Vec::new(),
            add_default_null_values: true,
        };
        annonars::tsv::cli::import::run(&common_args, &import_args)?;

        let contig_manager = Arc::new(ContigManager::new("GRCh38"));
        let annotator = CaddAnnotator::from_path(&output_path, contig_manager)?;

        let var1 = keys::Var {
            chrom: "chr1".to_string(),
            pos: 10000,
            reference: "C".to_string(),
            alternative: "T".to_string(),
        };
        let res1 = annotator.annotate(&var1)?.unwrap();
        assert_eq!(res1.raw_score, 1.5);
        assert_eq!(res1.phred, 15.0);

        let var2 = keys::Var {
            chrom: "chrX".to_string(),
            pos: 10005,
            reference: "A".to_string(),
            alternative: "G".to_string(),
        };
        let res2 = annotator.annotate(&var2)?.unwrap();
        assert_eq!(res2.raw_score, 2.3);
        assert_eq!(res2.phred, 20.5);

        let var3 = keys::Var {
            chrom: "1".to_string(),
            pos: 99999,
            reference: "C".to_string(),
            alternative: "G".to_string(),
        };
        assert!(annotator.annotate(&var3)?.is_none());

        Ok(())
    }
}
