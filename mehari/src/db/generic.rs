use crate::common::Args as CommonArgs;
use crate::db::PipelineConfig;
use crate::db::keys::Var;
use crate::pbs::seqvars::GenericLookupRecord;
use anyhow::{Error, anyhow};
use clap::Parser;
use prost::Message;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

/// Arguments for the generic database construction command.
#[derive(Parser, Debug, Clone)]
#[command(about = "Construct generic lookup RocksDB database", long_about = None)]
pub struct Args {
    /// Genome assembly version (e.g., "grch37" or "grch38").
    #[arg(long, required = true)]
    pub assembly: String,

    /// Path(s) to the input data file(s) (VCF or TSV).
    #[arg(long, required = true)]
    pub input: Vec<PathBuf>,

    /// Path to the output RocksDB database directory.
    #[arg(long, required = true)]
    pub output: PathBuf,

    /// Internal identifier or name for the generic database being built.
    #[arg(long, required = true)]
    pub db_name: String,

    /// Input file format specification: either "vcf" or "tsv".
    #[arg(long, required = true)]
    pub format: String,

    /// Column identifier for chromosome coordinates (TSV mode only).
    #[arg(long, default_value = "chrom")]
    pub col_chrom: String,

    /// Column identifier for 1-based genomic position coordinates (TSV mode only).
    #[arg(long, default_value = "pos")]
    pub col_pos: String,

    /// Column identifier for the reference allele sequence (TSV mode only).
    #[arg(long, default_value = "ref")]
    pub col_ref: String,

    /// Column identifier for the alternative allele sequence (TSV mode only).
    #[arg(long, default_value = "alt")]
    pub col_alt: String,

    /// Specific column headers to explicitly extract as metadata values (TSV mode only).
    /// If omitted, all non-coordinate data columns are automatically captured.
    #[arg(long)]
    pub col_values: Option<Vec<String>>,

    /// Specific INFO keys to extract from VCF records (VCF mode only).
    /// If omitted, all encountered INFO flags/keys are captured.
    #[arg(long)]
    pub vcf_info_fields: Option<Vec<String>>,

    /// Number of records to chunk and commit per database write batch.
    #[arg(long, default_value = "100000")]
    pub batch_size: usize,
}

pub mod cli {
    pub use super::Args;
}

pub fn run(_common: &CommonArgs, args: &Args) -> Result<(), Error> {
    let mut extra_meta = HashMap::new();
    extra_meta.insert("db_name".to_string(), args.db_name.clone());

    let config = PipelineConfig {
        assembly: &args.assembly,
        input: &args.input,
        output: &args.output,
        batch_size: args.batch_size,
        db_type: "generic",
        schema_version: "1.0",
        extra_meta,
    };

    if args.format.to_lowercase() == "vcf" {
        crate::db::run_vcf_pipeline(
            config,
            None::<fn(&mut noodles::vcf::Header)>,
            |record, contig_manager| {
                let mut kvs = Vec::new();
                let mut local_keys = HashSet::new();

                let chrom = record.reference_sequence_name();
                let pos = match record.variant_start() {
                    Some(start) => start.get() as i32,
                    None => return Ok((kvs, local_keys)),
                };

                let mut fields = HashMap::new();
                if let Some(keys) = &args.vcf_info_fields {
                    for k in keys {
                        if let Some(val) = record.info().get(k).flatten()
                            && let Some(v_str) = crate::db::get_info_string(val)
                        {
                            fields.insert(k.clone(), v_str);
                        }
                    }
                } else {
                    for (k, val) in record.info().as_ref() {
                        if let Some(val_inner) = val
                            && let Some(v_str) = crate::db::get_info_string(val_inner)
                        {
                            fields.insert(k.to_string(), v_str);
                        }
                    }
                }

                for k in fields.keys() {
                    local_keys.insert(k.clone());
                }
                let chrom_std = contig_manager
                    .get_primary_name(chrom)
                    .cloned()
                    .unwrap_or_else(|| chrom.to_string());
                let reference = record.reference_bases();

                for alt in record.alternate_bases().as_ref() {
                    let var = Var {
                        chrom: chrom_std.clone(),
                        pos,
                        reference: reference.to_string(),
                        alternative: alt.to_string(),
                    };
                    let key: Vec<u8> = var.clone().into();
                    let record_pb = GenericLookupRecord {
                        fields: fields.clone(),
                    };

                    let mut value = Vec::new();
                    record_pb.encode(&mut value)?;

                    let var_label = format!(
                        "{}:{}{}>{}",
                        var.chrom, var.pos, var.reference, var.alternative
                    );
                    kvs.push((key, value, var_label));
                }
                Ok((kvs, local_keys))
            },
        )
    } else {
        crate::db::run_tsv_pipeline(
            config,
            crate::db::open_tsv_reader,
            |record, headers_record, contig_manager| {
                let record_map: HashMap<String, String> =
                    record.deserialize(Some(headers_record))?;

                let chrom = record_map
                    .get(&args.col_chrom)
                    .ok_or_else(|| anyhow!("Missing Chromosome column"))?;
                let pos: i32 = record_map
                    .get(&args.col_pos)
                    .ok_or_else(|| anyhow!("Missing Position column"))?
                    .parse()?;
                let reference = record_map
                    .get(&args.col_ref)
                    .ok_or_else(|| anyhow!("Missing Reference column"))?;
                let alternative = record_map
                    .get(&args.col_alt)
                    .ok_or_else(|| anyhow!("Missing Alternative column"))?;

                let chrom_std = contig_manager
                    .get_primary_name(chrom)
                    .cloned()
                    .unwrap_or_else(|| chrom.to_string());
                let mut fields = HashMap::new();

                if let Some(vals) = &args.col_values {
                    for v in vals {
                        let val = record_map
                            .get(v)
                            .ok_or_else(|| anyhow!("Value column not found"))?;
                        fields.insert(v.clone(), val.clone());
                    }
                } else {
                    for (k, v) in &record_map {
                        if k != &args.col_chrom
                            && k != &args.col_pos
                            && k != &args.col_ref
                            && k != &args.col_alt
                        {
                            fields.insert(k.clone(), v.clone());
                        }
                    }
                }

                let local_keys: HashSet<String> = fields.keys().cloned().collect();
                let var = Var {
                    chrom: chrom_std,
                    pos,
                    reference: reference.to_string(),
                    alternative: alternative.to_string(),
                };
                let key: Vec<u8> = var.clone().into();
                let record_pb = GenericLookupRecord { fields };

                let mut value = Vec::new();
                record_pb.encode(&mut value)?;

                let var_label = format!(
                    "{}:{}{}>{}",
                    var.chrom, var.pos, var.reference, var.alternative
                );
                Ok((vec![(key, value, var_label)], local_keys))
            },
        )
    }
}
