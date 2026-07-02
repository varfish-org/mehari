use crate::common::Args as CommonArgs;
use crate::db::keys::Var;
use crate::db::{PipelineConfig, write_contig_dictionary};
use crate::pbs::seqvars::GenericLookupRecord;
use anyhow::{Error, anyhow};
use clap::Parser;
use prost::Message;
use rustc_hash::FxHashMap;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, RwLock};

/// Arguments for the generic database construction command.
#[derive(Parser, Debug, Clone)]
#[command(about = "Construct generic lookup RocksDB database", long_about = None)]
pub struct Args {
    #[command(flatten)]
    pub common: crate::db::CommonPipelineArgs,

    #[arg(long, required = true)]
    pub db_name: String,

    #[arg(long, required = true)]
    pub format: String,

    #[arg(long, default_value = "chrom")]
    pub col_chrom: String,

    #[arg(long, default_value = "pos")]
    pub col_pos: String,

    #[arg(long, default_value = "ref")]
    pub col_ref: String,

    #[arg(long, default_value = "alt")]
    pub col_alt: String,

    #[arg(long)]
    pub col_values: Option<Vec<String>>,

    #[arg(long)]
    pub vcf_info_fields: Option<Vec<String>>,
}

pub mod cli {
    pub use super::Args;
}

pub fn run(_common: &CommonArgs, args: &Args) -> Result<(), Error> {
    let mut extra_meta = HashMap::new();
    extra_meta.insert("db_name".to_string(), args.db_name.clone());

    let config = PipelineConfig {
        assembly: &args.common.assembly,
        input: &args.common.input,
        output: &args.common.output,
        batch_size: args.common.batch_size,
        no_progress: args.common.no_progress,
        threads: args.common.threads,
        db_type: "generic",
        schema_version: "1.0",
        extra_meta,
    };

    let chrom_to_id = Arc::new(RwLock::new(FxHashMap::default()));
    let chrom_to_id_vcf = Arc::clone(&chrom_to_id);
    let chrom_to_id_tsv = Arc::clone(&chrom_to_id);

    if args.format.to_lowercase() == "vcf" {
        crate::db::run_vcf_pipeline(
            config,
            None::<fn(&mut noodles::vcf::Header)>,
            move |record, contig_manager| {
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

                let (chrom_std, chrom_id) =
                    crate::db::get_or_intern_contig(chrom, contig_manager, &chrom_to_id_vcf);

                let reference = record.reference_bases();

                for alt in record.alternate_bases().as_ref() {
                    let var = Var::new(
                        chrom_std.clone(),
                        pos,
                        reference.to_string(),
                        alt.to_string(),
                    );
                    let key = var.encode_with_id(chrom_id);

                    let record_pb = GenericLookupRecord {
                        fields: fields.clone(),
                    };

                    let mut value = Vec::new();
                    record_pb.encode(&mut value)?;

                    let var_label = format!("{}:{}{}>{}", chrom_std, pos, reference, alt);
                    kvs.push((key, value, var_label));
                }
                Ok((kvs, local_keys))
            },
        )?;
    } else {
        crate::db::run_tsv_pipeline(
            config,
            crate::db::open_tsv_reader,
            move |record, headers_record, contig_manager| {
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

                let (chrom_std, chrom_id) =
                    crate::db::get_or_intern_contig(chrom, contig_manager, &chrom_to_id_tsv);

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

                let var = Var::new(
                    chrom_std.clone(),
                    pos,
                    reference.to_string(),
                    alternative.to_string(),
                );
                let key = var.encode_with_id(chrom_id);

                let record_pb = GenericLookupRecord { fields };

                let mut value = Vec::new();
                record_pb.encode(&mut value)?;

                let var_label = format!("{}:{}{}>{}", chrom_std, pos, reference, alternative);
                Ok((vec![(key, value, var_label)], local_keys))
            },
        )?;
    }

    tracing::info!("Writing generic database contig index metadata mapping...");
    write_contig_dictionary(&args.common.output, "generic", &chrom_to_id)?;

    Ok(())
}
