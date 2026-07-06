use anyhow::{Result, anyhow};
use clap::Parser;
use prost::Message;
use rocksdb::{Direction, IteratorMode, WriteBatch};
use rustc_hash::FxHashMap;
use std::collections::{BTreeMap, HashSet};
use std::path::PathBuf;
use std::sync::Arc;

use crate::common::Args as CommonArgs;
use crate::common::contig::ContigManager;
use crate::db::keys::Var;
use crate::db::{finalize_db, open_db_for_read, tune_options_for_build, write_contig_dictionary};
use crate::pbs::seqvars::{CaddRecord, DbsnpRecord, IntegratedVariantRecord, SpliceAiRecord};

#[derive(Parser, Debug, Clone)]
#[command(about = "Merge individual track databases into a unified WORM layout", long_about = None)]
pub struct Args {
    #[arg(long)]
    pub cadd: Option<PathBuf>,

    #[arg(long)]
    pub spliceai: Option<PathBuf>,

    #[arg(long)]
    pub dbsnp: Option<PathBuf>,

    #[arg(long, short)]
    pub output: PathBuf,

    #[arg(long, default_value = "GRCh38")]
    pub assembly: String,

    #[arg(long, default_value = "100000")]
    pub batch_size: usize,
}

pub mod cli {
    pub use super::Args;
}

struct SourceTrack {
    db: rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>,
    cf_name: String,
    contig_to_id: FxHashMap<String, u32>,
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct Coordinates {
    pos: i32,
    reference: String,
    alternative: String,
}

pub fn run(_common: &CommonArgs, args: &Args) -> Result<()> {
    tracing::info!(
        "Initializing Track Merger Pipeline for output: {:?}",
        args.output
    );
    let contig_manager = ContigManager::new(&args.assembly);

    let mut sources: BTreeMap<String, SourceTrack> = BTreeMap::new();
    let mut all_contigs: HashSet<String> = HashSet::new();

    // load active source databases and extract their contig layouts
    let configurations = [
        ("cadd", &args.cadd),
        ("spliceai", &args.spliceai),
        ("dbsnp", &args.dbsnp),
    ];

    for (name, path_opt) in configurations {
        if let Some(path) = path_opt {
            let db = open_db_for_read(path, name)?;

            let contig_to_id: FxHashMap<String, u32> = {
                let cf_meta = db
                    .cf_handle("meta")
                    .ok_or_else(|| anyhow!("Missing meta CF in {}", name))?;

                match db.get_cf(&cf_meta, b"contig_dictionary")? {
                    Some(bytes) => serde_json::from_slice(&bytes)?,
                    None => anyhow::bail!(
                        "Database [{}] is missing its contig_dictionary context map",
                        name
                    ),
                }
            };

            for contig in contig_to_id.keys() {
                all_contigs.insert(contig.clone());
            }

            sources.insert(
                name.to_string(),
                SourceTrack {
                    db,
                    cf_name: name.to_string(),
                    contig_to_id,
                },
            );
        }
    }

    if sources.is_empty() {
        anyhow::bail!(
            "No input track databases specified. Provide at least one track (--cadd, --spliceai, --dbsnp)."
        );
    }

    // open output database
    let mut out_options = rocksdb::Options::default();
    out_options = tune_options_for_build(out_options, None);
    let out_db = rocksdb::DB::open_cf(&out_options, &args.output, vec!["meta", "unified_tracks"])?;
    let cf_out = out_db.cf_handle("unified_tracks").unwrap();

    let output_contig_to_id = Arc::new(std::sync::RwLock::new(FxHashMap::default()));
    let mut sorted_contigs: Vec<String> = all_contigs.into_iter().collect();
    sorted_contigs.sort_by(|a, b| {
        let a_num = contig_manager.get_chrom_no(a);
        let b_num = contig_manager.get_chrom_no(b);

        match (a_num, b_num) {
            // both are canonical (1-22, X=23, Y=24, MT=25) -> sort numerically
            (Some(an), Some(bn)) => an.cmp(&bn),

            // canonical chromosomes take precedence over unplaced/custom scaffolds
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,

            // both are custom/unplaced scaffolds -> fall back to alphabetical sort
            (None, None) => a.cmp(b),
        }
    });

    let mut current_batch = WriteBatch::default();
    let mut batch_count = 0;

    // process contig by contig to guarantee sequential key formatting
    for contig_name in &sorted_contigs {
        tracing::info!(
            "Merging coordinate data streams for contig: {}",
            contig_name
        );

        let output_contig_id = {
            let mut map = output_contig_to_id.write().unwrap();
            let id = map.len() as u32;
            *map.entry(contig_name.clone()).or_insert(id)
        };

        // instantiate prefix iterators for every database containing this contig
        let mut active_iters = Vec::new();
        for (track_name, src) in &sources {
            if let Some(&local_id) = src.contig_to_id.get(contig_name) {
                let cf = src.db.cf_handle(&src.cf_name).unwrap();
                let prefix = &local_id.to_be_bytes()[1..4];

                let mut iter = src.db.iterator_cf_opt(
                    &cf,
                    rocksdb::ReadOptions::default(),
                    IteratorMode::From(prefix, Direction::Forward),
                );

                // verify prefix match
                if let Some(Ok((k, v))) = iter.next()
                    && k.starts_with(prefix)
                {
                    active_iters.push((track_name.clone(), (k, v), iter, prefix.to_vec()));
                }
            }
        }

        // multi-way merge
        while !active_iters.is_empty() {
            // find lowest genomic coordinate across all DBs
            let mut lowest_coords: Option<Coordinates> = None;
            for (_, (key, _), _, _) in &active_iters {
                let pos = i32::from_be_bytes(key[3..7].try_into()?);
                let alleles_buf = &key[7..];
                let null_idx = alleles_buf.iter().position(|&b| b == 0x00).unwrap();
                let reference = String::from_utf8_lossy(&alleles_buf[0..null_idx]).into_owned();
                let alternative =
                    String::from_utf8_lossy(&alleles_buf[null_idx + 1..]).into_owned();

                let item_coords = Coordinates {
                    pos,
                    reference,
                    alternative,
                };
                match &lowest_coords {
                    None => lowest_coords = Some(item_coords),
                    Some(lowest) if item_coords < *lowest => lowest_coords = Some(item_coords),
                    _ => {}
                }
            }

            let current_target = lowest_coords.unwrap();
            let mut unified_record = IntegratedVariantRecord::default();

            // get values matching this coordinate from all DBs
            let mut i = 0;
            while i < active_iters.len() {
                let track_match = {
                    let (ref _name, (ref key, ref _val), _, _) = active_iters[i];
                    let pos = i32::from_be_bytes(key[3..7].try_into()?);
                    let alleles_buf = &key[7..];
                    let null_idx = alleles_buf.iter().position(|&b| b == 0x00).unwrap();
                    let reference = String::from_utf8_lossy(&alleles_buf[0..null_idx]);
                    let alternative = String::from_utf8_lossy(&alleles_buf[null_idx + 1..]);

                    pos == current_target.pos
                        && reference == current_target.reference
                        && alternative == current_target.alternative
                };

                if track_match {
                    let (name, (_, val_bytes), mut iter, prefix) = active_iters.remove(i);

                    match name.as_str() {
                        "cadd" => unified_record.cadd = Some(CaddRecord::decode(&val_bytes[..])?),
                        "spliceai" => {
                            unified_record.splice_ai = Some(SpliceAiRecord::decode(&val_bytes[..])?)
                        }
                        "dbsnp" => {
                            unified_record.dbsnp = Some(DbsnpRecord::decode(&val_bytes[..])?)
                        }
                        _ => {}
                    }

                    if let Some(Ok((next_k, next_v))) = iter.next()
                        && next_k.starts_with(&prefix)
                    {
                        active_iters.insert(i, (name, (next_k, next_v), iter, prefix));
                        i += 1;
                        continue;
                    }
                } else {
                    i += 1;
                }
            }

            // append to batch
            let out_var = Var::new(
                contig_name.clone(),
                current_target.pos,
                current_target.reference,
                current_target.alternative,
            );
            let out_key = out_var.encode_with_id(output_contig_id);

            let mut out_val = Vec::with_capacity(unified_record.encoded_len());
            unified_record.encode(&mut out_val)?;

            current_batch.put_cf(&cf_out, &out_key, &out_val);
            batch_count += 1;

            if batch_count >= args.batch_size {
                out_db.write(current_batch)?;
                current_batch = WriteBatch::default();
                batch_count = 0;
            }
        }
    }

    if batch_count > 0 {
        out_db.write(current_batch)?;
    }
    drop(cf_out);

    // complete finalization and register structural maps
    finalize_db(&out_db, &["unified_tracks", "meta"], "unified")?;
    drop(out_db);

    write_contig_dictionary(&args.output, "unified_tracks", &output_contig_to_id)?;
    tracing::info!("Successfully generated unified track database directory.");
    Ok(())
}
