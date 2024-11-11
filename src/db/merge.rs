use crate::annotate::seqvars::load_tx_db;
use crate::common;
use crate::db::create::write_tx_db;
use crate::pbs::txs::TxSeqDatabase;
use clap::Parser;
use itertools::Itertools;
use std::collections::HashMap;

/// Merge two or more mehari transcript databases.
#[derive(Parser, Debug)]
#[command(about = "Merge two or more mehari transcript databases")]
pub struct Args {
    /// The input transcript databases to merge.
    #[clap(long, required = true)]
    pub database: Vec<String>,

    /// Output file to write the merged transcript database to.
    #[clap(long)]
    pub output: String,
}

pub fn merge_transcript_databases(
    mut databases: Vec<TxSeqDatabase>,
) -> anyhow::Result<TxSeqDatabase> {
    if let Some((first, others)) = databases.split_first_mut() {
        if !others.is_empty() {
            tracing::info!("Merging multiple transcript databases into one");
        }

        // Ensure that all databases use the same assembly.
        if !first.source_version.iter().map(|v| v.assembly).all_equal() {
            return Err(anyhow::anyhow!(
                "Inconsistent assembly versions in first database"
            ));
        }
        let assembly = first
            .source_version
            .first()
            .expect("At least one source_version entry expected")
            .assembly;

        if !others
            .iter()
            .all(|db| db.source_version.iter().all(|v| v.assembly == assembly))
        {
            return Err(anyhow::anyhow!(
                "All databases must use the same assembly version"
            ));
        }

        let seq_db = first.seq_db.as_mut().unwrap();
        let tx_db = first.tx_db.as_mut().unwrap();

        let mut prev_max_idx = *seq_db.aliases_idx.iter().max().unwrap();
        for other in others.iter_mut() {
            let mut other_seq_db = other.seq_db.take().unwrap();
            // Merge the sequence database.
            seq_db.seqs.append(&mut other_seq_db.seqs);
            // Merge the aliases.
            seq_db.aliases.append(&mut other_seq_db.aliases);
            // Merge the alias index, ensuring that they are unique by adding the max index.
            seq_db.aliases_idx.extend(
                other_seq_db
                    .aliases_idx
                    .drain(..)
                    .map(|idx| idx + prev_max_idx + 1),
            );
            prev_max_idx = *seq_db.aliases_idx.iter().max().unwrap();

            let mut other_tx_db = other.tx_db.take().unwrap();
            // Merge the transcript database.
            tx_db.transcripts.append(&mut other_tx_db.transcripts);

            // Merge the gene to transcript mapping:
            // 1. Create a map from gene ID to gene to transcript mapping.
            // 2. Iterate over the gene to transcript mappings of the first transcript database (target db).
            //    For each gene to transcript mapping, check if there is a corresponding mapping in the other database.
            //    If so, merge the transcript IDs.
            let mut other_gene_to_tx: HashMap<_, _> = other_tx_db
                .gene_to_tx
                .into_iter()
                .map(|g2t| (g2t.gene_id.clone(), g2t))
                .collect();
            for g2tx in tx_db.gene_to_tx.iter_mut() {
                if let Some(other_g2tx) = other_gene_to_tx.remove(&g2tx.gene_id) {
                    g2tx.tx_ids.extend(other_g2tx.tx_ids);
                    if let Some(f) = g2tx.filtered.as_mut() {
                        *f &= other_g2tx.filtered.unwrap_or(false);
                    }
                    if let Some(r) = g2tx.filter_reason.as_mut() {
                        *r |= other_g2tx.filter_reason.unwrap_or(0);
                    }
                }
            }
            // Add the remaining gene to transcript mappings from the other database.
            for other_g2tx in other_gene_to_tx.into_values() {
                tx_db.gene_to_tx.push(other_g2tx);
            }
        }

        let result = std::mem::take(first);
        Ok(result)
    } else {
        anyhow::bail!("No transcript databases provided");
    }
}

pub fn run(_common_args: &common::Args, args: &Args) -> Result<(), anyhow::Error> {
    if args.database.len() < 2 {
        anyhow::bail!("At least two transcript databases are required to merge");
    }
    tracing::info!("Loading transcript databases");
    let tx_dbs = args
        .database
        .iter()
        .map(load_tx_db)
        .collect::<anyhow::Result<Vec<_>>>()?;
    tracing::info!("Merging transcript databases");
    let tx_db = merge_transcript_databases(tx_dbs)?;
    tracing::info!("Writing merged transcript database");
    write_tx_db(tx_db, &args.output)?;
    tracing::info!("Done loading, merging and writing transcript databases");
    Ok(())
}
