use crate::annotate::seqvars::consequence::load_tx_db;
use crate::common;
use crate::db::transcripts::create::write_tx_db;
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

    /// ZSTD compression level to use.
    #[arg(long, default_value = "19")]
    compression_level: i32,
}

pub fn merge_transcript_databases(
    mut databases: Vec<TxSeqDatabase>,
) -> anyhow::Result<TxSeqDatabase> {
    if let Some((first, others)) = databases.split_first_mut() {
        if !others.is_empty() {
            tracing::info!("Merging multiple transcript databases into one");
        }

        // Ensure that all databases use the same assembly.
        if !first
            .source_version
            .iter()
            .map(|v| v.assembly.clone())
            .all_equal()
        {
            return Err(anyhow::anyhow!(
                "Inconsistent assembly versions in first database"
            ));
        }
        let assembly = first
            .source_version
            .first()
            .expect("At least one source_version entry expected")
            .assembly
            .clone();

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
                    // Keep the HGNC ID of the first database. Warn as the importer does if the
                    // other database names another one.
                    if let Some(other_hgnc_id) = &other_g2tx.hgnc_id
                        && g2tx.hgnc_id.as_ref() != Some(other_hgnc_id)
                    {
                        tracing::warn!(
                            "gene {} has HGNC ID {}, but another database names {}",
                            g2tx.gene_id,
                            g2tx.hgnc_id.as_deref().unwrap_or("none"),
                            other_hgnc_id
                        );
                    }
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
    write_tx_db(
        tx_db,
        &args.output,
        args.compression_level,
        &common::progress::NoProgress,
    )?;
    tracing::info!("Done loading, merging and writing transcript databases");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::merge_transcript_databases;
    use crate::pbs::txs::{GeneToTxId, SequenceDb, SourceVersion, TranscriptDb, TxSeqDatabase};

    /// A database with one `GeneToTxId` for each `(gene_id, legacy_gene_id, tx_id)`.
    #[allow(deprecated)]
    fn db(genes: &[(&str, &str, &str)]) -> TxSeqDatabase {
        TxSeqDatabase {
            tx_db: Some(TranscriptDb {
                gene_to_tx: genes
                    .iter()
                    .map(|&(gene_id, legacy_gene_id, tx_id)| GeneToTxId {
                        gene_id: gene_id.into(),
                        legacy_gene_id: legacy_gene_id.into(),
                        tx_ids: vec![tx_id.into()],
                        ..Default::default()
                    })
                    .collect(),
                ..Default::default()
            }),
            seq_db: Some(SequenceDb {
                aliases: vec!["NM_000001.1".into()],
                aliases_idx: vec![0],
                seqs: vec!["A".into()],
            }),
            version: Some("test".into()),
            source_version: vec![SourceVersion {
                assembly: "GRCh38".into(),
                ..Default::default()
            }],
        }
    }

    /// The `(gene_id, legacy_gene_id, tx_ids)` of the merged database, sorted.
    #[allow(deprecated)]
    fn merged_genes(databases: Vec<TxSeqDatabase>) -> Vec<(String, String, Vec<String>)> {
        let mut genes: Vec<_> = merge_transcript_databases(databases)
            .unwrap()
            .tx_db
            .unwrap()
            .gene_to_tx
            .into_iter()
            .map(|gene| (gene.gene_id, gene.legacy_gene_id, gene.tx_ids))
            .collect();
        genes.sort();
        genes
    }

    fn gene(gene_id: &str, legacy_gene_id: &str, tx_ids: &[&str]) -> (String, String, Vec<String>) {
        (
            gene_id.into(),
            legacy_gene_id.into(),
            tx_ids.iter().map(|&tx_id| tx_id.into()).collect(),
        )
    }

    /// The `hgnc_id` values of the merged database. Gene `11` has `first` in the first
    /// database and `other` in the other one.
    fn merged_hgnc_ids(first: Option<&str>, other: Option<&str>) -> Vec<Option<String>> {
        let databases = [(first, "NM_000001.1"), (other, "NM_000002.1")]
            .into_iter()
            .map(|(hgnc_id, tx_id)| {
                let mut db = db(&[("11", "HGNC:1", tx_id)]);
                db.tx_db.as_mut().unwrap().gene_to_tx[0].hgnc_id = hgnc_id.map(Into::into);
                db
            })
            .collect();
        merge_transcript_databases(databases)
            .unwrap()
            .tx_db
            .unwrap()
            .gene_to_tx
            .into_iter()
            .map(|gene| gene.hgnc_id)
            .collect()
    }

    /// A gene keeps the HGNC ID of the first database. The merge warns if the other
    /// database names another one, also if the first names none.
    #[tracing_test::traced_test]
    #[test]
    fn merge_warns_about_another_hgnc_id() {
        assert_eq!(
            merged_hgnc_ids(Some("HGNC:1"), Some("HGNC:2")),
            [Some("HGNC:1".to_string())]
        );
        assert!(logs_contain(
            "gene 11 has HGNC ID HGNC:1, but another database names HGNC:2"
        ));

        assert_eq!(merged_hgnc_ids(None, Some("HGNC:3")), [None]);
        assert!(logs_contain(
            "gene 11 has HGNC ID none, but another database names HGNC:3"
        ));
    }

    /// The merge does not warn if the other database names the same HGNC ID or none.
    #[tracing_test::traced_test]
    #[test]
    fn merge_accepts_the_same_or_no_hgnc_id() {
        for other in [Some("HGNC:1"), None] {
            assert_eq!(
                merged_hgnc_ids(Some("HGNC:1"), other),
                [Some("HGNC:1".to_string())]
            );
        }
        assert!(!logs_contain("HGNC ID"));
    }

    /// Two genes can share one `legacy_gene_id`. The merge keeps them apart by `gene_id`.
    #[test]
    fn genes_that_share_a_legacy_gene_id_stay_distinct() {
        let first = db(&[("11", "HGNC:1", "NM_000001.1")]);
        let other = db(&[
            ("11", "HGNC:1", "NM_000002.1"),
            ("22", "HGNC:1", "NM_000003.1"),
        ]);

        assert_eq!(
            merged_genes(vec![first, other]),
            [
                gene("11", "HGNC:1", &["NM_000001.1", "NM_000002.1"]),
                gene("22", "HGNC:1", &["NM_000003.1"]),
            ]
        );
    }
}
