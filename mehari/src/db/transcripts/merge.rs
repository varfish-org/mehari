use crate::annotate::seqvars::consequence::load_tx_db;
use crate::common;
use crate::db::transcripts::create::write_tx_db;
use crate::pbs::txs::TxSeqDatabase;
use clap::Parser;
use itertools::Itertools;
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};

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

        // A transcript can be in several databases, e.g., the Ensembl transcripts that RefSeq
        // databases contain for chrMT. The copy of the first database that has it wins, and
        // the others are skipped with their sequences, aliases and gene memberships.
        let mut known_tx_ids: HashSet<_> =
            tx_db.transcripts.iter().map(|tx| tx.id.clone()).collect();
        let mut skipped_count = 0;
        for other in others.iter_mut() {
            first.source_version.append(&mut other.source_version);

            let mut other_tx_db = other.tx_db.take().unwrap();
            let skipped: HashSet<_> = (other_tx_db.transcripts.iter())
                .filter(|tx| known_tx_ids.contains(&tx.id))
                .map(|tx| tx.id.clone())
                .collect();
            skipped_count += skipped.len();
            other_tx_db
                .transcripts
                .retain(|tx| !skipped.contains(&tx.id));
            known_tx_ids.extend(other_tx_db.transcripts.iter().map(|tx| tx.id.clone()));

            // Merge the sequence database, renumbering the sequences that are kept.
            let mut other_seq_db = other.seq_db.take().unwrap();
            let mut seq_idx = HashMap::new();
            for (alias, idx) in other_seq_db
                .aliases
                .into_iter()
                .zip(other_seq_db.aliases_idx)
            {
                if skipped.contains(&alias) {
                    continue;
                }
                let idx = match seq_idx.entry(idx) {
                    Entry::Occupied(entry) => *entry.get(),
                    Entry::Vacant(entry) => {
                        let seq = other_seq_db
                            .seqs
                            .get_mut(idx as usize)
                            .ok_or_else(|| anyhow::anyhow!("No sequence for alias {alias}"))?;
                        seq_db.seqs.push(std::mem::take(seq));
                        *entry.insert(seq_db.seqs.len() as u32 - 1)
                    }
                };
                seq_db.aliases.push(alias);
                seq_db.aliases_idx.push(idx);
            }

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
                .filter_map(|mut g2t| {
                    // Drop a gene whose transcripts were all skipped.
                    let tx_count = g2t.tx_ids.len();
                    g2t.tx_ids.retain(|tx_id| !skipped.contains(tx_id));
                    (tx_count == 0 || !g2t.tx_ids.is_empty()).then(|| (g2t.gene_id.clone(), g2t))
                })
                .collect();
            for g2tx in tx_db.gene_to_tx.iter_mut() {
                if let Some(other_g2tx) = other_gene_to_tx.remove(&g2tx.gene_id) {
                    // Keep the HGNC ID of the first database, else take the other one. Warn as
                    // the importer does if the two differ.
                    if g2tx.hgnc_id.is_none() {
                        g2tx.hgnc_id = other_g2tx.hgnc_id;
                    } else if let (Some(hgnc_id), Some(other_hgnc_id)) =
                        (&g2tx.hgnc_id, &other_g2tx.hgnc_id)
                        && hgnc_id != other_hgnc_id
                    {
                        tracing::warn!(
                            "gene {} has HGNC ID {hgnc_id}, but another database names \
                             {other_hgnc_id}",
                            g2tx.gene_id
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
        if skipped_count > 0 {
            tracing::info!("Skipped {skipped_count} transcripts that an earlier database has");
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
    use crate::pbs::txs::{
        GeneToTxId, SequenceDb, SourceVersion, Transcript, TranscriptDb, TranscriptTag,
        TxSeqDatabase,
    };

    /// A database with one `GeneToTxId` for each `(gene_id, hgnc_id, tx_id)`.
    fn db(genes: &[(&str, &str, &str)]) -> TxSeqDatabase {
        TxSeqDatabase {
            tx_db: Some(TranscriptDb {
                gene_to_tx: genes
                    .iter()
                    .map(|&(gene_id, hgnc_id, tx_id)| GeneToTxId {
                        gene_id: gene_id.into(),
                        hgnc_id: Some(hgnc_id.into()),
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

    /// The `(gene_id, hgnc_id, tx_ids)` of the merged database, sorted.
    fn merged_genes(databases: Vec<TxSeqDatabase>) -> Vec<(String, Option<String>, Vec<String>)> {
        let mut genes: Vec<_> = merge_transcript_databases(databases)
            .unwrap()
            .tx_db
            .unwrap()
            .gene_to_tx
            .into_iter()
            .map(|gene| (gene.gene_id, gene.hgnc_id, gene.tx_ids))
            .collect();
        genes.sort();
        genes
    }

    fn gene(
        gene_id: &str,
        hgnc_id: &str,
        tx_ids: &[&str],
    ) -> (String, Option<String>, Vec<String>) {
        (
            gene_id.into(),
            Some(hgnc_id.into()),
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
    /// database names another one.
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
    }

    /// A gene without HGNC ID in the first database takes the one of the other database.
    #[tracing_test::traced_test]
    #[test]
    fn merge_takes_the_hgnc_id_of_the_other_database() {
        assert_eq!(
            merged_hgnc_ids(None, Some("HGNC:3")),
            [Some("HGNC:3".to_string())]
        );
        assert!(!logs_contain("HGNC ID"));
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

    /// A database of `source_name` with one gene, transcript and sequence for each
    /// `(gene_id, tx_id, seq)`.
    fn db_with_txs(source_name: &str, txs: &[(&str, &str, &str)]) -> TxSeqDatabase {
        TxSeqDatabase {
            tx_db: Some(TranscriptDb {
                transcripts: txs
                    .iter()
                    .map(|&(gene_id, tx_id, _)| Transcript {
                        id: tx_id.into(),
                        gene_id: gene_id.into(),
                        ..Default::default()
                    })
                    .collect(),
                gene_to_tx: txs
                    .iter()
                    .map(|&(gene_id, tx_id, _)| GeneToTxId {
                        gene_id: gene_id.into(),
                        tx_ids: vec![tx_id.into()],
                        ..Default::default()
                    })
                    .collect(),
            }),
            seq_db: Some(SequenceDb {
                aliases: txs.iter().map(|&(_, tx_id, _)| tx_id.into()).collect(),
                aliases_idx: (0..txs.len() as u32).collect(),
                seqs: txs.iter().map(|&(_, _, seq)| seq.into()).collect(),
            }),
            version: Some("test".into()),
            source_version: vec![SourceVersion {
                assembly: "GRCh38".into(),
                source_name: source_name.into(),
                ..Default::default()
            }],
        }
    }

    /// A transcript whose ID an earlier database has is skipped, with its sequence and gene
    /// membership. Here, the RefSeq database has a copy of an Ensembl transcript.
    #[tracing_test::traced_test]
    #[test]
    fn merge_skips_transcripts_of_earlier_databases() {
        let first = db_with_txs("ensembl", &[("ENSG1", "ENST1.1", "ACGT")]);
        let mut other = db_with_txs(
            "refseq",
            &[("ENSG1", "ENST1.1", "TTTT"), ("2", "NM_000002.1", "GG")],
        );
        other.tx_db.as_mut().unwrap().transcripts[0].tags =
            vec![TranscriptTag::EnsemblGraft.into()];

        let merged = merge_transcript_databases(vec![first, other]).unwrap();

        let tx_db = merged.tx_db.unwrap();
        let transcripts: Vec<_> = (tx_db.transcripts.iter())
            .map(|tx| (tx.id.as_str(), tx.tags.as_slice()))
            .collect();
        assert_eq!(
            transcripts,
            [("ENST1.1", &[][..]), ("NM_000002.1", &[][..])]
        );
        let genes: Vec<_> = (tx_db.gene_to_tx.iter())
            .map(|gene| (gene.gene_id.as_str(), gene.tx_ids.clone()))
            .collect();
        assert_eq!(
            genes,
            [
                ("ENSG1", vec!["ENST1.1".to_string()]),
                ("2", vec!["NM_000002.1".to_string()])
            ]
        );
        let seq_db = merged.seq_db.unwrap();
        assert_eq!(seq_db.aliases, ["ENST1.1", "NM_000002.1"]);
        assert_eq!(seq_db.aliases_idx, [0, 1]);
        assert_eq!(seq_db.seqs, ["ACGT", "GG"]);
        assert!(logs_contain("Skipped 1 transcripts"));
    }

    /// A skipped transcript can be in a gene of another ID, e.g., of an NCBI Gene ID in a RefSeq
    /// database. It leaves that gene, and a gene without other transcripts is dropped. A gene
    /// that had no transcripts stays.
    #[test]
    fn merge_drops_genes_left_without_transcripts() {
        let first = db_with_txs("ensembl", &[("ENSG1", "ENST1.1", "ACGT")]);
        let mut other = db_with_txs(
            "refseq",
            &[("4567", "ENST1.1", "ACGT"), ("2", "NM_000002.1", "GG")],
        );
        (other.tx_db.as_mut().unwrap().gene_to_tx).push(GeneToTxId {
            gene_id: "3".into(),
            ..Default::default()
        });

        let mut genes: Vec<_> = merge_transcript_databases(vec![first, other])
            .unwrap()
            .tx_db
            .unwrap()
            .gene_to_tx
            .into_iter()
            .map(|gene| (gene.gene_id, gene.tx_ids))
            .collect();
        genes.sort();

        assert_eq!(
            genes,
            [
                ("2".to_string(), vec!["NM_000002.1".to_string()]),
                ("3".to_string(), vec![]),
                ("ENSG1".to_string(), vec!["ENST1.1".to_string()]),
            ]
        );
    }

    /// The merged database keeps the source versions of all databases.
    #[test]
    fn merge_keeps_all_source_versions() {
        let databases = vec![
            db_with_txs("ensembl", &[("ENSG1", "ENST1.1", "ACGT")]),
            db_with_txs("refseq", &[("2", "NM_000002.1", "GG")]),
        ];

        let merged = merge_transcript_databases(databases).unwrap();

        let source_names: Vec<_> = (merged.source_version.iter())
            .map(|sv| sv.source_name.as_str())
            .collect();
        assert_eq!(source_names, ["ensembl", "refseq"]);
    }

    /// Two genes can share one HGNC ID. The merge keeps them apart by `gene_id`.
    #[test]
    fn genes_that_share_an_hgnc_id_stay_distinct() {
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
