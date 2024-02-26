//! Subset transcript database.

use std::{io::Write, path::PathBuf};

use clap::Parser;
use indexmap::{IndexMap, IndexSet};
use prost::Message as _;

use crate::annotate::seqvars::load_tx_db;

/// Command line arguments for `db subset` sub command.
#[derive(Parser, Debug)]
#[command(about = "Subset transcript database", long_about = None)]
pub struct Args {
    /// Path to database file to read from.
    #[arg(long)]
    pub path_in: PathBuf,
    // Path to output file to write to.
    #[arg(long)]
    pub path_out: PathBuf,
    /// Limit transcript database to the following HGNC symbols.
    #[arg(long)]
    pub gene_symbols: Vec<String>,
}

fn subset_tx_db(
    container: &crate::pbs::txs::TxSeqDatabase,
    gene_symbols: &[String],
) -> Result<crate::pbs::txs::TxSeqDatabase, anyhow::Error> {
    let container_tx_db = container.tx_db.as_ref().expect("no tx_db");
    let (tx_idxs, tx_ids, gene_ids) = {
        let gene_symbols = IndexSet::<String>::from_iter(gene_symbols.iter().cloned());
        let mut tx_idxs = Vec::new();
        let mut tx_ids = Vec::new();
        let mut gene_ids = IndexSet::new();
        for (idx, tx) in container_tx_db.transcripts.iter().enumerate() {
            tracing::debug!("considering transcript: {:?}", &tx);
            if gene_symbols.contains(&tx.gene_symbol)
                || gene_symbols.contains(&tx.gene_id)
                || gene_symbols.contains(&tx.gene_id.replace("HGNC:", ""))
            {
                tracing::debug!("keeping transcript: {}", tx.id);
                tx_idxs.push(idx);
                tx_ids.push(tx.id.clone());
                gene_ids.insert(tx.gene_id.clone());
            }
        }
        (
            IndexSet::<usize>::from_iter(tx_idxs.into_iter()),
            IndexSet::<String>::from_iter(tx_ids.into_iter()),
            gene_ids,
        )
    };

    let tx_db = Some(crate::pbs::txs::TranscriptDb {
        transcripts: container_tx_db
            .transcripts
            .iter()
            .enumerate()
            .filter_map(|(idx, tx)| {
                if tx_idxs.contains(&idx) {
                    Some(tx.clone())
                } else {
                    None
                }
            })
            .collect(),
        gene_to_tx: container_tx_db
            .gene_to_tx
            .iter()
            .filter_map(|gene_to_tx| {
                if gene_ids.contains(gene_to_tx.gene_id.as_str()) {
                    Some(gene_to_tx.clone())
                } else {
                    None
                }
            })
            .collect(),
    });

    let container_seq_db = container.seq_db.as_ref().expect("no seq_db");
    let seq_db = {
        let mut old_to_new_idx = IndexMap::<u32, u32>::new();

        let mut aliases = Vec::new();
        let mut aliases_idx = Vec::new();
        let mut seqs = Vec::new();

        for (alias, old_alias_idx) in container_seq_db
            .aliases
            .iter()
            .zip(container_seq_db.aliases_idx.iter())
        {
            if tx_ids.contains(alias) {
                if let Some(new_alias_idx) = old_to_new_idx.get(old_alias_idx) {
                    aliases.push(alias.clone());
                    aliases_idx.push(*new_alias_idx);
                } else {
                    let new_alias_idx = seqs.len();
                    old_to_new_idx.insert(*old_alias_idx, new_alias_idx as u32);
                    aliases.push(alias.clone());
                    aliases_idx.push(new_alias_idx as u32);
                    seqs.push(container_seq_db.seqs[*old_alias_idx as usize].clone());
                }
            }
        }

        Some(crate::pbs::txs::SequenceDb {
            aliases,
            aliases_idx,
            seqs,
        })
    };

    Ok(crate::pbs::txs::TxSeqDatabase {
        tx_db,
        seq_db,
        version: container.version.clone(),
        genome_release: container.genome_release.clone(),
    })
}

/// Main entry point.
pub fn run(_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("Opening transcript database");
    let tx_db = load_tx_db(&format!("{}", args.path_in.display()))?;

    tracing::info!("Subsetting ...");
    let tx_db = subset_tx_db(&tx_db, &args.gene_symbols)?;

    tracing::info!("... and encoding ...");
    let mut buf = Vec::with_capacity(tx_db.encoded_len());
    tx_db
        .encode(&mut buf)
        .map_err(|e| anyhow::anyhow!("failed to encode: {}", e))?;

    tracing::info!("... and writing ...");
    // Open file and if necessary, wrap in a decompressor.
    let file = std::fs::File::create(&args.path_out)
        .map_err(|e| anyhow::anyhow!("failed to create file {}: {}", args.path_out.display(), e))?;
    let ext = args.path_out.extension().map(|s| s.to_str());
    let mut writer: Box<dyn Write> = if ext == Some(Some("gz")) {
        Box::new(flate2::write::GzEncoder::new(
            file,
            flate2::Compression::default(),
        ))
    } else if ext == Some(Some("zst")) {
        Box::new(
            zstd::Encoder::new(file, 0)
                .map_err(|e| {
                    anyhow::anyhow!(
                        "failed to open zstd encoder for {}: {}",
                        args.path_out.display(),
                        e
                    )
                })?
                .auto_finish(),
        )
    } else {
        Box::new(file)
    };
    writer
        .write_all(&buf)
        .map_err(|e| anyhow::anyhow!("failed to write to {}: {}", args.path_out.display(), e))?;

    tracing::info!("... done");
    Ok(())
}

#[cfg(test)]
mod tests {
    use temp_testdir::TempDir;

    #[test]
    fn test_subset_tx_db() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("txs.bin.zst");

        super::run(
            &crate::common::Args::default(),
            &super::Args {
                path_in: "tests/data/annotate/db/grch37/txs.bin.zst".into(),
                path_out: path_out.clone(),
                gene_symbols: vec!["HGNC:12403".into()],
            },
        )?;

        let mut buf = Vec::new();
        super::super::dump::run_with_write(
            &crate::common::Args::default(),
            &super::super::dump::Args { path_db: path_out },
            &mut buf,
        )?;

        insta::assert_display_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }
}
