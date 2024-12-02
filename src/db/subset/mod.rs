//! Subset transcript database.

use crate::annotate::seqvars::load_tx_db;
use crate::annotate::seqvars::provider::TxIntervalTrees;
use crate::db::TranscriptDatabase;
use crate::pbs::txs::{TranscriptDb, TxSeqDatabase};
use anyhow::{Error, Result};
use biocommons_bioutils::assemblies::{Assembly, ASSEMBLY_INFOS};
use clap::Parser;
use indexmap::{IndexMap, IndexSet};
use noodles::vcf::variant::Record;
use prost::Message as _;
use std::{io::Write, path::PathBuf};

/// Command line arguments for `db subset` sub command.
#[derive(Parser, Debug)]
#[command(about = "Subset transcript database", long_about = None)]
pub struct Args {
    /// Path to database file to read from.
    #[arg(long)]
    pub path_in: PathBuf,

    /// Path to output file to write to.
    #[arg(long)]
    pub path_out: PathBuf,

    #[command(flatten)]
    pub selection: Selection,
}

#[derive(Debug, Default, Clone, clap::Args)]
#[group(required = true, multiple = false)]
pub struct Selection {
    /// Limit transcript database to the transcripts affected by the variants described in the specified VCF file.
    #[arg(long)]
    vcf: Option<PathBuf>,

    /// Limit transcript database to the specified HGNC ID. Can be specified multiple times.
    #[arg(long)]
    hgnc_id: Option<Vec<String>>,

    /// Limit transcript database to the specified transcript ID. Can be specified multiple times.
    #[arg(long)]
    transcript_id: Option<Vec<String>>,
}

/// Main entry point.
pub fn run(_common: &crate::common::Args, args: &Args) -> Result<(), Error> {
    tracing::info!("Opening transcript database");
    let tx_db = load_tx_db(format!("{}", args.path_in.display()))?;

    tracing::info!("Subsetting ...");
    let tx_db = subset_tx_db(&tx_db, &args.selection)?;

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

fn subset_tx_db(container: &TxSeqDatabase, selection: &Selection) -> Result<TxSeqDatabase> {
    let container_tx_db = container.tx_db.as_ref().expect("no tx_db");

    let Selection {
        vcf,
        hgnc_id: hgnc_ids,
        transcript_id: transcript_ids,
    } = selection;

    let (tx_idxs, tx_ids, gene_ids) = match (hgnc_ids, transcript_ids, vcf) {
        (Some(hgnc_ids), _, _) => _extract_transcripts_by_hgnc_id(container_tx_db, hgnc_ids)?,
        (_, Some(transcript_ids), _) => {
            _extract_transcripts_by_txid(container_tx_db, transcript_ids)
        }
        (_, _, Some(vcf)) => _extract_transcripts_by_region(container, vcf)?,
        _ => unreachable!(),
    };
    tracing::info!(
        "Keeping {} transcripts (of {} gene(s))",
        tx_idxs.len(),
        gene_ids.len()
    );

    let tx_db = Some(TranscriptDb {
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

    Ok(TxSeqDatabase {
        tx_db,
        seq_db,
        version: container.version.clone(),
        source_version: container.source_version.clone(),
    })
}

fn _extract_transcripts_by_region(
    container: &TxSeqDatabase,
    vcf: &PathBuf,
) -> Result<(IndexSet<usize>, IndexSet<String>, IndexSet<String>)> {
    let container_tx_db = container.tx_db.as_ref().expect("no tx_db");
    let trees = TxIntervalTrees::new(container);

    let chrom_to_acc = _chrom_to_acc(container.assembly());

    let mut transcript_ids = IndexSet::<String>::new();
    let mut reader = noodles::vcf::io::reader::Builder::default().build_from_path(vcf)?;
    let header = reader.read_header()?;
    for result in reader.records() {
        let record = result?;
        if let (Some(Ok(start)), Ok(end)) = (record.variant_start(), record.variant_end(&header)) {
            let chrom = record.reference_sequence_name().to_string();
            let accession = chrom_to_acc.get(&chrom).unwrap_or(&chrom);

            let txs = trees.get_tx_for_region(
                container,
                accession,
                "splign",
                usize::from(start).try_into()?,
                usize::from(end).try_into()?,
            )?;
            transcript_ids.extend(txs.into_iter().map(|tx| tx.tx_ac));
        }
    }
    if transcript_ids.is_empty() {
        return Err(anyhow::anyhow!(
            "no transcripts overlapping VCF variants found in transcript database"
        ));
    }
    let (tx_idxs, gene_ids) = __extract_transcripts_from_db(container_tx_db, &transcript_ids);
    Ok((
        IndexSet::<usize>::from_iter(tx_idxs),
        transcript_ids,
        gene_ids,
    ))
}

fn _extract_transcripts_by_txid(
    container_tx_db: &TranscriptDb,
    transcript_ids: &Vec<String>,
) -> (IndexSet<usize>, IndexSet<String>, IndexSet<String>) {
    let transcript_ids = IndexSet::<String>::from_iter(transcript_ids.iter().cloned());
    let (tx_idxs, gene_ids) = __extract_transcripts_from_db(container_tx_db, &transcript_ids);
    (
        IndexSet::<usize>::from_iter(tx_idxs),
        transcript_ids,
        gene_ids,
    )
}

fn _extract_transcripts_by_hgnc_id(
    container_tx_db: &TranscriptDb,
    gene_symbols: &Vec<String>,
) -> Result<(IndexSet<usize>, IndexSet<String>, IndexSet<String>)> {
    let hgnc_ids = IndexSet::<String>::from_iter(gene_symbols.iter().cloned());
    let mut tx_idxs = Vec::new();
    let mut tx_ids = Vec::new();
    let mut gene_ids = IndexSet::new();
    for (idx, tx) in container_tx_db.transcripts.iter().enumerate() {
        tracing::debug!("considering transcript: {:?}", &tx);
        if hgnc_ids.contains(&tx.gene_symbol)
            || hgnc_ids.contains(&tx.gene_id)
            || hgnc_ids.contains(&tx.gene_id.replace("HGNC:", ""))
        {
            tracing::debug!("keeping transcript: {}", tx.id);
            tx_idxs.push(idx);
            tx_ids.push(tx.id.clone());
            gene_ids.insert(tx.gene_id.clone());
        }
    }
    if tx_ids.is_empty() {
        return Err(anyhow::anyhow!(
            "no transcripts found for HGNC IDs: {:?}",
            hgnc_ids
        ));
    }
    Ok((
        IndexSet::<usize>::from_iter(tx_idxs),
        IndexSet::<String>::from_iter(tx_ids),
        gene_ids,
    ))
}

fn __extract_transcripts_from_db(
    container_tx_db: &TranscriptDb,
    transcript_ids: &IndexSet<String>,
) -> (Vec<usize>, IndexSet<String>) {
    let mut tx_idxs = Vec::new();
    let mut gene_ids = IndexSet::new();
    for (idx, tx) in container_tx_db.transcripts.iter().enumerate() {
        tracing::debug!("considering transcript: {:?}", &tx);
        if transcript_ids.contains(&tx.id) {
            tracing::debug!("keeping transcript: {}", tx.id);
            tx_idxs.push(idx);
            gene_ids.insert(tx.gene_id.clone());
        }
    }
    (tx_idxs, gene_ids)
}

fn _acc_to_chrom(assembly: Assembly) -> IndexMap<String, String> {
    IndexMap::from_iter(
        ASSEMBLY_INFOS[assembly]
            .sequences
            .iter()
            .map(|record| (record.refseq_ac.clone(), record.name.clone())),
    )
}

fn _chrom_to_acc(assembly: Assembly) -> IndexMap<String, String> {
    let acc_to_chrom = _acc_to_chrom(assembly);
    let mut chrom_to_acc = IndexMap::new();
    for (acc, chrom) in &acc_to_chrom {
        let chrom = if chrom.starts_with("chr") {
            chrom.strip_prefix("chr").unwrap()
        } else {
            chrom
        };
        chrom_to_acc.insert(chrom.to_string(), acc.clone());
        chrom_to_acc.insert(format!("chr{}", chrom), acc.clone());
    }
    chrom_to_acc
}

#[cfg(test)]
mod tests {
    use crate::db::subset::Selection;
    use temp_testdir::TempDir;

    #[tracing_test::traced_test]
    #[test]
    fn test_subset_tx_db() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("txs.bin.zst");

        super::run(
            &crate::common::Args::default(),
            &super::Args {
                path_in: "tests/data/annotate/db/grch37/txs.bin.zst".into(),
                path_out: path_out.clone(),
                selection: Selection {
                    hgnc_id: Some(vec!["HGNC:12403".into()]),
                    ..Default::default()
                },
            },
        )?;

        let mut buf = Vec::new();
        super::super::dump::run_with_write(
            &crate::common::Args::default(),
            &super::super::dump::Args { path_db: path_out },
            &mut buf,
        )?;

        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }
}
