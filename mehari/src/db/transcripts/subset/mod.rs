//! Subset transcript database.

use crate::annotate::seqvars::consequence::load_tx_db;
use crate::annotate::seqvars::provider::TxIntervalTrees;
use crate::common::contig::ContigManager;
use crate::db::TranscriptDatabase;
use crate::pbs::txs::{TranscriptDb, TxSeqDatabase};
use anyhow::{Error, Result};
use clap::Parser;
use indexmap::{IndexMap, IndexSet};
use noodles::core::Region;
use noodles::vcf::variant::Record;
use prost::Message as _;
use std::collections::Bound;
use std::{io::Write, path::PathBuf};

/// Command line arguments for `db subset` sub command.
///
/// The subsetting process works in two stages:
/// 1. First, an initial set of transcripts is selected based on zero or one criterion
///    provided in the optional `Selection` group (e.g. by HGNC ID, transcript ID, region, or VCF).
/// 2. Then, the resulting transcripts can be further filtered using the optional
///    `--include-transcripts` and `--exclude-transcripts` regular expressions. If an
///    include regex is provided, only transcripts matching it are kept. If an exclude
///    regex is provided, any matching transcripts are discarded. If a transcript matches
///    both the include and exclude regex, an error is returned.
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

    /// Compression level to use when writing the output file.
    #[arg(long)]
    pub compression_level: Option<u32>,

    /// Regular expression to include transcripts by their ID. Only transcripts matching this
    /// expression will be kept in the final subset. Applied after the initial selection.
    #[arg(long)]
    pub include_transcripts: Option<String>,

    /// Regular expression to exclude transcripts by their ID. Transcripts matching this
    /// expression will be dropped from the final subset. Applied after the initial selection.
    #[arg(long)]
    pub exclude_transcripts: Option<String>,
}

#[derive(Debug, Default, Clone, clap::Args)]
#[group(required = false, multiple = false)]
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

    /// Limit transcript database to the transcripts overlapping the specified region. Can be specified multiple times.
    #[arg(long)]
    region: Option<Vec<Region>>,
}

/// Main entry point.
pub fn run(_common: &crate::common::Args, args: &Args) -> Result<(), Error> {
    tracing::info!("Opening transcript database");
    let tx_db = load_tx_db(format!("{}", args.path_in.display()))?;

    tracing::info!("Subsetting ...");
    let tx_db = subset_tx_db(&tx_db, args)?;

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
            args.compression_level
                .map(flate2::Compression::new)
                .unwrap_or_default(),
        ))
    } else if ext == Some(Some("zst")) {
        let mut encoder =
            zstd::Encoder::new(file, args.compression_level.map(|l| l as i32).unwrap_or(0))
                .map_err(|e| {
                    anyhow::anyhow!(
                        "failed to open zstd encoder for {}: {}",
                        args.path_out.display(),
                        e
                    )
                })?;
        // Record the decompressed size in the frame header, so `load_tx_db` can allocate it up front.
        encoder
            .set_pledged_src_size(Some(buf.len() as u64))
            .map_err(|e| {
                anyhow::anyhow!(
                    "failed to set zstd content size for {}: {}",
                    args.path_out.display(),
                    e
                )
            })?;
        Box::new(encoder.auto_finish())
    } else {
        Box::new(file)
    };
    writer
        .write_all(&buf)
        .map_err(|e| anyhow::anyhow!("failed to write to {}: {}", args.path_out.display(), e))?;

    tracing::info!("... done");
    Ok(())
}

#[allow(deprecated)]
fn subset_tx_db(container: &TxSeqDatabase, args: &Args) -> Result<TxSeqDatabase> {
    let container_tx_db = container.tx_db.as_ref().expect("no tx_db");

    let Selection {
        vcf,
        hgnc_id: hgnc_ids,
        transcript_id: transcript_ids,
        region: regions,
    } = &args.selection;

    let include_re = match &args.include_transcripts {
        Some(re) => Some(regex::Regex::new(re)?),
        None => None,
    };
    let exclude_re = match &args.exclude_transcripts {
        Some(re) => Some(regex::Regex::new(re)?),
        None => None,
    };

    let (mut tx_idxs, mut tx_ids, mut gene_ids) = match (hgnc_ids, transcript_ids, vcf, regions) {
        (Some(hgnc_ids), _, _, _) => _extract_transcripts_by_hgnc_id(container_tx_db, hgnc_ids)?,
        (_, Some(transcript_ids), _, _) => {
            _extract_transcripts_by_txid(container_tx_db, transcript_ids)
        }
        (_, _, Some(vcf), _) => _extract_transcripts_by_vcf(container, vcf)?,
        (_, _, _, Some(regions)) => _extract_transcripts_by_regions(container, regions)?,
        (None, None, None, None) => {
            let mut tx_idxs = IndexSet::new();
            let mut tx_ids = IndexSet::new();
            let mut gene_ids = IndexSet::new();
            for (idx, tx) in container_tx_db.transcripts.iter().enumerate() {
                tx_idxs.insert(idx);
                tx_ids.insert(tx.id.clone());
                gene_ids.insert(tx.legacy_gene_id.clone());
            }
            (tx_idxs, tx_ids, gene_ids)
        }
    };

    if include_re.is_some() || exclude_re.is_some() {
        let mut new_tx_idxs = IndexSet::new();
        let mut new_tx_ids = IndexSet::new();
        let mut new_gene_ids = IndexSet::new();

        for (idx, tx) in container_tx_db.transcripts.iter().enumerate() {
            if tx_idxs.contains(&idx) {
                let tx_id = &tx.id;

                let included = include_re.as_ref().map(|re| re.is_match(tx_id));
                let excluded = exclude_re.as_ref().map(|re| re.is_match(tx_id));

                if included == Some(true) && excluded == Some(true) {
                    anyhow::bail!(
                        "Transcript {} matches both include and exclude regular expressions. Please adjust the regular expressions to not overlap.",
                        tx_id
                    );
                }

                let keep = included.unwrap_or(true) && !excluded.unwrap_or(false);

                if keep {
                    new_tx_idxs.insert(idx);
                    new_tx_ids.insert(tx_id.clone());
                    new_gene_ids.insert(tx.legacy_gene_id.clone());
                }
            }
        }

        tx_idxs = new_tx_idxs;
        tx_ids = new_tx_ids;
        gene_ids = new_gene_ids;
    }

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
                let mut mapping = gene_to_tx.clone();
                mapping.tx_ids.retain(|tx_id| tx_ids.contains(tx_id));
                if mapping.tx_ids.is_empty() {
                    None
                } else {
                    Some(mapping)
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
                match old_to_new_idx.get(old_alias_idx) {
                    Some(new_alias_idx) => {
                        aliases.push(alias.clone());
                        aliases_idx.push(*new_alias_idx);
                    }
                    _ => {
                        let new_alias_idx = seqs.len();
                        old_to_new_idx.insert(*old_alias_idx, new_alias_idx as u32);
                        aliases.push(alias.clone());
                        aliases_idx.push(new_alias_idx as u32);
                        seqs.push(container_seq_db.seqs[*old_alias_idx as usize].clone());
                    }
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

fn _extract_transcripts_by_vcf(
    container: &TxSeqDatabase,
    vcf: &PathBuf,
) -> Result<(IndexSet<usize>, IndexSet<String>, IndexSet<String>)> {
    let container_tx_db = container.tx_db.as_ref().expect("no tx_db");
    let trees = TxIntervalTrees::new(container);

    let contig_manager = ContigManager::new(&container.assembly());

    let mut transcript_ids = IndexSet::<String>::new();
    let mut reader = noodles::vcf::io::reader::Builder::default().build_from_path(vcf)?;
    let header = reader.read_header()?;
    for result in reader.records() {
        let record = result?;
        if let (Some(Ok(start)), Ok(end)) = (record.variant_start(), record.variant_end(&header)) {
            let chrom = record.reference_sequence_name().to_string();
            if let Some(accession) = contig_manager.get_accession(&chrom) {
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

fn _extract_transcripts_by_regions(
    container: &TxSeqDatabase,
    regions: &[Region],
) -> Result<(IndexSet<usize>, IndexSet<String>, IndexSet<String>)> {
    let container_tx_db = container.tx_db.as_ref().expect("no tx_db");
    let trees = TxIntervalTrees::new(container);

    let contig_manager = ContigManager::new(&container.assembly());

    let mut transcript_ids = IndexSet::<String>::new();
    for region in regions {
        let (start, end) = (region.start(), region.end());
        let chrom = region.name().to_string();
        if let Some(accession) = contig_manager.get_accession(&chrom) {
            let txs = trees.get_tx_for_region(
                container,
                accession,
                "splign",
                match start {
                    Bound::Included(i) => usize::from(i).try_into()?,
                    Bound::Excluded(i) => usize::from(i).try_into().map(|i: i32| i + 1)?,
                    Bound::Unbounded => 0,
                },
                match end {
                    Bound::Included(i) => usize::from(i).try_into()?,
                    Bound::Excluded(i) => usize::from(i).try_into().map(|i: i32| i - 1)?,
                    Bound::Unbounded => i32::MAX,
                },
            )?;
            transcript_ids.extend(txs.into_iter().map(|tx| tx.tx_ac));
        }
    }
    if transcript_ids.is_empty() {
        return Err(anyhow::anyhow!(
            "no transcripts overlapping the specified regions found in transcript database"
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
    transcript_ids: &[String],
) -> (IndexSet<usize>, IndexSet<String>, IndexSet<String>) {
    let transcript_ids = IndexSet::<String>::from_iter(transcript_ids.iter().cloned());
    let (tx_idxs, gene_ids) = __extract_transcripts_from_db(container_tx_db, &transcript_ids);
    (
        IndexSet::<usize>::from_iter(tx_idxs),
        transcript_ids,
        gene_ids,
    )
}

#[allow(deprecated)]
fn _extract_transcripts_by_hgnc_id(
    container_tx_db: &TranscriptDb,
    gene_symbols: &[String],
) -> Result<(IndexSet<usize>, IndexSet<String>, IndexSet<String>)> {
    let hgnc_ids = IndexSet::<String>::from_iter(gene_symbols.iter().cloned());
    let mut tx_idxs = Vec::new();
    let mut tx_ids = Vec::new();
    let mut gene_ids = IndexSet::new();
    for (idx, tx) in container_tx_db.transcripts.iter().enumerate() {
        tracing::debug!("considering transcript: {:?}", &tx);
        if hgnc_ids.contains(&tx.gene_symbol)
            || tx.hgnc_id_or_legacy().is_some_and(|hgnc_id| {
                hgnc_ids.contains(hgnc_id)
                    || hgnc_id
                        .strip_prefix("HGNC:")
                        .is_some_and(|number| hgnc_ids.contains(number))
            })
        {
            tracing::debug!("keeping transcript: {}", tx.id);
            tx_idxs.push(idx);
            tx_ids.push(tx.id.clone());
            gene_ids.insert(tx.legacy_gene_id.clone());
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

#[allow(deprecated)]
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
            gene_ids.insert(tx.legacy_gene_id.clone());
        }
    }
    (tx_idxs, gene_ids)
}

#[cfg(test)]
mod tests {
    use crate::annotate::seqvars::consequence::load_tx_db;
    use crate::db::transcripts::subset::Selection;
    use crate::pbs::txs::{GeneToTxId, SequenceDb, Transcript, TranscriptDb, TxSeqDatabase};
    use prost::Message;
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
                compression_level: None,
                include_transcripts: None,
                exclude_transcripts: None,
            },
        )?;

        let tx_db = load_tx_db(format!("{}", path_out.display()))?;
        let compressed = std::fs::read(&path_out)?;
        assert_eq!(
            zstd::zstd_safe::get_frame_content_size(&compressed)
                .ok()
                .flatten(),
            Some(tx_db.encoded_len() as u64)
        );

        let mut buf = Vec::new();
        super::super::dump::run_with_write(
            &crate::common::Args::default(),
            &super::super::dump::Args { path_db: path_out },
            &mut buf,
        )?;

        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[test]
    fn test_subset_tx_db_with_regex_only() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("txs_regex_only.bin.zst");

        super::run(
            &crate::common::Args::default(),
            &super::Args {
                path_in: "tests/data/annotate/db/grch37/txs.bin.zst".into(),
                path_out: path_out.clone(),
                selection: Selection::default(),
                compression_level: None,
                include_transcripts: Some("NM_.*".into()),
                exclude_transcripts: None,
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

    #[tracing_test::traced_test]
    #[test]
    fn test_subset_tx_db_with_regex_filters() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("txs_filtered.bin.zst");

        super::run(
            &crate::common::Args::default(),
            &super::Args {
                path_in: "tests/data/annotate/db/grch37/txs.bin.zst".into(),
                path_out: path_out.clone(),
                selection: Selection {
                    hgnc_id: Some(vec!["HGNC:12403".into()]),
                    ..Default::default()
                },
                compression_level: None,
                // Only include NM_ and exclude a specific one if it existed (we'll just use a pattern)
                include_transcripts: Some("NM_.*".into()),
                exclude_transcripts: Some(".*_999999.*".into()),
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

    #[tracing_test::traced_test]
    #[test]
    fn test_subset_tx_db_with_regex_contradiction() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("txs_err.bin.zst");

        let result = super::run(
            &crate::common::Args::default(),
            &super::Args {
                path_in: "tests/data/annotate/db/grch37/txs.bin.zst".into(),
                path_out: path_out.clone(),
                selection: Selection {
                    hgnc_id: Some(vec!["HGNC:12403".into()]),
                    ..Default::default()
                },
                compression_level: None,
                // These regexes will overlap on all NM_ transcripts
                include_transcripts: Some("NM_.*".into()),
                exclude_transcripts: Some("NM_.*".into()),
            },
        );

        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("matches both include and exclude regular expressions"));

        Ok(())
    }

    /// A database with transcript `NM_007294.4` of BRCA1.
    #[allow(deprecated)]
    fn brca1_tx_db(gene_id: &str, legacy_gene_id: &str, hgnc_id: Option<&str>) -> TranscriptDb {
        TranscriptDb {
            transcripts: vec![Transcript {
                id: "NM_007294.4".into(),
                gene_symbol: "BRCA1".into(),
                legacy_gene_id: legacy_gene_id.into(),
                gene_id: gene_id.into(),
                hgnc_id: hgnc_id.map(String::from),
                ..Default::default()
            }],
            gene_to_tx: vec![],
        }
    }

    /// `--hgnc-id` selects by HGNC ID, with or without `HGNC:`. A database built before
    /// `gene_id` and `hgnc_id` existed names the HGNC ID in `legacy_gene_id`.
    #[rstest::rstest]
    #[case::hgnc_id("672", "GENE:672", Some("HGNC:1100"), "HGNC:1100")]
    #[case::hgnc_id_number("672", "GENE:672", Some("HGNC:1100"), "1100")]
    #[case::old_db("", "HGNC:1100", None, "HGNC:1100")]
    #[case::old_db_number("", "HGNC:1100", None, "1100")]
    fn hgnc_id_selects_transcripts(
        #[case] gene_id: &str,
        #[case] legacy_gene_id: &str,
        #[case] hgnc_id: Option<&str>,
        #[case] query: &str,
    ) -> Result<(), anyhow::Error> {
        let tx_db = brca1_tx_db(gene_id, legacy_gene_id, hgnc_id);
        let (_, tx_ids, _) = super::_extract_transcripts_by_hgnc_id(&tx_db, &[query.into()])?;
        assert_eq!(tx_ids.into_iter().collect::<Vec<_>>(), ["NM_007294.4"]);
        Ok(())
    }

    /// A gene keeps its entry if it holds a selected transcript, also if its
    /// `legacy_gene_id` differs from that of the transcript.
    #[test]
    #[allow(deprecated)]
    fn gene_with_selected_transcript_is_kept() -> Result<(), anyhow::Error> {
        let container = TxSeqDatabase {
            tx_db: Some(TranscriptDb {
                transcripts: vec![Transcript {
                    id: "NM_000014.6".into(),
                    gene_symbol: "A2M".into(),
                    legacy_gene_id: "GENE:2".into(),
                    gene_id: "2".into(),
                    hgnc_id: Some("HGNC:7".into()),
                    ..Default::default()
                }],
                gene_to_tx: vec![GeneToTxId {
                    legacy_gene_id: "HGNC:7".into(),
                    tx_ids: vec!["NM_000014.6".into()],
                    gene_id: "2".into(),
                    hgnc_id: Some("HGNC:7".into()),
                    ..Default::default()
                }],
            }),
            seq_db: Some(SequenceDb::default()),
            ..Default::default()
        };
        let args = super::Args {
            path_in: Default::default(),
            path_out: Default::default(),
            selection: Selection {
                transcript_id: Some(vec!["NM_000014.6".into()]),
                ..Default::default()
            },
            compression_level: None,
            include_transcripts: None,
            exclude_transcripts: None,
        };

        let subset = super::subset_tx_db(&container, &args)?;
        let gene_to_tx = subset
            .tx_db
            .map(|tx_db| tx_db.gene_to_tx)
            .unwrap_or_default();
        assert_eq!(
            gene_to_tx,
            container
                .tx_db
                .map(|tx_db| tx_db.gene_to_tx)
                .unwrap_or_default()
        );

        Ok(())
    }
}
