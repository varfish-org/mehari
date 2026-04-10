//! Transcript database.

use crate::annotate::seqvars::ann::FeatureTag;
use crate::common::trace_rss_now;
use crate::pbs::txs::{SourceVersion, TxSeqDatabase};
use anyhow::{Error, anyhow};
use cli::Args;
use hgvs::data::cdot::json::models as cdot_models;
use itertools::Itertools;
use models::{LabelEntry, ReportEntry};
use prost::Message;
use serde_json::json;
use serde_with::DisplayFromStr;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::{io::Write, time::Instant};
use thousands::Separable;

mod build;
pub(crate) mod cdot;
pub mod cli;
mod filter;
mod gff3;
pub mod models;
pub mod reference;

use crate::db::create::build::build_protobuf;
use crate::db::create::cdot::load_cdot;
use crate::db::create::filter::{
    filter_empty_gene_id_mappings, filter_genes, filter_initial_gene_id_entries,
    filter_transcripts, filter_transcripts_with_sequence,
};
use crate::db::create::gff3::load_gff3;
use models::*;

fn txid_to_label(
    label_tsv_path: impl AsRef<Path>,
) -> Result<HashMap<TranscriptId, Vec<FeatureTag>>, Error> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_path(label_tsv_path.as_ref())?;

    rdr.deserialize()
        .map(|result| {
            result
                .map_err(anyhow::Error::from)
                .and_then(|entry: LabelEntry| {
                    TranscriptId::try_new(entry.transcript_id)
                        .map(|txid| {
                            (
                                txid,
                                entry
                                    .label
                                    .split(',')
                                    .map(|s| {
                                        let cdot_tag = cdot_models::str_to_tag(s);
                                        FeatureTag::from(cdot_tag)
                                    })
                                    .collect::<Vec<_>>(),
                            )
                        })
                        .map_err(anyhow::Error::from)
                })
        })
        .collect()
}

/// Load the annotations (JSON or GFF3).
fn load_annotations(args: &Args) -> Result<TranscriptLoader, Error> {
    tracing::info!("Loading annotations …");
    let start = Instant::now();
    let labels = args
        .mane_transcripts
        .as_ref()
        .map(txid_to_label)
        .transpose()?;
    let loaders = args
        .annotation
        .iter()
        .map(|path| {
            let mut loader = TranscriptLoader::new(args.assembly.clone(), args.disable_filters);

            let ext = path.extension().unwrap_or_default().to_string_lossy();
            let file_stem = path.file_stem().unwrap_or_default().to_string_lossy();
            let is_gff3 = ext.ends_with("gff")
                || ext.ends_with("gff3")
                || (ext == "gz" && (file_stem.ends_with("gff3") || file_stem.ends_with("gff")));

            if is_gff3 {
                load_gff3(&mut loader, path).map(|_| loader)
            } else {
                load_cdot(&mut loader, path).map(|_| loader)
            }
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let mut merged = loaders
        .into_iter()
        .reduce(|mut a, mut b| {
            a.merge(&mut b);
            a
        })
        .unwrap();
    merged.apply_fixes(&labels);

    tracing::info!(
        "… done loading annotations in {:?} -- #transcripts = {}, #gene_ids = {}",
        start.elapsed(),
        merged
            .transcript_id_to_transcript
            .len()
            .separate_with_underscores(),
        merged
            .gene_id_to_transcript_ids
            .len()
            .separate_with_underscores()
    );

    Ok(merged)
}

/// Main entry point for `db create txs` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), Error> {
    fn _run(common: &crate::common::Args, args: &Args) -> Result<(), Error> {
        let mut report_file =
            File::create(format!("{}.report.jsonl", args.output.display())).map(BufWriter::new)?;
        let mut report = |r: ReportEntry| -> Result<(), Error> {
            writeln!(report_file, "{}", serde_json::to_string(&r)?)?;
            Ok(())
        };
        tracing::info!(
            "Building transcript and sequence database file\ncommon args: {:#?}\nargs: {:#?}",
            common,
            args
        );

        let mut tx_data = load_annotations(args)?;
        for (id, fix) in tx_data.fixes.iter() {
            report(ReportEntry::Fix(LogFix {
                source: "annotations".into(),
                fix: *fix,
                id: id.clone(),
                gene_name: tx_data.gene_name(id),
                tags: tx_data.tags(id),
            }))?;
        }

        let raw_tx_data = tx_data.clone();
        trace_rss_now();
        report(ReportEntry::Log(json!({
            "source": "annotations",
            "total_transcripts": tx_data.transcript_id_to_transcript.len(),
            "total_gene_ids": tx_data.gene_id_to_transcript_ids.len()
        })))?;

        let tx_data_ = &mut tx_data;

        // … then filter gene id entries with no transcripts to boot …
        filter_initial_gene_id_entries(tx_data_)?;
        // … then filter genes (missing gene id and/or symbol) …
        filter_genes(tx_data_)?;
        // … then filter transcripts …
        filter_transcripts(tx_data_)?;
        // … ensure there are no gene keys without associated transcripts left …
        filter_empty_gene_id_mappings(tx_data_)?;

        // Open seqrepo / FASTA …
        let mut seq_provider = reference::open_sequence_provider(args)?;
        // … and filter transcripts based on their sequences,
        // e.g. checking whether their translation contains a stop codon …
        let mut sequence_map = filter_transcripts_with_sequence(tx_data_, &mut seq_provider)?;
        filter_empty_gene_id_mappings(tx_data_)?;
        // … if there are genes with no transcripts left, check whether they are pseudogenes …
        tx_data_.update_pseudogene_status()?;

        // … trigger the discard routine, but do not remove anything, just make sure they are consistent,
        // and report the stats.
        let remove = false;
        tx_data_.discard(remove)?;

        // … and update all discard annotations …
        tx_data_.propagate_discard_reasons(&raw_tx_data)?;

        trace_rss_now();

        report(ReportEntry::Log(json!({
            "source": "annotations_filtered",
            "total_transcripts": tx_data_.transcript_id_to_transcript.len(),
            "total_gene_ids": tx_data_.gene_id_to_transcript_ids.len()
        })))?;

        // For some stats, count number of chrMt, MANE Select and MANE Plus Clinical transcripts.
        let (n_mt, n_mane_select, n_mane_plus_clinical) = tx_data_.gather_transcript_stats()?;

        report(ReportEntry::Log(json!({
            "source": "annotations_filtered",
            "n_mt": n_mt,
            "n_mane_select": n_mane_select,
            "n_mane_plus_clinical": n_mane_plus_clinical
        })))?;

        let source_version = bundle_source_version_information(args);

        // … and finally construct protobuf txdb data structures.
        let tx_db = build_protobuf(tx_data_, &mut sequence_map, source_version)?;

        // List all discarded transcripts and genes.
        for (id, reason) in tx_data.discards.into_iter().sorted_unstable() {
            if reason.intersects(Reason::hard()) {
                // if it has at least one hard reason → actually discarded.
                report(ReportEntry::Discard(Discard {
                    source: "protobuf".into(),
                    reason,
                    id: id.clone(),
                    gene_name: raw_tx_data.gene_name(&id),
                    tags: raw_tx_data.tags(&id),
                }))?;
            } else {
                // only soft reasons → kept but flagged
                report(ReportEntry::SoftFilter(SoftFilter {
                    source: "protobuf".into(),
                    reason,
                    id: id.clone(),
                    gene_name: raw_tx_data.gene_name(&id),
                    tags: raw_tx_data.tags(&id),
                }))?;
            }
        }
        trace_rss_now();

        write_tx_db(tx_db, &args.output, args.compression_level)?;

        tracing::info!("Done building transcript and sequence database file");
        Ok(())
    }

    fn bundle_source_version_information(args: &Args) -> SourceVersion {
        let assembly = args.assembly.clone();

        let assembly_version = args.assembly_version.clone();

        let source_name = args.transcript_source.clone();

        let source_version = args.transcript_source_version.clone().unwrap_or("".into());
        let annotation_version = args.annotation_version.clone().unwrap_or("".into());
        let annotation_name = args
            .annotation
            .iter()
            .map(|p| p.file_stem().unwrap().to_string_lossy())
            .collect::<Vec<_>>()
            .join(",");

        SourceVersion {
            mehari_version: crate::common::version().to_string(),
            assembly: assembly.to_string(),
            assembly_version,
            source_name: source_name.to_string(),
            source_version,
            annotation_version,
            annotation_name,
        }
    }

    let threadpool = rayon::ThreadPoolBuilder::default()
        .num_threads(args.threads)
        .build()?;

    threadpool.install(|| _run(common, args))
}

pub(crate) fn write_tx_db(
    tx_db: TxSeqDatabase,
    path: impl AsRef<Path>,
    compression_level: i32,
) -> Result<(), Error> {
    tracing::info!("Writing out final database …");
    let path = path.as_ref();
    let mut buf = prost::bytes::BytesMut::with_capacity(tx_db.encoded_len());
    tx_db
        .encode(&mut buf)
        .map_err(|e| anyhow!("failed to encode: {}", e))?;
    tracing::info!("  … done constructing final tx and seq database");

    // Write out the final transcript and sequence database.
    tracing::info!("  Writing out final database …");
    // Open file and if necessary, wrap in a decompressor.
    let file = std::fs::File::create(path)
        .map_err(|e| anyhow!("failed to create file {}: {}", path.display(), e))?;
    let ext = path.extension().map(|s| s.to_str());
    let mut writer: Box<dyn Write> = if ext == Some(Some("zst")) {
        Box::new(
            zstd::Encoder::new(file, compression_level)
                .map_err(|e| anyhow!("failed to open zstd encoder for {}: {}", path.display(), e))?
                .auto_finish(),
        )
    } else {
        Box::new(file)
    };
    writer
        .write_all(&buf)
        .map_err(|e| anyhow!("failed to write to {}: {}", path.display(), e))?;
    tracing::info!("  … done writing out final database");

    Ok(())
}

#[cfg(test)]
pub mod test {
    use std::path::{Path, PathBuf};

    use clap_verbosity_flag::Verbosity;
    use itertools::Itertools;
    use rstest::rstest;
    use temp_testdir::TempDir;

    use crate::common::Args as CommonArgs;
    use crate::db::create::cdot::load_cdot;
    use crate::db::create::cli::Args;
    use crate::db::create::filter::filter_transcripts;
    use crate::db::create::models::GeneId;
    use crate::db::create::models::TranscriptLoader;
    use crate::db::dump;

    use super::run;

    #[test]
    fn filter_transcripts_brca1() -> Result<(), anyhow::Error> {
        let path_tsv = Path::new("tests/data/db/create/txs/txs_main.tsv");
        let mut tx_data = TranscriptLoader::new("GRCh37".to_string(), false);
        let labels = super::txid_to_label(path_tsv)?;
        load_cdot(
            &mut tx_data,
            Path::new("tests/data/db/create/txs/cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json"),
        )?;
        tx_data.apply_fixes(&Some(labels));

        eprintln!("{:#?}", &tx_data.gene_id_to_transcript_ids);
        insta::assert_yaml_snapshot!(
            tx_data
                .gene_id_to_transcript_ids
                .get(&GeneId::Hgnc(1100))
                .unwrap()
                .iter()
                .map(|s| s.as_str())
                .sorted_unstable()
                .collect::<Vec<_>>()
        );

        filter_transcripts(&mut tx_data)?;
        insta::assert_yaml_snapshot!(
            tx_data
                .gene_id_to_transcript_ids
                .get(&GeneId::Hgnc(1100))
                .unwrap()
                .iter()
                .map(|s| s.as_str())
                .sorted_unstable()
                .collect::<Vec<_>>()
        );

        insta::assert_snapshot!(&tx_data.annotation_version);

        Ok(())
    }

    #[rstest]
    #[case("grch37")]
    #[case("grch38")]
    fn run_smoke_brca1_opa1(#[case] assembly: String) -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            output: tmp_dir.join("out.bin.zst"),
            annotation: vec![PathBuf::from(
                "tests/data/db/create/txs/cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json",
            )],
            mane_transcripts: Some(PathBuf::from("tests/data/db/create/txs/txs_main.tsv")),
            seqrepo: Some(PathBuf::from("tests/data/db/create/txs/latest")),
            transcript_sequences: None,
            assembly: assembly.clone(),
            assembly_version: None,
            transcript_source: "refseq".to_string(),
            transcript_source_version: None,
            disable_filters: false,
            threads: 1,
            annotation_version: Some("0.2.22".to_string()),
            compression_level: 19,
        };

        run(&common_args, &args)?;

        let mut buf: Vec<u8> = Vec::new();
        dump::run_with_write(
            &Default::default(),
            &dump::Args {
                path_db: tmp_dir.join("out.bin.zst"),
            },
            &mut buf,
        )?;
        crate::common::set_snapshot_suffix!("{}", assembly.to_lowercase());
        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }

    #[test]
    fn run_smoke_selenoproteins() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            output: tmp_dir.join("out.bin.zst"),
            annotation: vec![PathBuf::from(
                "tests/data/db/create/seleonoproteins/cdot-0.2.22.refseq.grch38.selenon.json",
            )],
            mane_transcripts: Some(PathBuf::from("tests/data/db/create/txs/txs_main.tsv")),
            seqrepo: Some(PathBuf::from("tests/data/db/create/seleonoproteins/latest")),
            transcript_sequences: None,
            assembly: "grch38".to_string(),
            assembly_version: None,
            transcript_source: "refseq".to_string(),
            transcript_source_version: None,
            disable_filters: false,
            threads: 1,
            annotation_version: Some("0.2.22".to_string()),
            compression_level: 19,
        };

        run(&common_args, &args)?;

        let mut buf: Vec<u8> = Vec::new();
        dump::run_with_write(
            &Default::default(),
            &dump::Args {
                path_db: tmp_dir.join("out.bin.zst"),
            },
            &mut buf,
        )?;
        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[test]
    fn run_smoke_mitochondrial() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(5, 0),
        };
        let args = Args {
            output: tmp_dir.join("out.bin.zst"),
            annotation: vec![PathBuf::from(
                "tests/data/db/create/mitochondrial/cdot-0.2.23.ensembl.chrMT.grch37.gff3.json",
            )],
            mane_transcripts: None,
            seqrepo: Some(PathBuf::from("tests/data/db/create/mitochondrial/latest")),
            transcript_sequences: None,
            assembly: "grch37".to_string(),
            assembly_version: None,
            transcript_source: "ensembl".to_string(),
            transcript_source_version: Some("98".into()),
            disable_filters: false,
            threads: 1,
            annotation_version: Some("0.2.23".to_string()),
            compression_level: 19,
        };

        run(&common_args, &args)?;

        let mut buf: Vec<u8> = Vec::new();
        dump::run_with_write(
            &Default::default(),
            &dump::Args {
                path_db: tmp_dir.join("out.bin.zst"),
            },
            &mut buf,
        )?;
        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }
}
