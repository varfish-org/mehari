//! Transcript database.

use std::collections::HashSet;
use std::{collections::HashMap, fs::File, io::Write, path::PathBuf, time::Instant};

use clap::Parser;
use flatbuffers::FlatBufferBuilder;
use hgvs::data::cdot::json::models;
use seqrepo::SeqRepo;

use crate::common::trace_rss_now;
use crate::world_flatbuffers::mehari::{
    ExonAlignment, ExonAlignmentArgs, GeneToTxId, GeneToTxIdArgs, GenomeAlignment,
    GenomeAlignmentArgs, GenomeBuild, SequenceDb, SequenceDbArgs, Strand, Transcript,
    TranscriptArgs, TranscriptBiotype, TranscriptDb, TranscriptDbArgs, TranscriptTag,
    TxSeqDatabase, TxSeqDatabaseArgs,
};

/// Command line arguments for `db create txs` sub command.
#[derive(Parser, Debug)]
#[command(about = "Construct mehari transcripts and sequence database", long_about = None)]
pub struct Args {
    /// Path to output flatbuffers file to write to.
    #[arg(long)]
    pub path_out: String,
    /// Paths to the cdot JSON transcripts to import.
    #[arg(long, required = true)]
    pub path_cdot_json: Vec<String>,
    /// Path to the seqrepo instance directory to use.
    #[arg(long)]
    pub path_seqrepo_instance: String,
}

/// Load and extract from cdot JSON.
fn load_and_extract(
    json_path: &str,
    transcript_ids_for_gene: &mut HashMap<String, Vec<String>>,
    genes: &mut HashMap<String, models::Gene>,
    transcripts: &mut HashMap<String, models::Transcript>,
) -> Result<(), anyhow::Error> {
    log::debug!("Loading cdot transcripts from {:?}", json_path);
    let start = Instant::now();
    let models::Container {
        genes: c_genes,
        transcripts: c_txs,
        ..
    } = if json_path.ends_with(".gz") {
        serde_json::from_reader(flate2::bufread::GzDecoder::new(std::io::BufReader::new(
            std::fs::File::open(json_path)?,
        )))?
    } else {
        serde_json::from_reader(std::io::BufReader::new(std::fs::File::open(json_path)?))?
    };
    log::debug!(
        "loading / deserializing {} genes and {} transcripts from cdot took {:?}",
        c_genes.len(),
        c_txs.len(),
        start.elapsed()
    );

    let start = Instant::now();
    c_genes
        .values()
        .filter(|gene| {
            gene.gene_symbol.is_some()
                && !gene.gene_symbol.as_ref().unwrap().is_empty()
                && gene.map_location.is_some()
                && !gene.map_location.as_ref().unwrap().is_empty()
        })
        .for_each(|gene| {
            let gene_symbol = gene.gene_symbol.as_ref().unwrap().clone();
            transcript_ids_for_gene
                .entry(gene_symbol.clone())
                .or_insert(Vec::new());
            genes.insert(gene_symbol, gene.clone());
        });
    c_txs
        .values()
        .filter(|tx| tx.gene_name.is_some() && !tx.gene_name.as_ref().unwrap().is_empty())
        .for_each(|tx| {
            let gene_name = tx.gene_name.as_ref().unwrap();
            transcript_ids_for_gene
                .get_mut(gene_name)
                .unwrap_or_else(|| panic!("tx {:?} for unknown gene {:?}", tx.id, gene_name))
                .push(tx.id.clone());
            transcripts.insert(tx.id.clone(), tx.clone());
        });
    log::debug!("extracting datastructures took {:?}", start.elapsed());
    Ok(())
}

/// Perform flatbuffers file construction.
fn build_flatbuffers(
    path_out: &str,
    _seqrepo: SeqRepo,
    _genes: HashMap<String, models::Gene>,
    _transcripts: HashMap<String, models::Transcript>,
    transcript_ids_for_gene: HashMap<String, Vec<String>>,
) -> Result<(), anyhow::Error> {
    tracing::info!("Constructing flatbuffers file ...");
    trace_rss_now();

    // First Pass: Select transcripts for each gene.

    // Second Pass: Build flatbuffers.

    let mut builder = FlatBufferBuilder::new();

    let cigar = Some(builder.create_shared_string("x"));

    let exon = ExonAlignment::create(
        &mut builder,
        &ExonAlignmentArgs {
            alt_start_i: -1,
            alt_end_i: -1,
            ord: -1,
            alt_cds_start_i: -1,
            alt_cds_end_i: -1,
            cigar,
        },
    );

    let genome_build = GenomeBuild::Grch37;
    let contig = Some(builder.create_shared_string("x"));
    let cds_start = -1;
    let cds_end = -1;
    let strand = Strand::Plus;
    let exons = Some(builder.create_vector(&[exon]));

    let aln = GenomeAlignment::create(
        &mut builder,
        &GenomeAlignmentArgs {
            genome_build,
            contig,
            cds_start,
            cds_end,
            strand,
            exons,
        },
    );

    let id = Some(builder.create_shared_string("x"));
    let gene_name = Some(builder.create_shared_string("x"));
    let gene_id = Some(builder.create_shared_string("x"));
    let biotype = TranscriptBiotype::Coding;
    let tags = TranscriptTag::Basic.0 as u8;
    let protein = Some(builder.create_shared_string("x"));
    let start_codon = -1;
    let stop_codon = -1;
    let genome_alignments = Some(builder.create_vector(&[aln]));

    let tx = Transcript::create(
        &mut builder,
        &TranscriptArgs {
            id,
            gene_name,
            gene_id,
            biotype,
            tags,
            protein,
            start_codon,
            stop_codon,
            genome_alignments,
        },
    );

    let transcripts = builder.create_vector(&[tx]);

    // Build mapping of gene HGNC symbol to transcript IDs.
    let gene_to_tx = {
        let items = transcript_ids_for_gene
            .iter()
            .map(|(gene_name, tx_ids)| {
                let gene_name = Some(builder.create_shared_string(gene_name));
                let tx_ids = tx_ids
                    .iter()
                    .map(|s| builder.create_shared_string(s))
                    .collect::<Vec<_>>();
                let tx_ids = Some(builder.create_vector(&tx_ids));
                GeneToTxId::create(&mut builder, &GeneToTxIdArgs { gene_name, tx_ids })
            })
            .collect::<Vec<_>>();
        builder.create_vector(&items)
    };

    let tx_db = TranscriptDb::create(
        &mut builder,
        &TranscriptDbArgs {
            transcripts: Some(transcripts),
            gene_to_tx: Some(gene_to_tx),
        },
    );

    let s = builder.create_shared_string("x");
    let aliases = builder.create_vector(&[s]);
    let aliases_idx = builder.create_vector(&[1u32]);
    let seqs = builder.create_vector(&[s]);

    let seq_db = SequenceDb::create(
        &mut builder,
        &SequenceDbArgs {
            aliases: Some(aliases),
            aliases_idx: Some(aliases_idx),
            seqs: Some(seqs),
        },
    );

    let tx_seq_db = TxSeqDatabase::create(
        &mut builder,
        &TxSeqDatabaseArgs {
            tx_db: Some(tx_db),
            seq_db: Some(seq_db),
        },
    );

    builder.finish_minimal(tx_seq_db);
    let mut output_file = File::create(path_out)?;
    output_file.write_all(builder.finished_data())?;
    output_file.flush()?;

    trace_rss_now();
    tracing::info!("... done with constructing flatbuffers file");
    Ok(())
}

/// Filter transcripts for gene.
///
/// We employ the following rules:
///
/// - Remove redundant transcripts with the same identifier and pick only the
///   transcripts that have the highest version number for one assembly.
/// - Do not pick any `XM_`/`XR_` (NCBI predicted only) transcripts.
/// - Do not pick any `NR_` transcripts when there are coding `NM_` transcripts.
fn filter_transcripts(
    transcripts: &HashMap<String, models::Transcript>,
    transcript_ids_for_gene: HashMap<String, Vec<String>>,
) -> HashMap<String, Vec<String>> {
    transcript_ids_for_gene
        .into_iter()
        .map(|(gene_symbol, prev_tx_ids)| {
            // Split off transcript versions from accessions and look for NM transcript.
            let mut seen_nm = false;
            let mut versioned: Vec<_> = prev_tx_ids
                .iter()
                .map(|tx_id| {
                    if tx_id.starts_with("NM_") {
                        seen_nm = true;
                    }
                    let s: Vec<_> = tx_id.split('.').collect();
                    (s[0], s[1].parse::<u32>().expect("invalid version"))
                })
                .collect();
            // Sort descendingly by version.
            versioned.sort_by(|a, b| b.1.cmp(&a.1));

            // Build `next_tx_ids`.
            let mut seen_ac = HashSet::new();
            let mut next_tx_ids = Vec::new();
            for (ac, version) in versioned {
                let full_ac = format!("{}.{}", &ac, version);
                let ac = ac.to_string();

                let releases = transcripts
                    .get(&full_ac)
                    .map(|tx| tx.genome_builds.keys().cloned().collect::<Vec<_>>())
                    .unwrap_or_default();

                for release in releases {
                    #[allow(clippy::if_same_then_else)]
                    if seen_ac.contains(&(ac.clone(), release.clone())) {
                        continue; // skip, already have later version
                    } else if ac.starts_with("NR_") && seen_nm {
                        continue; // skip NR transcript as we have NM one
                    } else if ac.starts_with('X') {
                        continue; // skip XR/XM transcript
                    } else {
                        next_tx_ids.push(full_ac.clone());
                        seen_ac.insert((ac.clone(), release));
                    }
                }
            }

            next_tx_ids.sort();
            next_tx_ids.dedup();

            (gene_symbol, next_tx_ids)
        })
        .collect()
}

/// Main entry point for `db create txs` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!(
        "Building transcript and sequence database file\ncommon args: {:#?}\nargs: {:#?}",
        common,
        args
    );

    tracing::info!("Opening seqrepo...");
    let start = Instant::now();
    let seqrepo = PathBuf::from(&args.path_seqrepo_instance);
    let path = seqrepo
        .parent()
        .ok_or(anyhow::anyhow!(
            "Could not get parent from {}",
            &args.path_seqrepo_instance
        ))?
        .to_str()
        .unwrap()
        .to_string();
    let instance = seqrepo
        .file_name()
        .ok_or(anyhow::anyhow!(
            "Could not get basename from {}",
            &args.path_seqrepo_instance
        ))?
        .to_str()
        .unwrap()
        .to_string();
    let seqrepo = SeqRepo::new(path, &instance)?;
    tracing::info!("... seqrepo opened in {:?}", start.elapsed());

    tracing::info!("Loading cdot JSON files ...");
    let start = Instant::now();
    let mut genes = HashMap::new();
    let mut transcripts = HashMap::new();
    let mut transcript_ids_for_gene = HashMap::new();
    for json_path in &args.path_cdot_json {
        load_and_extract(
            json_path,
            &mut transcript_ids_for_gene,
            &mut genes,
            &mut transcripts,
        )?;
    }
    tracing::info!(
        "... done loading cdot JSON files in {:?} -- #genes = {}, #transcripts = {}, #transcript_ids_for_gene = {}",
        start.elapsed(),
        genes.len(),
        transcripts.len(),
        transcript_ids_for_gene.len()
    );

    // Filter out redundant (NR if we have NM, older versions) and predicted only transcripts.
    tracing::info!("Filtering transcripts ...");
    let start = Instant::now();
    transcript_ids_for_gene = filter_transcripts(&transcripts, transcript_ids_for_gene);
    tracing::info!("... done filtering transcripts in {:?}", start.elapsed());

    build_flatbuffers(
        &args.path_out,
        seqrepo,
        genes,
        transcripts,
        transcript_ids_for_gene,
    )?;

    tracing::info!("Done building transcript and sequence database file");
    Ok(())
}

#[cfg(test)]
pub mod test {
    use std::collections::HashMap;

    use pretty_assertions::assert_eq;

    use super::{filter_transcripts, load_and_extract};

    #[test]
    fn filter_transcripts_brca1() -> Result<(), anyhow::Error> {
        let mut genes = HashMap::new();
        let mut transcripts = HashMap::new();
        let mut transcript_ids_for_gene = HashMap::new();
        load_and_extract(
            "tests/data/db/create/txs/cdot-0.2.12.refseq.grch37_grch38.brca1.json",
            &mut transcript_ids_for_gene,
            &mut genes,
            &mut transcripts,
        )?;

        assert_eq!(
            &transcript_ids_for_gene
                .get("BRCA1")
                .unwrap()
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>(),
            &vec![
                "NM_007294.3",
                "NM_007294.4",
                "NM_007297.3",
                "NM_007297.4",
                "NM_007298.3",
                "NM_007299.3",
                "NM_007299.4",
                "NM_007300.3",
                "NM_007300.4",
                "NR_027676.1",
                "NR_027676.2",
                "XM_006722029.1",
                "XM_006722030.1",
                "XM_006722031.1",
                "XM_006722032.1",
                "XM_006722033.1",
                "XM_006722034.1",
                "XM_006722035.1",
                "XM_006722036.1",
                "XM_006722037.1",
                "XM_006722038.1",
                "XM_006722039.1",
                "XM_006722040.1",
                "XM_006722041.1"
            ]
        );
        let filtered = filter_transcripts(&transcripts, transcript_ids_for_gene);
        assert_eq!(
            &filtered
                .get("BRCA1")
                .unwrap()
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>(),
            &vec![
                "NM_007294.4",
                "NM_007297.4",
                "NM_007298.3",
                "NM_007299.4",
                "NM_007300.4"
            ]
        );
        Ok(())
    }
}
