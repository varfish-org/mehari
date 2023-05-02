//! Transcript database.

use std::collections::HashSet;
use std::fs::File;
use std::path::Path;
use std::{collections::HashMap, io::Write, path::PathBuf, time::Instant};

use anyhow::anyhow;
use clap::Parser;
use hgvs::data::cdot::json::models;
use hgvs::sequences::{translate_cds, TranslationTable};
use indicatif::{ProgressBar, ProgressStyle};
use prost::Message;
use seqrepo::{AliasOrSeqId, Interface, SeqRepo};
use thousands::Separable;

use crate::common::{trace_rss_now, GenomeRelease};

lazy_static::lazy_static! {
    /// Progress bar style to use.
    pub static ref PROGRESS_STYLE: ProgressStyle = ProgressStyle::with_template(
        "[{elapsed_precise}] [{wide_bar:.cyan/blue}] {human_pos}/{human_len} ({eta})",
    )
    .unwrap();
}

/// Data structures for (de-)serialization as generated by `prost-build`.
pub mod data {
    include!(concat!(env!("OUT_DIR"), "/mehari.db.create.txs.data.rs"));
}

/// Command line arguments for `db create txs` sub command.
#[derive(Parser, Debug)]
#[command(about = "Construct mehari transcripts and sequence database", long_about = None)]
pub struct Args {
    /// Genome release to extract transcripts for.
    #[arg(long)]
    pub genome_release: GenomeRelease,
    /// Path to output protobuf file to write to.
    #[arg(long)]
    pub path_out: PathBuf,
    /// Paths to the cdot JSON transcripts to import.
    #[arg(long, required = true)]
    pub path_cdot_json: Vec<PathBuf>,
    /// Path to the seqrepo instance directory to use.
    #[arg(long)]
    pub path_seqrepo_instance: PathBuf,
    /// Maximal number of transcripts to process.
    #[arg(long)]
    pub max_txs: Option<u32>,
    /// Limit transcript database to the following HGNC symbols.  Useful for
    /// building test databases.
    #[arg(long)]
    pub gene_symbols: Option<Vec<String>>,
}

/// Load and extract from cdot JSON.
fn load_and_extract(
    json_path: &Path,
    transcript_ids_for_gene: &mut HashMap<String, Vec<String>>,
    genes: &mut HashMap<String, models::Gene>,
    transcripts: &mut HashMap<String, models::Transcript>,
    genome_release: GenomeRelease,
    report_file: &mut File,
) -> Result<(), anyhow::Error> {
    writeln!(report_file, "genome_release\t{:?}", genome_release)?;

    tracing::info!("Loading cdot transcripts from {:?}", json_path);
    writeln!(report_file, "cdot_json_path\t{:?}", json_path)?;
    let start = Instant::now();
    let models::Container {
        genes: c_genes,
        transcripts: c_txs,
        ..
    } = if json_path.extension().unwrap_or_default() == "gz" {
        tracing::info!("(from gzip compressed file)");
        serde_json::from_reader(std::io::BufReader::new(flate2::read::GzDecoder::new(
            File::open(json_path)?,
        )))?
    } else {
        tracing::info!("(from uncompressed file)");
        serde_json::from_reader(std::io::BufReader::new(File::open(json_path)?))?
    };
    tracing::info!(
        "loading / deserializing {} genes and {} transcripts from cdot took {:?}",
        c_genes.len().separate_with_commas(),
        c_txs.len().separate_with_commas(),
        start.elapsed()
    );

    let start = Instant::now();
    writeln!(report_file, "total_genes\t{}", c_genes.len())?;
    c_genes
        .values()
        .filter(|gene| {
            gene.gene_symbol.is_some()
                && !gene.gene_symbol.as_ref().unwrap().is_empty()
                && gene.map_location.is_some()
                && !gene.map_location.as_ref().unwrap().is_empty()
                && gene.hgnc.is_some()
                && !gene.hgnc.as_ref().unwrap().is_empty()
        })
        .for_each(|gene| {
            let gene_symbol = gene.gene_symbol.as_ref().unwrap().clone();
            transcript_ids_for_gene
                .entry(gene_symbol.clone())
                .or_insert(Vec::new());
            genes.insert(gene_symbol, gene.clone());
        });
    writeln!(
        report_file,
        "genes with gene_symbol, map_location, hgnc\t{}",
        genes.len()
    )?;
    tracing::info!(
        "Processed {} genes; total gene count: {}",
        c_genes.len().separate_with_commas(),
        genes.len()
    );
    tracing::debug!(
        "some 10 genes: {:?}",
        genes.keys().take(10).collect::<Vec<_>>()
    );

    tracing::info!("Processing transcripts");
    writeln!(report_file, "total_transcripts\t{}", c_txs.len())?;
    c_txs
        .values()
        .map(|tx| models::Transcript {
            genome_builds: tx
                .genome_builds
                .iter()
                .filter(|(key, _)| {
                    matches!(
                        (key.as_str(), genome_release),
                        ("GRCh37", GenomeRelease::Grch37) | ("GRCh38", GenomeRelease::Grch38)
                    )
                })
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect(),
            ..tx.clone()
        })
        .filter(|tx| {
            tx.gene_name.is_some()
                && !tx.gene_name.as_ref().unwrap().is_empty()
                && genes.contains_key(tx.gene_name.as_ref().unwrap())
                && !tx.genome_builds.is_empty()
        })
        .for_each(|tx| {
            let gene_name = tx.gene_name.as_ref().unwrap();
            transcript_ids_for_gene
                .get_mut(gene_name)
                .unwrap_or_else(|| panic!("tx {:?} for unknown gene {:?}", tx.id, gene_name))
                .push(tx.id.clone());
            transcripts.insert(tx.id.clone(), tx.clone());
        });
    writeln!(
        report_file,
        "transcripts with alignment on genome and link to selected gene\t{}",
        transcripts.len()
    )?;
    tracing::info!(
        "Processed {} genes; total transcript count: {}",
        c_txs.len().separate_with_commas(),
        transcripts.len().separate_with_commas()
    );
    tracing::info!("extracting datastructures took {:?}", start.elapsed());
    Ok(())
}

/// Perform protobuf file construction.
///
/// This can be done by simply converting the models from HGVS to the prost generated data structures.
fn build_protobuf(
    path_out: &Path,
    seqrepo: SeqRepo,
    tx_data: TranscriptData,
    is_silent: bool,
    report_file: &mut File,
) -> Result<(), anyhow::Error> {
    let TranscriptData {
        genes,
        transcripts,
        transcript_ids_for_gene,
    } = tx_data;

    tracing::info!("Constructing protobuf data structures ...");
    trace_rss_now();

    // Construct sequence database.
    tracing::info!("  Constructing sequence database ...");
    let mut tx_skipped_noseq = HashSet::new(); // skipped because of missing sequence
    let mut tx_skipped_nostop = HashSet::new(); // skipped because of missing stop codon
    let seq_db = {
        // Insert into flatbuffer and keep track of pointers in `Vec`s.
        let mut aliases = Vec::new();
        let mut aliases_idx = Vec::new();
        let mut seqs = Vec::new();
        let pb = if is_silent {
            ProgressBar::hidden()
        } else {
            ProgressBar::new(transcripts.len() as u64)
        };
        pb.set_style(PROGRESS_STYLE.clone());
        for (tx_id, tx) in &transcripts {
            pb.inc(1);
            let res_seq = seqrepo.fetch_sequence(&AliasOrSeqId::Alias {
                value: tx_id.clone(),
                namespace: None,
            });
            let seq = if let Ok(seq) = res_seq {
                seq
            } else {
                tracing::debug!("Skipping transcript {} because of missing sequence", tx_id);
                writeln!(
                    report_file,
                    "skip transcript because it has no sequence\t{}",
                    tx_id
                )?;
                tx_skipped_noseq.insert(tx_id.clone());
                continue;
            };

            // Skip transcript if it is coding and the translated CDS does not have a stop codon.
            if let Some(cds_start) = tx.start_codon {
                let cds_start = cds_start as usize;
                let cds_end = tx.stop_codon.expect("must be some if start_codon is some") as usize;
                if cds_end > seq.len() {
                    tracing::error!(
                        "CDS end {} is larger than sequence length {} for {}",
                        cds_end,
                        seq.len(),
                        tx_id
                    );
                    writeln!(
                        report_file,
                        "skip transcript CDS end {} is longer than sequence length {} for\t{}",
                        cds_end,
                        seq.len(),
                        tx_id
                    )?;
                    continue;
                }
                let tx_seq_to_translate = &seq[cds_start..cds_end];
                let aa_sequence =
                    translate_cds(tx_seq_to_translate, true, "*", TranslationTable::Standard)?;
                if !aa_sequence.ends_with('*') {
                    tracing::debug!(
                        "Skipping transcript {} because of missing stop codon in translated CDS",
                        tx_id
                    );
                    writeln!(
                        report_file,
                        "Skipping transcript {} because of missing stop codon in translated CDS",
                        tx_id
                    )?;
                    tx_skipped_nostop.insert(tx_id.clone());
                    continue;
                }
            }

            // Register sequence into flatbuffer.
            aliases.push(tx_id.clone());
            aliases_idx.push(seqs.len() as u32);
            seqs.push(seq.clone());
        }
        pb.finish_and_clear();
        // Finalize by creating `SequenceDb`.
        data::SequenceDb {
            aliases,
            aliases_idx,
            seqs,
        }
    };
    tracing::info!(
        "  ... done constructing sequence database (no seq for {} transcripts, \
        no stop codon for {}, will be skipped)",
        tx_skipped_noseq.len().separate_with_commas(),
        tx_skipped_nostop.len().separate_with_commas(),
    );

    trace_rss_now();

    tracing::info!("  Creating transcript records for each gene...");
    let data_transcripts = {
        let gene_symbols = {
            let mut gene_symbols: Vec<_> = genes.keys().cloned().collect();
            gene_symbols.sort();
            gene_symbols
        };
        let mut data_transcripts = Vec::new();
        // For each gene (in lexicographic symbol order) ...
        for gene_symbol in &gene_symbols {
            let gene = genes.get(gene_symbol).unwrap();
            let tx_ids = transcript_ids_for_gene
                .get(gene_symbol.as_str())
                .unwrap_or_else(|| panic!("No transcripts for gene {:?}", &gene_symbol));
            let tx_ids = tx_ids
                .iter()
                .filter(|tx_id| {
                    !tx_skipped_noseq.contains(*tx_id) && !tx_skipped_nostop.contains(*tx_id)
                })
                .collect::<Vec<_>>();
            if tx_ids.is_empty() {
                tracing::debug!(
                    "Skipping gene {} as all transcripts have been removed.",
                    gene_symbol
                );
                writeln!(
                    report_file,
                    "skip gene from flatbuffers because all transcripts have been removed\t{}",
                    gene_symbol
                )?;
                continue;
            }

            // ... for each transcript of the gene ...
            for tx_id in tx_ids {
                let tx_model = transcripts
                    .get(tx_id)
                    .unwrap_or_else(|| panic!("No transcript model for id {:?}", tx_id));
                // ... build genome alignment for selected:
                let mut genome_alignments = Vec::new();
                for (genome_build, alignment) in &tx_model.genome_builds {
                    // obtain basic properties
                    let genome_build = match genome_build.as_ref() {
                        "GRCh37" => data::GenomeBuild::Grch37,
                        "GRCh38" => data::GenomeBuild::Grch38,
                        _ => panic!("Unknown genome build {:?}", genome_build),
                    };
                    let models::GenomeAlignment {
                        contig,
                        cds_start,
                        cds_end,
                        ..
                    } = alignment.clone();
                    let strand = match alignment.strand {
                        models::Strand::Plus => data::Strand::Plus,
                        models::Strand::Minus => data::Strand::Minus,
                    };
                    // and construct vector of all exons
                    let exons: Vec<_> = alignment
                        .exons
                        .iter()
                        .map(|exon| {
                            let models::Exon {
                                alt_start_i,
                                alt_end_i,
                                ord,
                                alt_cds_start_i,
                                alt_cds_end_i,
                                cigar,
                            } = exon.clone();
                            data::ExonAlignment {
                                alt_start_i,
                                alt_end_i,
                                ord,
                                alt_cds_start_i: if alt_cds_start_i == -1 {
                                    None
                                } else {
                                    Some(alt_cds_start_i)
                                },
                                alt_cds_end_i: if alt_cds_end_i == -1 {
                                    None
                                } else {
                                    Some(alt_cds_end_i)
                                },
                                cigar,
                            }
                        })
                        .collect();
                    // and finally push the genome alignment
                    genome_alignments.push(data::GenomeAlignment {
                        genome_build: genome_build.into(),
                        contig,
                        cds_start,
                        cds_end,
                        strand: strand.into(),
                        exons,
                    });
                }

                // Now, just obtain the basic properties and create a new `data::Transcript`.
                let models::Gene {
                    biotype,
                    hgnc,
                    gene_symbol,
                    ..
                } = gene.clone();
                let biotype = if biotype.unwrap().contains(&models::BioType::ProteinCoding) {
                    data::TranscriptBiotype::Coding.into()
                } else {
                    data::TranscriptBiotype::NonCoding.into()
                };
                let mut tags = Vec::new();
                if let Some(tag) = tx_model.tag.as_ref() {
                    for t in tag {
                        let elem = match t {
                            models::Tag::Basic => data::TranscriptTag::Basic.into(),
                            models::Tag::EnsemblCanonical => {
                                data::TranscriptTag::EnsemblCanonical.into()
                            }
                            models::Tag::ManeSelect => data::TranscriptTag::ManeSelect.into(),
                            models::Tag::ManePlusClinical => {
                                data::TranscriptTag::ManePlusClinical.into()
                            }
                            models::Tag::RefSeqSelect => data::TranscriptTag::RefSeqSelect.into(),
                        };
                        tags.push(elem);
                    }
                }
                let models::Transcript {
                    protein,
                    start_codon,
                    stop_codon,
                    ..
                } = tx_model.clone();

                data_transcripts.push(data::Transcript {
                    id: tx_id.clone(),
                    gene_name: gene_symbol.expect("missing gene symbol"),
                    gene_id: hgnc.expect("missing HGNC ID"),
                    biotype,
                    tags,
                    protein,
                    start_codon,
                    stop_codon,
                    genome_alignments,
                });
            }
        }

        data_transcripts
    };
    tracing::info!(" ... done creating transcripts");

    trace_rss_now();

    // Build mapping of gene HGNC symbol to transcript IDs.
    tracing::info!("  Build gene symbol to transcript ID mapping ...");
    let gene_to_tx = transcript_ids_for_gene
        .into_iter()
        .map(|(gene_name, tx_ids)| data::GeneToTxId { gene_name, tx_ids })
        .collect::<Vec<_>>();
    tracing::info!(" ... done building gene symbol to transcript ID mapping");

    trace_rss_now();

    // Compose transcript database from transcripts and gene to transcript mapping.
    tracing::info!("  Composing transcript database ...");
    let tx_db = data::TranscriptDb {
        transcripts: data_transcripts,
        gene_to_tx,
    };
    tracing::info!(" ... done composing transcript database");

    trace_rss_now();

    // Compose the final transcript and sequence database.
    tracing::info!("  Constructing final tx and seq database ...");
    let tx_seq_db = data::TxSeqDatabase {
        tx_db: Some(tx_db),
        seq_db: Some(seq_db),
    };
    let mut buf = Vec::new();
    buf.reserve(tx_seq_db.encoded_len());
    tx_seq_db
        .encode(&mut buf)
        .map_err(|e| anyhow!("failed to encode: {}", e))?;
    tracing::info!("  ... done constructing final tx and seq database");

    trace_rss_now();

    // Write out the final transcript and sequence database.
    tracing::info!("  Writing out final database ...");
    // Open file and if necessary, wrap in a decompressor.
    let file = std::fs::File::create(path_out)
        .map_err(|e| anyhow!("failed to create file {}: {}", path_out.display(), e))?;
    let ext = path_out.extension().map(|s| s.to_str());
    let mut writer: Box<dyn Write> = if ext == Some(Some("gz")) {
        Box::new(flate2::write::GzEncoder::new(
            file,
            flate2::Compression::default(),
        ))
    } else if ext == Some(Some("zst")) {
        Box::new(zstd::Encoder::new(file, 0).map_err(|e| {
            anyhow!(
                "failed to open zstd enoder for {}: {}",
                path_out.display(),
                e
            )
        })?)
    } else {
        Box::new(file)
    };
    writer
        .write_all(&buf)
        .map_err(|e| anyhow!("failed to write to {}: {}", path_out.display(), e))?;
    tracing::info!("  ... done writing out final database");

    trace_rss_now();

    tracing::info!("... done with constructing protobuf file");
    Ok(())
}

/// Data as loaded from cdot after processing.
struct TranscriptData {
    pub genes: HashMap<String, models::Gene>,
    pub transcripts: HashMap<String, models::Transcript>,
    pub transcript_ids_for_gene: HashMap<String, Vec<String>>,
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
    tx_data: TranscriptData,
    max_genes: Option<u32>,
    gene_symbols: &Option<Vec<String>>,
    report_file: &mut File,
) -> Result<TranscriptData, anyhow::Error> {
    tracing::info!("Filtering transcripts ...");
    let start = Instant::now();
    let gene_symbols = gene_symbols.clone().unwrap_or_default();

    let TranscriptData {
        genes,
        transcripts,
        transcript_ids_for_gene,
    } = tx_data;

    // Potentially limit number of genes.
    let transcript_ids_for_gene = if let Some(max_genes) = max_genes {
        tracing::warn!("Limiting to {} genes!", max_genes);
        transcript_ids_for_gene
            .into_iter()
            .take(max_genes as usize)
            .collect()
    } else {
        transcript_ids_for_gene
    };

    // We keep track of the chosen transcript identifiers.
    let mut chosen = HashSet::new();
    // Filter map from gene symbol to Vec of chosen transcript identifiers.
    let transcript_ids_for_gene = {
        let mut tmp = HashMap::new();

        for (gene_symbol, tx_ids) in &transcript_ids_for_gene {
            // Skip transcripts where the gene symbol is not contained in `gene_symbols`.
            if !gene_symbols.is_empty() && !gene_symbols.contains(gene_symbol) {
                continue;
            }

            // Only select the highest version of each transcript.
            //
            // First, split off transcript versions from accessions and look for NM transcript.
            let mut seen_nm = false;
            let mut versioned: Vec<_> = tx_ids
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
                        writeln!(
                            report_file,
                            "skipped transcript {} because we have a later version already",
                            &full_ac
                        )?;
                        continue; // skip, already have later version
                    } else if ac.starts_with("NR_") && seen_nm {
                        writeln!(
                            report_file,
                            "skipped transcript {} because we have a NR transcript",
                            &full_ac
                        )?;
                        continue; // skip NR transcript as we have NM one
                    } else if ac.starts_with('X') {
                        writeln!(
                            report_file,
                            "skipped transcript {} because it is an XR/XM transcript",
                            &full_ac
                        )?;
                        continue; // skip XR/XM transcript
                    } else {
                        // Check transcript's CDS length for being multiple of 3 and skip unless it is.
                        let tx = transcripts
                            .get(&full_ac)
                            .expect("must exist; accession taken from map earlier");
                        if let Some(cds_start) = tx.start_codon {
                            let cds_end =
                                tx.stop_codon.expect("must be some if start_codon is some");
                            let cds_len = cds_end - cds_start;
                            if cds_len % 3 != 0 {
                                tracing::debug!("skipping transcript {} because its CDS length is not a multiple of 3", &full_ac);
                                writeln!(report_file, "skipped transcript {} because its CDS length {} is not a multiple of 3", &full_ac, cds_len)?;
                                continue;
                            }
                        }

                        // Otherwise, mark transcript as included by storing its accession.
                        next_tx_ids.push(full_ac.clone());
                        seen_ac.insert((ac.clone(), release));
                    }
                }
            }

            next_tx_ids.sort();
            next_tx_ids.dedup();
            chosen.extend(next_tx_ids.iter().cloned());

            if !next_tx_ids.is_empty() {
                tmp.insert(gene_symbol.clone(), next_tx_ids);
            } else {
                writeln!(
                    report_file,
                    "skipped gene {} because we have no transcripts left",
                    gene_symbol
                )?;
            }
        }

        tmp
    };

    let transcripts: HashMap<_, _> = transcripts
        .into_iter()
        .filter(|(tx_id, _)| chosen.contains(tx_id))
        .collect();
    tracing::debug!(
        "  => {} transcripts left",
        transcripts.len().separate_with_commas()
    );
    writeln!(report_file, "total transcripts\t{}", transcripts.len())?;

    let genes: HashMap<_, _> = genes
        .into_iter()
        .filter(|(gene_id, _)| transcript_ids_for_gene.contains_key(gene_id))
        .collect();
    tracing::debug!("  => {} genes left", genes.len().separate_with_commas());

    tracing::info!("... done filtering transcripts in {:?}", start.elapsed());
    Ok(TranscriptData {
        genes,
        transcripts,
        transcript_ids_for_gene,
    })
}

/// Create file-backed `SeqRepo`.
fn open_seqrepo(args: &Args) -> Result<SeqRepo, anyhow::Error> {
    tracing::info!("Opening seqrepo...");
    let start = Instant::now();
    let seqrepo = PathBuf::from(&args.path_seqrepo_instance);
    let path = seqrepo
        .parent()
        .ok_or(anyhow::anyhow!(
            "Could not get parent from {:?}",
            &args.path_seqrepo_instance
        ))?
        .to_str()
        .unwrap()
        .to_string();
    let instance = seqrepo
        .file_name()
        .ok_or(anyhow::anyhow!(
            "Could not get basename from {:?}",
            &args.path_seqrepo_instance
        ))?
        .to_str()
        .unwrap()
        .to_string();
    let seqrepo = SeqRepo::new(path, &instance)?;
    tracing::info!("... seqrepo opened in {:?}", start.elapsed());
    Ok(seqrepo)
}

/// Load the cdot JSON files.
fn load_cdot_files(args: &Args, report_file: &mut File) -> Result<TranscriptData, anyhow::Error> {
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
            args.genome_release,
            report_file,
        )?;
    }
    tracing::info!(
        "... done loading cdot JSON files in {:?} -- #genes = {}, #transcripts = {}, #transcript_ids_for_gene = {}",
        start.elapsed(),
        genes.len().separate_with_commas(),
        transcripts.len().separate_with_commas(),
        transcript_ids_for_gene.len().separate_with_commas()
    );
    writeln!(
        report_file,
        "total genes\t{}\ntotal transcripts\t{}",
        transcripts.len(),
        transcript_ids_for_gene.len()
    )?;

    Ok(TranscriptData {
        genes,
        transcripts,
        transcript_ids_for_gene,
    })
}

/// Main entry point for `db create txs` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let mut report_file = File::create(format!("{}.report", args.path_out.display()))?;
    tracing::info!(
        "Building transcript and sequence database file\ncommon args: {:#?}\nargs: {:#?}",
        common,
        args
    );

    // Open seqrepo,
    let seqrepo = open_seqrepo(args)?;
    // then load cdot files,
    let tx_data = load_cdot_files(args, &mut report_file)?;
    // then remove redundant onces, and
    let tx_data = filter_transcripts(tx_data, args.max_txs, &args.gene_symbols, &mut report_file)?;
    // finally build protobuf file.
    build_protobuf(
        &args.path_out,
        seqrepo,
        tx_data,
        common.verbose.is_silent(),
        &mut report_file,
    )?;

    tracing::info!("Done building transcript and sequence database file");
    Ok(())
}

#[cfg(test)]
pub mod test {
    use std::collections::HashMap;
    use std::fs::File;
    use std::path::{Path, PathBuf};

    use clap_verbosity_flag::Verbosity;
    use pretty_assertions::assert_eq;
    use temp_testdir::TempDir;

    use crate::common::{Args as CommonArgs, GenomeRelease};
    use crate::db::create::txs::TranscriptData;

    use super::{filter_transcripts, load_and_extract, run, Args};

    #[test]
    fn filter_transcripts_brca1() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();
        let mut report_file = File::create(tmp_dir.join("report"))?;

        let mut genes = HashMap::new();
        let mut transcripts = HashMap::new();
        let mut transcript_ids_for_gene = HashMap::new();
        load_and_extract(
            Path::new("tests/data/db/create/txs/cdot-0.2.12.refseq.grch37_grch38.brca1_opa1.json"),
            &mut transcript_ids_for_gene,
            &mut genes,
            &mut transcripts,
            GenomeRelease::Grch37,
            &mut report_file,
        )?;

        let tx_data = TranscriptData {
            genes,
            transcripts,
            transcript_ids_for_gene,
        };

        assert_eq!(
            &tx_data
                .transcript_ids_for_gene
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
            ]
        );
        let filtered = filter_transcripts(tx_data, None, &None, &mut report_file)?;
        assert_eq!(
            &filtered
                .transcript_ids_for_gene
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

    #[test]
    fn run_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            path_out: tmp_dir.join("out.bin.zst"),
            path_cdot_json: vec![PathBuf::from(
                "tests/data/db/create/txs/cdot-0.2.12.refseq.grch37_grch38.brca1_opa1.json",
            )],
            path_seqrepo_instance: PathBuf::from("tests/data/db/create/txs/latest"),
            genome_release: GenomeRelease::Grch38,
            max_txs: None,
            gene_symbols: None,
        };

        run(&common_args, &args)?;

        Ok(())
    }
}
