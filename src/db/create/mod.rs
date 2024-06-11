//! Transcript database.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::{io::Write, path::PathBuf, time::Instant};

use anyhow::{anyhow, Error};
use clap::Parser;
use hgvs::data::cdot::json::models;
use hgvs::data::cdot::json::models::{Gene, Tag, Transcript};
use hgvs::sequences::{translate_cds, TranslationTable};
use indexmap::{IndexMap, IndexSet};
use indicatif::{ProgressBar, ProgressStyle};
use itertools::{Either, Itertools};
use prost::Message;
use seqrepo::{AliasOrSeqId, Interface, SeqRepo};
use serde::Serialize;
use thousands::Separable;

use crate::common::{trace_rss_now, GenomeRelease};

lazy_static::lazy_static! {
    /// Progress bar style to use.
    pub static ref PROGRESS_STYLE: ProgressStyle = ProgressStyle::with_template(
        "[{elapsed_precise}] [{wide_bar:.cyan/blue}] {human_pos}/{human_len} ({eta})",
    )
    .unwrap();
}

/// Mitochondrial accessions.
const MITOCHONDRIAL_ACCESSIONS: &[&str] = &["NC_012920.1"];

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
    /// Path to TSV file for label transfer of transcripts.  Columns are
    /// transcript id (without version), (unused) gene symbol, and label.
    #[arg(long)]
    pub path_mane_txs_tsv: Option<PathBuf>,
    /// Maximal number of transcripts to process.
    #[arg(long)]
    pub max_txs: Option<u32>,
    /// Limit transcript database to the following HGNC symbols.  Useful for
    /// building test databases.
    #[arg(long)]
    pub gene_symbols: Option<Vec<String>>,
}

/// Helper struct for parsing the label TSV file.
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize)]
struct LabelEntry {
    /// Transcript identifier without version.
    transcript_id: String,
    /// Gene symbol (unused).
    _gene_symbol: String,
    /// Label to transfer.
    label: String,
}

/// Load and extract from cdot JSON.
#[allow(clippy::too_many_arguments)]
fn load_and_extract_into(
    json_path: &Path,
    label_tsv_path: &Option<&Path>,
    transcript_ids_for_gene: &mut IndexMap<String, Vec<String>>,
    genes: &mut IndexMap<String, Gene>,
    transcripts: &mut IndexMap<String, Transcript>,
    genome_release: GenomeRelease,
    cdot_version: &mut String,
    report_file: &mut impl Write,
    mt_tx_ids: &mut IndexSet<String>,
) -> Result<(), Error> {
    writeln!(
        report_file,
        "{}",
        serde_json::to_string(
            &serde_json::json!({"source": json_path, "genome_release": genome_release})
        )?
    )?;
    writeln!(
        report_file,
        "{}",
        serde_json::to_string(
            &serde_json::json!({"source": json_path, "label_tsv_path": label_tsv_path})
        )?
    )?;

    let txid_to_label = label_tsv_path.map(txid_to_label).transpose()?;
    let source = json_path.to_str().expect("Invalid path");
    let (c_genes, c_txs, c_version) = load_cdot_transcripts(json_path)?;
    *cdot_version = c_version;

    // Count number of MANE Select and MANE Plus Clinical transcripts, collect
    // chrMT gene names.
    let (genes_chrmt, n_mane_select, n_mane_plus_clinical) =
        gather_transcript_stats(mt_tx_ids, &c_txs);
    writeln!(
        report_file,
        "{}",
        serde_json::to_string(
            &serde_json::json!({"source": source, "n_chr_mt": genes_chrmt.len(), "n_mane_select": n_mane_select, "n_mane_plus_clinical": n_mane_plus_clinical})
        )?
    )?;

    let (keep, discard) = filter_genes(&source, &c_genes, &genes_chrmt);
    writeln!(
        report_file,
        "{}",
        serde_json::to_string(
            &serde_json::json!({"source": source, "kept": keep.len(), "discarded": discard.len()})
        )?
    )?;

    for (_gene_id, gene) in keep {
        let hgnc_id = format!("HGNC:{}", gene.hgnc.as_ref().unwrap());
        transcript_ids_for_gene.entry(hgnc_id.clone()).or_default();
        genes.insert(hgnc_id, gene.clone());
    }
    for d in discard {
        writeln!(report_file, "{}", serde_json::to_string(&d)?)?;
    }
    writeln!(
        report_file,
        "{}",
        serde_json::to_string(
            &serde_json::json!({"source": source, "total_genes": c_genes.len(), "genes_kept": genes.len()})
        )?,
    )?;
    let total_transcripts = c_txs.len();
    process_transcripts(
        source,
        transcript_ids_for_gene,
        genes,
        transcripts,
        genome_release,
        c_txs,
        txid_to_label,
        report_file,
    );
    writeln!(
        report_file,
        "{}",
        serde_json::to_string(
            &serde_json::json!({"source": source, "total_transcripts": total_transcripts, "transcripts_kept": transcripts.len()})
        )?,
    )?;
    Ok(())
}

fn process_transcripts(
    source: &str,
    transcript_ids_for_gene: &mut IndexMap<String, Vec<String>>,
    genes: &mut IndexMap<String, Gene>,
    transcripts: &mut IndexMap<String, Transcript>,
    genome_release: GenomeRelease,
    c_txs: IndexMap<String, Transcript>,
    txid_to_label: Option<IndexMap<String, Vec<Tag>>>,
    report_file: &mut impl Write,
) {
    let missing_hgnc =
        |tx: &Transcript| -> bool { tx.hgnc.is_none() || tx.hgnc.as_ref().unwrap().is_empty() };
    let deselected_gene = |tx: &Transcript| -> bool {
        !genes.contains_key(&format!("HGNC:{}", tx.hgnc.as_ref().unwrap()))
    };
    let empty_genome_builds = |tx: &Transcript| -> bool { tx.genome_builds.is_empty() };
    let filters: [(&dyn Fn(&Transcript) -> bool, Reason); 3] = [
        (&missing_hgnc, Reason::MissingHgncId),
        (&deselected_gene, Reason::DeselectedGene),
        (&empty_genome_builds, Reason::EmptyGenomeBuilds),
    ];
    c_txs
        .values()
        .map(|tx| Transcript {
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
            for (filter, reason) in &filters {
                if filter(tx) {
                    let d = Discard {
                        source: source.into(),
                        kind: GeneOrTranscript::Transcript,
                        reason: *reason,
                        id: tx.id.clone(),
                        gene_name: tx.gene_name.clone(),
                    };
                    writeln!(
                        report_file,
                        "{}",
                        serde_json::to_string(&d).expect("Failed serializing")
                    )
                    .expect("Failed writing to report file");
                    return false;
                }
            }
            true
        })
        .for_each(|tx| {
            let hgnc_id = &format!("HGNC:{}", tx.hgnc.as_ref().unwrap());
            transcript_ids_for_gene
                .get_mut(hgnc_id)
                .unwrap_or_else(|| panic!("tx {:?} for unknown gene {:?}", tx.id, hgnc_id))
                .push(tx.id.clone());
            // build output transcripts
            let mut tx_out = tx.clone();
            // transfer MANE-related labels from TSV file
            if let Some(txid_to_tags) = txid_to_label.as_ref() {
                let tx_id_no_version = tx.id.split('.').next().unwrap();
                if let Some(tags) = txid_to_tags.get(tx_id_no_version) {
                    tx_out.genome_builds.iter_mut().for_each(|(_, alignment)| {
                        if let Some(alignment_tag) = &mut alignment.tag {
                            alignment_tag.extend(tags.iter().cloned());
                            alignment_tag.sort();
                            alignment_tag.dedup();
                        }
                    });
                }
            }
            // fix coding mitochondrial transcripts that have a CDS that is not a multiple of 3
            if let Some(cds_start) = tx_out.start_codon {
                let cds_end = tx_out
                    .stop_codon
                    .expect("must be some if start_codon is some");
                let cds_len = cds_end - cds_start;
                if cds_len % 3 != 0 {
                    assert_eq!(
                        tx.genome_builds.len(),
                        1,
                        "only one genome build expected at this point"
                    );
                    let gb = tx_out.genome_builds.iter_mut().next().unwrap().1;
                    if MITOCHONDRIAL_ACCESSIONS.contains(&gb.contig.as_ref()) {
                        assert_eq!(gb.exons.len(), 1, "only single-exon genes assumed on chrMT");
                        let delta = 3 - cds_len % 3;
                        tx_out.stop_codon = Some(cds_end + delta);
                        let exon = gb.exons.iter_mut().next().unwrap();
                        exon.alt_cds_end_i += delta;
                        exon.cigar.push_str(&format!("{}I", delta));
                    }
                }
            }
            // finally, insert into transcripts
            transcripts.insert(tx.id.clone(), tx_out);
        });
}

#[derive(Debug, Clone, Copy, Serialize)]
enum Reason {
    MissingHgncId,
    // NotMTandMissingMapLocation,
    DeselectedGene,
    EmptyGenomeBuilds,
    OldVersion,
    UseNmTranscriptInsteadOfNr,
    PredictedTranscript,
    InvalidCdsLength,
    NoTranscript,
    MissingSequence,
    CdsEndAfterSequenceEnd,
    MissingStopCodon,
}

#[derive(Debug, Clone, Copy, Serialize)]
enum GeneOrTranscript {
    Gene,
    Transcript,
}

#[derive(Debug, Clone, Serialize)]
struct Discard {
    source: String,
    kind: GeneOrTranscript,
    reason: Reason,
    id: String,
    gene_name: Option<String>,
}

fn filter_genes(
    source: &str,
    c_genes: &IndexMap<String, Gene>,
    _genes_chrmt: &IndexSet<String>,
) -> (Vec<(String, Gene)>, Vec<Discard>) {
    let missing_hgnc = |_gene_id: &str, gene: &Gene| -> bool {
        gene.hgnc.is_none() || gene.hgnc.as_ref().unwrap().is_empty()
    };
    // let not_mitochondrion_and_missing_map_location = |gene_id: &str, gene: &Gene| -> bool {
    //     !genes_chrmt.contains(gene_id)
    //         && (gene.map_location.is_none() || gene.map_location.as_ref().unwrap().is_empty())
    // };
    let filters: [(&dyn Fn(&str, &Gene) -> bool, Reason); 1] = [
        (&missing_hgnc, Reason::MissingHgncId),
        // (
        //     &not_mitochondrion_and_missing_map_location,
        //     Reason::NotMTandMissingMapLocation,
        // ),
    ];

    c_genes.iter().partition_map(|(gene_id, gene)| {
        for (filter, reason) in &filters {
            if filter(gene_id, gene) {
                return Either::Right(Discard {
                    source: source.into(),
                    kind: GeneOrTranscript::Gene,
                    reason: *reason,
                    id: gene_id.to_string(),
                    gene_name: gene.gene_symbol.clone(),
                });
            }
        }
        return Either::Left((gene_id.to_string(), gene.clone()));
    })
}

fn gather_transcript_stats(
    mt_tx_ids: &mut IndexSet<String>,
    c_txs: &IndexMap<String, Transcript>,
) -> (IndexSet<String>, i32, i32) {
    let mut genes_chrmt = IndexSet::new();
    let mut n_mane_select = 0;
    let mut n_mane_plus_clinical = 0;
    for tx in c_txs.values() {
        let mut is_mane_select = false;
        let mut is_mane_plus_clinical = false;
        for gb in tx.genome_builds.values() {
            if MITOCHONDRIAL_ACCESSIONS.contains(&gb.contig.as_str()) {
                genes_chrmt.insert(tx.gene_version.clone());
                mt_tx_ids.insert(tx.id.clone());
            }
            if let Some(tag) = &gb.tag {
                if tag.contains(&models::Tag::ManeSelect) {
                    is_mane_select = true;
                }
                if tag.contains(&models::Tag::ManePlusClinical) {
                    is_mane_plus_clinical = true;
                }
            }
        }
        if is_mane_select {
            n_mane_select += 1;
        }
        if is_mane_plus_clinical {
            n_mane_plus_clinical += 1;
        }
    }
    (genes_chrmt, n_mane_select, n_mane_plus_clinical)
}

fn txid_to_label(label_tsv_path: impl AsRef<Path>) -> Result<IndexMap<String, Vec<Tag>>, Error> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_path(label_tsv_path.as_ref())?;

    rdr.deserialize()
        .map(|result| {
            result
                .map(|entry: LabelEntry| {
                    (
                        entry.transcript_id,
                        entry
                            .label
                            .split(',')
                            .map(models::str_to_tag)
                            .collect::<Vec<_>>(),
                    )
                })
                .map_err(anyhow::Error::from)
        })
        .collect()
}
fn load_cdot_transcripts(
    json_path: &Path,
) -> Result<(IndexMap<String, Gene>, IndexMap<String, Transcript>, String), Error> {
    let models::Container {
        genes: c_genes,
        transcripts: c_txs,
        cdot_version: c_version,
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
    Ok((c_genes, c_txs, c_version))
}

/// Perform protobuf file construction.
///
/// This can be done by simply converting the models from ``hvs-rs`` to the prost generated data structures.
fn build_protobuf(
    path_out: &Path,
    seqrepo: SeqRepo,
    mt_tx_ids: IndexSet<String>,
    tx_data: TranscriptData,
    is_silent: bool,
    genome_release: GenomeRelease,
    report_file: &mut impl Write,
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
    let mut tx_skipped_noseq = IndexSet::new(); // skipped because of missing sequence
    let mut tx_skipped_nostop = IndexSet::new(); // skipped because of missing stop codon
    let seq_db = {
        // Insert into protobuf and keep track of pointers in `Vec`s.
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
            let is_mt = mt_tx_ids.contains(tx_id);
            let namespace: Option<String> = if tx_id.starts_with("ENST") {
                Some(String::from("Ensembl"))
            } else {
                Some(String::from("NCBI"))
            };
            let res_seq = seqrepo.fetch_sequence(&AliasOrSeqId::Alias {
                value: tx_id.clone(),
                namespace,
            });
            let seq = if let Ok(seq) = res_seq {
                // Append poly-A for chrMT transcripts (which are from ENSEMBL).
                // This also potentially fixes the stop codon.
                if is_mt {
                    let mut seq = seq.into_bytes();
                    seq.extend_from_slice(b"A".repeat(300).as_slice());
                    String::from_utf8(seq).expect("must be valid UTF-8")
                } else {
                    seq
                }
            } else {
                let d = Discard {
                    source: "protobuf".into(),
                    kind: GeneOrTranscript::Transcript,
                    reason: Reason::MissingSequence,
                    id: tx_id.clone(),
                    gene_name: None,
                };
                writeln!(report_file, "{}", serde_json::to_string(&d)?)?;
                tx_skipped_noseq.insert(tx_id.clone());
                continue;
            };

            // Skip transcript if it is coding and the translated CDS does not have a stop codon.
            if let Some(cds_start) = tx.start_codon {
                let cds_start = cds_start as usize;
                let cds_end = tx.stop_codon.expect("must be some if start_codon is some") as usize;
                if cds_end > seq.len() {
                    let d = Discard {
                        source: "protobuf".into(),
                        kind: GeneOrTranscript::Transcript,
                        reason: Reason::CdsEndAfterSequenceEnd, // (cds_end, seq.len())
                        id: tx_id.clone(),
                        gene_name: None,
                    };
                    writeln!(report_file, "{}", serde_json::to_string(&d)?)?;
                    continue;
                }
                let tx_seq_to_translate = &seq[cds_start..cds_end];
                let aa_sequence =
                    translate_cds(tx_seq_to_translate, true, "*", TranslationTable::Standard)?;
                if (!is_mt && !aa_sequence.ends_with('*')) || (is_mt && !aa_sequence.contains('*'))
                {
                    let d = Discard {
                        source: "protobuf".into(),
                        kind: GeneOrTranscript::Transcript,
                        reason: Reason::MissingStopCodon,
                        id: tx_id.clone(),
                        gene_name: None,
                    };
                    writeln!(report_file, "{}", serde_json::to_string(&d)?)?;
                    tx_skipped_nostop.insert(tx_id.clone());
                    continue;
                }
            }

            // Register sequence into protobuf.
            aliases.push(tx_id.clone());
            aliases_idx.push(seqs.len() as u32);
            seqs.push(seq.clone());
        }
        pb.finish_and_clear();
        // Finalize by creating `SequenceDb`.
        crate::pbs::txs::SequenceDb {
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
                let d = Discard {
                    source: "protobuf".into(),
                    kind: GeneOrTranscript::Gene,
                    reason: Reason::NoTranscript,
                    id: gene_symbol.clone(),
                    gene_name: gene.gene_symbol.clone(),
                };
                writeln!(report_file, "{}", serde_json::to_string(&d)?)?;
                continue;
            }

            // ... for each transcript of the gene ...
            for tx_id in tx_ids {
                let mut tags: Vec<i32> = Vec::new();
                let tx_model = transcripts
                    .get(tx_id)
                    .unwrap_or_else(|| panic!("No transcript model for id {:?}", tx_id));
                // ... build genome alignment for selected:
                let mut genome_alignments = Vec::new();
                for (genome_build, alignment) in &tx_model.genome_builds {
                    // obtain basic properties
                    let genome_build = match genome_build.as_ref() {
                        "GRCh37" => crate::pbs::txs::GenomeBuild::Grch37,
                        "GRCh38" => crate::pbs::txs::GenomeBuild::Grch38,
                        _ => panic!("Unknown genome build {:?}", genome_build),
                    };
                    let models::GenomeAlignment {
                        contig,
                        cds_start,
                        cds_end,
                        ..
                    } = alignment.clone();
                    let strand = match alignment.strand {
                        models::Strand::Plus => crate::pbs::txs::Strand::Plus,
                        models::Strand::Minus => crate::pbs::txs::Strand::Minus,
                    };
                    if let Some(tag) = alignment.tag.as_ref() {
                        for t in tag {
                            let elem = match t {
                                models::Tag::Basic => crate::pbs::txs::TranscriptTag::Basic.into(),
                                models::Tag::EnsemblCanonical => {
                                    crate::pbs::txs::TranscriptTag::EnsemblCanonical.into()
                                }
                                models::Tag::ManeSelect => {
                                    crate::pbs::txs::TranscriptTag::ManeSelect.into()
                                }
                                models::Tag::ManePlusClinical => {
                                    crate::pbs::txs::TranscriptTag::ManePlusClinical.into()
                                }
                                models::Tag::RefSeqSelect => {
                                    crate::pbs::txs::TranscriptTag::RefSeqSelect.into()
                                }
                            };
                            if !tags.contains(&elem) {
                                tags.push(elem);
                            }
                        }
                    }
                    // Look into any "note" string for a selenoprotein marker and
                    // add this as a tag.
                    if let Some(note) = alignment.note.as_ref() {
                        let needle = "UGA stop codon recoded as selenocysteine";
                        if note.contains(needle) {
                            tags.push(crate::pbs::txs::TranscriptTag::Selenoprotein.into());
                        }
                    }
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
                            crate::pbs::txs::ExonAlignment {
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
                    genome_alignments.push(crate::pbs::txs::GenomeAlignment {
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
                    crate::pbs::txs::TranscriptBiotype::Coding.into()
                } else {
                    crate::pbs::txs::TranscriptBiotype::NonCoding.into()
                };
                let models::Transcript {
                    protein,
                    start_codon,
                    stop_codon,
                    ..
                } = tx_model.clone();

                tags.sort();
                tags.dedup();

                data_transcripts.push(crate::pbs::txs::Transcript {
                    id: tx_id.clone(),
                    gene_symbol: gene_symbol.expect("missing gene symbol"),
                    gene_id: format!("HGNC:{}", hgnc.expect("missing HGNC ID")),
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
        .map(|(gene_id, tx_ids)| crate::pbs::txs::GeneToTxId { gene_id, tx_ids })
        .collect::<Vec<_>>();
    tracing::info!(" ... done building gene symbol to transcript ID mapping");

    trace_rss_now();

    // Compose transcript database from transcripts and gene to transcript mapping.
    tracing::info!("  Composing transcript database ...");
    let tx_db = crate::pbs::txs::TranscriptDb {
        transcripts: data_transcripts,
        gene_to_tx,
    };
    tracing::info!(" ... done composing transcript database");

    trace_rss_now();

    // Compose the final transcript and sequence database.
    tracing::info!("  Constructing final tx and seq database ...");
    let tx_seq_db = crate::pbs::txs::TxSeqDatabase {
        tx_db: Some(tx_db),
        seq_db: Some(seq_db),
        version: Some(crate::common::version().to_string()),
        genome_release: Some(genome_release.name()),
    };
    let mut buf = Vec::with_capacity(tx_seq_db.encoded_len());
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
        Box::new(
            zstd::Encoder::new(file, 0)
                .map_err(|e| {
                    anyhow!(
                        "failed to open zstd encoder for {}: {}",
                        path_out.display(),
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
        .map_err(|e| anyhow!("failed to write to {}: {}", path_out.display(), e))?;
    tracing::info!("  ... done writing out final database");

    trace_rss_now();

    tracing::info!("... done with constructing protobuf file");
    Ok(())
}

/// Data as loaded from cdot after processing.
#[derive(Debug)]
struct TranscriptData {
    pub genes: IndexMap<String, Gene>,
    pub transcripts: IndexMap<String, Transcript>,
    pub transcript_ids_for_gene: IndexMap<String, Vec<String>>,
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
    report_file: &mut impl Write,
) -> Result<TranscriptData, anyhow::Error> {
    tracing::info!("Filtering transcripts ...");
    let start = Instant::now();
    let selected_hgnc_ids = gene_symbols.as_ref().map(|gene_symbols| {
        let symbol_to_hgnc: IndexMap<_, _> =
            IndexMap::from_iter(tx_data.genes.iter().flat_map(|(hgnc_id, g)| {
                g.gene_symbol
                    .as_ref()
                    .map(|gene_symbol| (gene_symbol.clone(), hgnc_id.clone()))
            }));
        let result = gene_symbols
            .iter()
            .map(|s| symbol_to_hgnc.get(s).unwrap_or(s).clone())
            .collect::<Vec<_>>();
        tracing::info!("Will limit to {:?}", &result);
        result
    });

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
    let mut chosen = IndexSet::new();
    // Filter map from gene symbol to Vec of chosen transcript identifiers.
    let transcript_ids_for_gene = {
        let mut tmp = IndexMap::new();

        for (hgnc_id, tx_ids) in &transcript_ids_for_gene {
            // Skip transcripts where the gene symbol is not contained in `selected_hgnc_ids`.
            if !selected_hgnc_ids
                .as_ref()
                .map(|ids| ids.contains(hgnc_id))
                .unwrap_or(true)
            {
                tracing::trace!("skipping {} / {:?}, because not selected", hgnc_id, tx_ids);
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
            versioned.sort_unstable_by_key(|(_, v)| std::cmp::Reverse(*v));

            // Build `next_tx_ids`.
            let mut seen_ac = IndexSet::new();
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
                        let d = Discard {
                            source: "aggregated_cdot".into(),
                            kind: GeneOrTranscript::Transcript,
                            reason: Reason::OldVersion,
                            id: full_ac.clone(),
                            gene_name: None,
                        };
                        writeln!(report_file, "{}", serde_json::to_string(&d)?)?;
                        continue; // skip, already have later version
                    } else if ac.starts_with("NR_") && seen_nm {
                        let d = Discard {
                            source: "aggregated_cdot".into(),
                            kind: GeneOrTranscript::Transcript,
                            reason: Reason::UseNmTranscriptInsteadOfNr,
                            id: full_ac.clone(),
                            gene_name: None,
                        };
                        writeln!(report_file, "{}", serde_json::to_string(&d)?)?;
                        continue; // skip NR transcript as we have NM one
                    } else if ac.starts_with('X') {
                        let d = Discard {
                            source: "aggregated_cdot".into(),
                            kind: GeneOrTranscript::Transcript,
                            reason: Reason::PredictedTranscript,
                            id: full_ac.clone(),
                            gene_name: None,
                        };
                        writeln!(report_file, "{}", serde_json::to_string(&d)?)?;
                        continue; // skip XR/XM transcript
                    } else {
                        // Check transcript's CDS length for being multiple of 3 and skip unless
                        // it is.
                        //
                        // Note that the chrMT transcripts have been fixed earlier already to
                        // accomodate for how they are fixed by poly-A tailing.
                        let tx = transcripts
                            .get(&full_ac)
                            .expect("must exist; accession taken from map earlier");
                        if let Some(cds_start) = tx.start_codon {
                            let cds_end =
                                tx.stop_codon.expect("must be some if start_codon is some");
                            let cds_len = cds_end - cds_start;
                            if cds_len % 3 != 0 {
                                let d = Discard {
                                    source: "aggregated_cdot".into(),
                                    kind: GeneOrTranscript::Transcript,
                                    reason: Reason::InvalidCdsLength,
                                    id: full_ac.clone(),
                                    gene_name: None,
                                };
                                writeln!(report_file, "{}", serde_json::to_string(&d)?)?;
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
                tmp.insert(hgnc_id.clone(), next_tx_ids);
            } else {
                let d = Discard {
                    source: "aggregated_cdot".into(),
                    kind: GeneOrTranscript::Gene,
                    reason: Reason::NoTranscript,
                    id: hgnc_id.clone(),
                    gene_name: None,
                };
                writeln!(report_file, "{}", serde_json::to_string(&d)?)?;
            }
        }

        tmp
    };

    let transcripts: IndexMap<_, _> = transcripts
        .into_iter()
        .filter(|(tx_id, _)| chosen.contains(tx_id))
        .collect();
    tracing::debug!(
        "  => {} transcripts left",
        transcripts.len().separate_with_commas()
    );
    writeln!(
        report_file,
        "{}",
        serde_json::to_string(
            &serde_json::json!({"source": "aggregated_cdot", "total_transcripts": transcripts.len()})
        )?
    )?;

    let genes: indexmap::IndexMap<_, _> = genes
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
fn load_cdot_files(
    args: &Args,
    report_file: &mut impl Write,
) -> Result<(IndexSet<String>, TranscriptData), Error> {
    tracing::info!("Loading cdot JSON files ...");
    let start = Instant::now();
    let mut genes = IndexMap::new();
    let mut transcripts = IndexMap::new();
    let mut transcript_ids_for_gene = IndexMap::new();
    let mut cdot_version = String::new();
    let mut mt_tx_ids = IndexSet::new();
    for json_path in &args.path_cdot_json {
        load_and_extract_into(
            json_path,
            &args.path_mane_txs_tsv.as_ref().map(|p| p.as_ref()),
            &mut transcript_ids_for_gene,
            &mut genes,
            &mut transcripts,
            args.genome_release,
            &mut cdot_version,
            report_file,
            &mut mt_tx_ids,
        )?;
    }
    tracing::info!(
        "... done loading cdot JSON files in {:?} -- #genes = {}, #transcripts = {}, #transcript_ids_for_gene = {}",
        start.elapsed(),
        genes.len().separate_with_underscores(),
        transcripts.len().separate_with_underscores(),
        transcript_ids_for_gene.len().separate_with_underscores()
    );
    writeln!(
        report_file,
        "{}",
        serde_json::to_string(
            &serde_json::json!({"source": "aggregated_cdot", "total_genes": transcripts.len(), "total_transcripts": transcripts.len(), "total_transcript_ids_for_gene": transcript_ids_for_gene.len()})
        )?
    )?;

    Ok((
        mt_tx_ids,
        TranscriptData {
            genes,
            transcripts,
            transcript_ids_for_gene,
        },
    ))
}

/// Main entry point for `db create txs` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let mut report_file =
        File::create(format!("{}.report.jsonl", args.path_out.display())).map(BufWriter::new)?;
    tracing::info!(
        "Building transcript and sequence database file\ncommon args: {:#?}\nargs: {:#?}",
        common,
        args
    );

    // Open seqrepo,
    let seqrepo = open_seqrepo(args)?;
    // then load cdot files,
    let (mt_tx_ids, tx_data) = load_cdot_files(args, &mut report_file)?;
    // then remove redundant onces, and
    let tx_data = filter_transcripts(tx_data, args.max_txs, &args.gene_symbols, &mut report_file)?;
    // finally build protobuf file.
    build_protobuf(
        &args.path_out,
        seqrepo,
        mt_tx_ids,
        tx_data,
        common.verbose.is_silent(),
        args.genome_release,
        &mut report_file,
    )?;

    tracing::info!("Done building transcript and sequence database file");
    Ok(())
}

#[cfg(test)]
pub mod test {
    use std::fs::File;
    use std::path::{Path, PathBuf};

    use clap_verbosity_flag::Verbosity;
    use temp_testdir::TempDir;

    use crate::common::{Args as CommonArgs, GenomeRelease};
    use crate::db::create::TranscriptData;
    use crate::db::dump;

    use super::{filter_transcripts, load_and_extract_into, run, Args};

    #[test]
    fn filter_transcripts_brca1() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();
        let mut report_file = File::create(tmp_dir.join("report"))?;

        let mut genes = indexmap::IndexMap::new();
        let mut transcripts = indexmap::IndexMap::new();
        let mut transcript_ids_for_gene = indexmap::IndexMap::new();
        let mut cdot_version = String::new();
        let path_tsv = Path::new("tests/data/db/create/txs/txs_main.tsv");
        let mut mt_tx_ids = indexmap::IndexSet::new();
        load_and_extract_into(
            Path::new("tests/data/db/create/txs/cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json"),
            &Some(path_tsv),
            &mut transcript_ids_for_gene,
            &mut genes,
            &mut transcripts,
            GenomeRelease::Grch37,
            &mut cdot_version,
            &mut report_file,
            &mut mt_tx_ids,
        )?;

        let tx_data = TranscriptData {
            genes,
            transcripts,
            transcript_ids_for_gene,
        };

        eprintln!("{:#?}", &tx_data.transcript_ids_for_gene);
        insta::assert_yaml_snapshot!(tx_data
            .transcript_ids_for_gene
            .get("HGNC:1100")
            .unwrap()
            .iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>());

        let filtered = filter_transcripts(tx_data, None, &None, &mut report_file)?;
        insta::assert_yaml_snapshot!(filtered
            .transcript_ids_for_gene
            .get("HGNC:1100")
            .unwrap()
            .iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>());

        insta::assert_snapshot!(&cdot_version);

        Ok(())
    }

    #[test]
    fn run_smoke_brca1_opa1() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            path_out: tmp_dir.join("out.bin.zst"),
            path_cdot_json: vec![PathBuf::from(
                "tests/data/db/create/txs/cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json",
            )],
            path_mane_txs_tsv: Some(PathBuf::from("tests/data/db/create/txs/txs_main.tsv")),
            path_seqrepo_instance: PathBuf::from("tests/data/db/create/txs/latest"),
            genome_release: GenomeRelease::Grch38,
            max_txs: None,
            gene_symbols: None,
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

    #[test]
    fn run_smoke_selenoproteins() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();

        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 1),
        };
        let args = Args {
            path_out: tmp_dir.join("out.bin.zst"),
            path_cdot_json: vec![PathBuf::from(
                "tests/data/db/create/seleonoproteins/cdot-0.2.22.refseq.grch38.selenon.json",
            )],
            path_mane_txs_tsv: Some(PathBuf::from("tests/data/db/create/txs/txs_main.tsv")),
            path_seqrepo_instance: PathBuf::from("tests/data/db/create/seleonoproteins/latest"),
            genome_release: GenomeRelease::Grch38,
            max_txs: None,
            gene_symbols: None,
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
            path_out: tmp_dir.join("out.bin.zst"),
            path_cdot_json: vec![PathBuf::from(
                "tests/data/db/create/mitochondrial/cdot-0.2.23.ensembl.chrMT.grch37.gff3.json",
            )],
            path_mane_txs_tsv: None,
            path_seqrepo_instance: PathBuf::from("tests/data/db/create/mitochondrial/latest"),
            genome_release: GenomeRelease::Grch37,
            max_txs: None,
            gene_symbols: None,
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
