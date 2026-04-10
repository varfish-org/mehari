use crate::db::create::cdot_models;
use crate::db::create::models::{GeneId, TranscriptId, TranscriptLoader};
use anyhow::Error;
use hgvs::data::cdot::json::models::{Gene, GenomeAlignment, Transcript};
use indexmap::IndexMap;
use noodles::gff::feature::record::Strand;
use noodles::gff::feature::record_buf::attributes::field::tag;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Load and extract from standard generic GFF3 using noodles::gff.
pub fn load_gff3(loader: &mut TranscriptLoader, path: impl AsRef<Path>) -> Result<(), Error> {
    let file = File::open(path.as_ref())?;
    let reader: Box<dyn std::io::Read> = if path.as_ref().extension().is_some_and(|e| e == "gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let reader = BufReader::new(reader);
    let mut gff_reader = noodles::gff::io::Reader::new(reader);

    let mut tx_exons: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
    let mut tx_cds: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
    let mut tx_to_gene: HashMap<String, String> = HashMap::new();
    let mut tx_info: HashMap<String, (String, Strand)> = HashMap::new();
    let mut gene_symbols: HashMap<String, String> = HashMap::new();

    let mut raw_id_to_gene_id: HashMap<String, String> = HashMap::new();
    let mut raw_id_to_tx_id: HashMap<String, String> = HashMap::new();

    // Phase 1: Keep raw parent IDs during parsing
    let mut tx_exons_raw: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
    let mut tx_cds_raw: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
    let mut tx_to_gene_raw: HashMap<String, String> = HashMap::new();

    for result in gff_reader.record_bufs() {
        let record = result?;

        let contig = record.reference_sequence_name().to_string();
        let feature = record.ty().to_string();
        let strand = record.strand();

        let start = usize::from(record.start()) as i32 - 1;
        let end = usize::from(record.end()) as i32;

        let attrs = record.attributes();
        let get_attr = |key: &str| {
            attrs
                .get(key.as_bytes())
                .and_then(|v| v.as_string())
                .map(|s| s.to_string())
        };

        let raw_id = get_attr(tag::ID);
        let raw_parent = get_attr(tag::PARENT);
        let name = get_attr(tag::NAME).or_else(|| get_attr("gene_name"));

        let resolve_id =
            |id: Option<String>, version: Option<String>, prefixes: &[&str]| -> String {
                match (id, version) {
                    (Some(i), Some(v)) => format!("{i}.{v}"),
                    (Some(i), None) => i,
                    _ => {
                        let mut s = raw_id.clone().unwrap_or_default();
                        for prefix in prefixes {
                            s = s.replace(prefix, "");
                        }
                        s
                    }
                }
            };

        match feature.as_str() {
            f if f.contains("gene") => {
                let resolved_gene_id = resolve_id(
                    get_attr("gene_id"),
                    get_attr("version").or_else(|| get_attr("gene_version")),
                    &["gene:"],
                );

                if let Some(rid) = &raw_id {
                    raw_id_to_gene_id.insert(rid.clone(), resolved_gene_id.clone());
                }

                if !resolved_gene_id.is_empty() {
                    let gene_id = GeneId::Gene(resolved_gene_id.clone());
                    loader.gene_id_to_gene.insert(
                        gene_id.clone(),
                        Gene {
                            hgnc: Some(gene_id.to_string()),
                            gene_symbol: name.clone(),
                            aliases: None,
                            biotype: None,
                            description: None,
                            map_location: None,
                            summary: None,
                            url: String::new(),
                        },
                    );
                    if let Some(n) = name {
                        gene_symbols.insert(resolved_gene_id, n);
                    }
                }
            }
            f if f.contains("transcript") || f.contains("mRNA") || f.ends_with("RNA") => {
                let resolved_tx_id = resolve_id(
                    get_attr("transcript_id"),
                    get_attr("version").or_else(|| get_attr("transcript_version")),
                    &["transcript:", "rna:", "rna-"],
                );

                if let Some(rid) = &raw_id {
                    raw_id_to_tx_id.insert(rid.clone(), resolved_tx_id.clone());
                }

                if !resolved_tx_id.is_empty() {
                    if let Some(p) = raw_parent {
                        let first_parent = p.split(',').next().unwrap().to_string();
                        // Store raw parent ID for later resolution
                        tx_to_gene_raw.insert(resolved_tx_id.clone(), first_parent);
                    }
                    tx_info.insert(resolved_tx_id, (contig, strand));
                }
            }
            "exon" | "CDS" => {
                if let Some(p) = raw_parent {
                    let target_map = if feature == "exon" {
                        &mut tx_exons_raw
                    } else {
                        &mut tx_cds_raw
                    };

                    for parent_id in p.split(',') {
                        // Store with raw parent ID
                        target_map
                            .entry(parent_id.to_string())
                            .or_default()
                            .push((start, end));
                    }
                }
            }
            _ => {}
        }
    }

    // Phase 2: After parsing, resolve all raw parent IDs
    for (raw_parent, gene_id) in tx_to_gene_raw {
        if let Some(resolved_gene) = raw_id_to_gene_id.get(&gene_id) {
            tx_to_gene.insert(raw_parent, resolved_gene.clone());
        } else {
            // Fallback: use the raw parent as-is
            tx_to_gene.insert(raw_parent, gene_id);
        }
    }

    for (raw_parent, exons_list) in tx_exons_raw {
        if let Some(resolved_tx) = raw_id_to_tx_id.get(&raw_parent) {
            tx_exons.entry(resolved_tx.clone())
                .or_default()
                .extend(exons_list);
        } else {
            // Fallback: use the raw parent as-is
            tx_exons.entry(raw_parent)
                .or_default()
                .extend(exons_list);
        }
    }

    for (raw_parent, cds_list) in tx_cds_raw {
        if let Some(resolved_tx) = raw_id_to_tx_id.get(&raw_parent) {
            tx_cds.entry(resolved_tx.clone())
                .or_default()
                .extend(cds_list);
        } else {
            // Fallback: use the raw parent as-is
            tx_cds.entry(raw_parent)
                .or_default()
                .extend(cds_list);
        }
    }

    // Finalize transcripts by resolving genomic-to-transcript coordinates
    for (tx_id, (contig, gff_strand)) in tx_info {
        let mut exons = tx_exons.remove(&tx_id).unwrap_or_default();
        let cds_fragments = tx_cds.remove(&tx_id).unwrap_or_default();

        if exons.is_empty() {
            continue;
        }

        // Sort exons by genomic position
        exons.sort_by_key(|e| e.0);

        let is_reverse = matches!(gff_strand, Strand::Reverse);
        if is_reverse {
            exons.reverse();
        }

        let tx_strand = if is_reverse {
            cdot_models::Strand::Minus
        } else {
            cdot_models::Strand::Plus
        };

        let mut current_tx_pos = 0;
        let mut tx_cds_start = None;
        let mut tx_cds_end = None;

        let final_exons: Vec<_> = exons
            .into_iter()
            .enumerate()
            .map(|(i, (start, end))| {
                let e_len = end - start;

                // Compute overlaps against actual per-exon CDS fragments
                for cds_fragment in &cds_fragments {
                    let (cds_start, cds_end) = *cds_fragment;
                    let overlap_start = cds_start.max(start);
                    let overlap_end = cds_end.min(end);

                    if overlap_start < overlap_end {
                        // Calculate offset within this exon based on strand
                        let (offset_start, offset_end) = if !is_reverse {
                            (overlap_start - start, overlap_end - start)
                        } else {
                            (end - overlap_end, end - overlap_start)
                        };

                        if tx_cds_start.is_none() {
                            tx_cds_start = Some((current_tx_pos + offset_start) as u32);
                        }
                        tx_cds_end = Some((current_tx_pos + offset_end) as u32);
                    }
                }

                let exon_record = cdot_models::Exon {
                    alt_start_i: start,
                    alt_end_i: end,
                    ord: i as i32,
                    alt_cds_start_i: current_tx_pos + 1,
                    alt_cds_end_i: current_tx_pos + e_len,
                    cigar: format!("{}M", e_len),
                };

                current_tx_pos += e_len;
                exon_record
            })
            .collect();

        let cds_start_genomic = cds_fragments.iter().map(|c| c.0).min();
        let cds_end_genomic = cds_fragments.iter().map(|c| c.1).max();

        let alignment = GenomeAlignment {
            contig,
            strand: tx_strand,
            cds_start: cds_start_genomic,
            cds_end: cds_end_genomic,
            exons: final_exons,
            tag: None,
            note: None,
        };

        let gene_ref = tx_to_gene
            .get(&tx_id)
            .cloned()
            .unwrap_or_else(|| tx_id.clone());
        let gene_name = gene_symbols.get(&gene_ref).cloned();
        let fake_gene_id = GeneId::Gene(gene_ref.clone());

        let transcript = Transcript {
            id: tx_id.clone(),
            hgnc: Some(fake_gene_id.to_string()),
            gene_name,
            gene_version: "".to_string(),
            biotype: None,
            protein: tx_cds_start.map(|_| "unspecified_protein".to_string()),
            start_codon: tx_cds_start.map(i32::try_from).transpose()?,
            stop_codon: tx_cds_end.map(i32::try_from).transpose()?,
            partial: None,
            genome_builds: IndexMap::from([(loader.genome_release.clone(), alignment)]),
        };

        let t_id = TranscriptId::try_new(tx_id)?;
        loader
            .transcript_id_to_transcript
            .insert(t_id.clone(), transcript);
        loader
            .gene_id_to_transcript_ids
            .entry(fake_gene_id.clone())
            .or_default()
            .push(t_id);

        // Ensure the gene entry exists even if no explicit 'gene' feature was in GFF
        loader
            .gene_id_to_gene
            .entry(fake_gene_id)
            .or_insert_with(|| Gene {
                hgnc: Some(gene_ref),
                gene_symbol: None,
                aliases: None,
                biotype: None,
                description: None,
                map_location: None,
                summary: None,
                url: "".into(),
            });
    }
    Ok(())
}