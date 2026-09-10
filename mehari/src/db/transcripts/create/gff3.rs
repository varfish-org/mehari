use crate::db::transcripts::create::cdot_models;
use crate::db::transcripts::create::models::{GeneId, TranscriptId, TranscriptLoader};
use anyhow::Error;
use hgvs::data::cdot::json::models::{Gene, GenomeAlignment, Transcript};
use indexmap::IndexMap;
use noodles::gff::feature::record::{Phase, Strand};
use noodles::gff::feature::record_buf::attributes::field::tag;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Load and extract from standard generic GFF3 using noodles::gff.
pub fn load_gff3(loader: &mut TranscriptLoader, path: impl AsRef<Path>) -> Result<(), Error> {
    let file = File::open(path.as_ref())?;
    let reader: Box<dyn std::io::Read> = if path.as_ref().extension().is_some_and(|e| e == "gz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let reader = BufReader::new(reader);
    let mut gff_reader = noodles::gff::io::Reader::new(reader);

    let mut tx_exons: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
    // Keep each CDS fragment's GFF3 phase (column 8) next to its (start, end). Below, the
    // first fragment in transcript direction is shifted by its phase so the CDS starts on a
    // codon boundary; after that the phase is dropped again.
    let mut tx_cds: HashMap<String, Vec<(i32, i32, u8)>> = HashMap::new();
    let mut tx_to_gene: HashMap<String, String> = HashMap::new();
    let mut tx_info: HashMap<String, (String, Strand)> = HashMap::new();
    let mut gene_symbols: HashMap<String, String> = HashMap::new();
    let mut tx_to_gene_name: HashMap<String, String> = HashMap::new();

    let mut raw_id_to_gene_id: HashMap<String, String> = HashMap::new();
    let mut raw_id_to_tx_id: HashMap<String, String> = HashMap::new();

    // Phase 1: Keep raw parent IDs during parsing
    let mut tx_exons_raw: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
    let mut tx_cds_raw: HashMap<String, Vec<(i32, i32, u8)>> = HashMap::new();
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
                    // Capture transcript-level gene_name
                    if let Some(gene_name) = name.clone() {
                        tx_to_gene_name.insert(resolved_tx_id.clone(), gene_name);
                    }
                    tx_info.insert(resolved_tx_id, (contig, strand));
                }
            }
            "exon" => {
                if let Some(p) = raw_parent {
                    for parent_id in p.split(',') {
                        // Store with raw parent ID
                        tx_exons_raw
                            .entry(parent_id.to_string())
                            .or_default()
                            .push((start, end));
                    }
                }
            }
            "CDS" => {
                if let Some(p) = raw_parent {
                    // Phase (GFF3 column 8) is required for CDS records; default to 0
                    // (in-frame) for malformed input rather than failing the whole file.
                    let phase = match record.phase() {
                        Some(Phase::Zero) | None => 0u8,
                        Some(Phase::One) => 1,
                        Some(Phase::Two) => 2,
                    };

                    for parent_id in p.split(',') {
                        // Store with raw parent ID
                        tx_cds_raw
                            .entry(parent_id.to_string())
                            .or_default()
                            .push((start, end, phase));
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
            tx_exons
                .entry(resolved_tx.clone())
                .or_default()
                .extend(exons_list);
        } else {
            // Fallback: use the raw parent as-is
            tx_exons.entry(raw_parent).or_default().extend(exons_list);
        }
    }

    for (raw_parent, cds_list) in tx_cds_raw {
        if let Some(resolved_tx) = raw_id_to_tx_id.get(&raw_parent) {
            tx_cds
                .entry(resolved_tx.clone())
                .or_default()
                .extend(cds_list);
        } else {
            // Fallback: use the raw parent as-is
            tx_cds.entry(raw_parent).or_default().extend(cds_list);
        }
    }

    // Finalize transcripts by resolving genomic-to-transcript coordinates
    for (tx_id, (contig, gff_strand)) in tx_info {
        let mut exons = tx_exons.remove(&tx_id).unwrap_or_default();
        let mut cds_fragments = tx_cds.remove(&tx_id).unwrap_or_default();

        if exons.is_empty() {
            continue;
        }

        // Sort exons by genomic position
        exons.sort_by_key(|e| e.0);

        let is_reverse = matches!(gff_strand, Strand::Reverse);
        if is_reverse {
            exons.reverse();
        }

        // Honor the GFF3 CDS `phase`: it counts how many bases of the *previous* codon
        // are already consumed at the first base of a CDS fragment. For a 5'-incomplete
        // transcript (e.g. GENCODE's `cds_start_NF` tag) the first CDS fragment in
        // transcript direction has phase 1 or 2, and ignoring it shifts the translated
        // frame by that many bases. Advance that fragment's genomic start (on `+`) or
        // pull back its genomic end (on `-`) by its phase so the CDS -- and thus
        // translation -- start in the right frame. (3'-incomplete CDS lengths, i.e.
        // `cds_end_NF`, are unrelated to phase and stay flagged as InvalidCdsLength.)
        let first_fragment = if is_reverse {
            cds_fragments.iter_mut().max_by_key(|(_, end, _)| *end)
        } else {
            cds_fragments.iter_mut().min_by_key(|(start, _, _)| *start)
        };
        if let Some((start, end, phase)) = first_fragment {
            let phase = i32::from(*phase);
            if is_reverse {
                *end -= phase;
            } else {
                *start += phase;
            }
        }
        let cds_fragments: Vec<(i32, i32)> = cds_fragments
            .into_iter()
            .map(|(start, end, _phase)| (start, end))
            .collect();

        let tx_strand = if is_reverse {
            cdot_models::Strand::Minus
        } else {
            cdot_models::Strand::Plus
        };

        let mut current_tx_pos = 0;
        let mut tx_cds_start = None;
        let mut tx_cds_end = None;

        let mut final_exons: Vec<_> = exons
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

        // Store exons in ascending genomic order (as cdot does), regardless of strand;
        // `ord` above already reflects transcript direction and decreases along this list for minus-strand transcripts.
        final_exons.sort_by_key(|e| e.alt_start_i);

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
        let gene_name = gene_symbols
            .get(&gene_ref)
            .cloned()
            .or_else(|| tx_to_gene_name.get(&tx_id).cloned());
        let fake_gene_id = GeneId::Gene(gene_ref.clone());

        let transcript = Transcript {
            id: tx_id.clone(),
            hgnc: Some(fake_gene_id.to_string()),
            gene_name: gene_name.clone(),
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
                gene_symbol: gene_name,
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

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    /// A plus-strand transcript whose first (only) CDS fragment has phase 1, and a
    /// minus-strand transcript whose first (only) CDS fragment has phase 2 -- as GENCODE
    /// emits for 5'-incomplete transcripts tagged `cds_start_NF`.
    const GFF3: &str = "\
##gff-version 3
chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=gene:G1P;Name=G1P
chr1\ttest\ttranscript\t1\t1000\t.\t+\t.\tID=transcript:T1P;Parent=gene:G1P
chr1\ttest\texon\t1\t1000\t.\t+\t.\tID=exon:T1P.1;Parent=transcript:T1P
chr1\ttest\tCDS\t101\t400\t.\t+\t1\tID=cds:T1P.1;Parent=transcript:T1P
chr1\ttest\tgene\t2001\t3000\t.\t-\t.\tID=gene:G2M;Name=G2M
chr1\ttest\ttranscript\t2001\t3000\t.\t-\t.\tID=transcript:T2M;Parent=gene:G2M
chr1\ttest\texon\t2001\t3000\t.\t-\t.\tID=exon:T2M.1;Parent=transcript:T2M
chr1\ttest\tCDS\t2301\t2600\t.\t-\t2\tID=cds:T2M.1;Parent=transcript:T2M
";

    /// One plus-strand and one minus-strand transcript, each with three exons.
    const GFF3_THREE_EXONS: &str = "\
##gff-version 3
chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=gene:G1;Name=G1
chr1\ttest\ttranscript\t1\t1000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1
chr1\ttest\texon\t1\t100\t.\t+\t.\tID=exon:T1.1;Parent=transcript:T1
chr1\ttest\texon\t301\t400\t.\t+\t.\tID=exon:T1.2;Parent=transcript:T1
chr1\ttest\texon\t601\t700\t.\t+\t.\tID=exon:T1.3;Parent=transcript:T1
chr1\ttest\tgene\t2001\t3000\t.\t-\t.\tID=gene:G2;Name=G2
chr1\ttest\ttranscript\t2001\t3000\t.\t-\t.\tID=transcript:T2;Parent=gene:G2
chr1\ttest\texon\t2001\t2100\t.\t-\t.\tID=exon:T2.1;Parent=transcript:T2
chr1\ttest\texon\t2301\t2400\t.\t-\t.\tID=exon:T2.2;Parent=transcript:T2
chr1\ttest\texon\t2601\t2700\t.\t-\t.\tID=exon:T2.3;Parent=transcript:T2
";

    fn load(gff3: &str) -> TranscriptLoader {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(gff3.as_bytes()).unwrap();
        let mut loader = TranscriptLoader::new("GRCh38".to_string(), false);
        load_gff3(&mut loader, file.path()).unwrap();
        loader
    }

    #[test]
    fn cds_start_is_advanced_by_phase_on_plus_strand() {
        let loader = load(GFF3);

        let tx = loader
            .transcript_id_to_transcript
            .get(&TranscriptId::try_new("T1P").unwrap())
            .unwrap();
        let alignment = tx.genome_builds.get("GRCh38").unwrap();

        // Phase 1 on the (0-based) fragment (100, 400) moves the genomic CDS start
        // one base to the right; the CDS end is untouched.
        assert_eq!(alignment.cds_start, Some(101));
        assert_eq!(alignment.cds_end, Some(400));
    }

    #[test]
    fn cds_end_is_pulled_back_by_phase_on_minus_strand() {
        let loader = load(GFF3);

        let tx = loader
            .transcript_id_to_transcript
            .get(&TranscriptId::try_new("T2M").unwrap())
            .unwrap();
        let alignment = tx.genome_builds.get("GRCh38").unwrap();

        // Phase 2 on the (0-based) fragment (2300, 2600) moves the genomic CDS end
        // two bases to the left (the transcript-direction CDS start, since this
        // transcript is on the `-` strand); the CDS start is untouched.
        assert_eq!(alignment.cds_start, Some(2300));
        assert_eq!(alignment.cds_end, Some(2598));
    }

    /// `bgzip` output is a multi-member gzip stream (one gzip member per block). A plain
    /// `GzDecoder` only reads the first member, so make sure `load_gff3` reads all of them.
    #[test]
    fn load_gff3_reads_all_members_of_multi_member_gzip() -> Result<(), anyhow::Error> {
        let half_a = "##gff-version 3\n\
            chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=gene1;gene_id=GENE1;Name=GENE1\n\
            chr1\ttest\tmRNA\t1\t1000\t.\t+\t.\tID=tx1;Parent=gene1;transcript_id=TX1\n\
            chr1\ttest\texon\t1\t500\t.\t+\t.\tID=exon1;Parent=tx1\n\
            chr1\ttest\texon\t600\t1000\t.\t+\t.\tID=exon2;Parent=tx1\n";
        let half_b = "chr1\ttest\tgene\t2000\t3000\t.\t+\t.\tID=gene2;gene_id=GENE2;Name=GENE2\n\
            chr1\ttest\tmRNA\t2000\t3000\t.\t+\t.\tID=tx2;Parent=gene2;transcript_id=TX2\n\
            chr1\ttest\texon\t2000\t2500\t.\t+\t.\tID=exon3;Parent=tx2\n\
            chr1\ttest\texon\t2600\t3000\t.\t+\t.\tID=exon4;Parent=tx2\n";

        // Two independently gzip-compressed halves concatenated, like `bgzip` produces.
        let mut gzipped = Vec::new();
        for half in [half_a, half_b] {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(half.as_bytes())?;
            gzipped.extend(encoder.finish()?);
        }

        let dir = tempfile::tempdir()?;
        let path = dir.path().join("annotation.gff3.gz");
        std::fs::write(&path, &gzipped)?;

        let mut loader = TranscriptLoader::new("GRCh38".to_string(), false);
        load_gff3(&mut loader, &path)?;

        let mut ids = loader
            .transcript_id_to_transcript
            .keys()
            .map(|id| id.to_string())
            .collect::<Vec<_>>();
        ids.sort();
        assert_eq!(ids, vec!["TX1".to_string(), "TX2".to_string()]);

        Ok(())
    }

    #[test]
    fn exons_are_stored_in_ascending_genomic_order() {
        let loader = load(GFF3_THREE_EXONS);

        let plus_tx = loader
            .transcript_id_to_transcript
            .get(&TranscriptId::try_new("T1").unwrap())
            .unwrap();
        let plus_exons = &plus_tx.genome_builds.get("GRCh38").unwrap().exons;
        assert_eq!(
            plus_exons.iter().map(|e| e.alt_start_i).collect::<Vec<_>>(),
            vec![0, 300, 600],
            "plus-strand exons must be in ascending genomic order"
        );
        assert_eq!(
            plus_exons.iter().map(|e| e.ord).collect::<Vec<_>>(),
            vec![0, 1, 2],
            "plus-strand ord must increase along the (ascending) exon list"
        );

        let minus_tx = loader
            .transcript_id_to_transcript
            .get(&TranscriptId::try_new("T2").unwrap())
            .unwrap();
        let minus_exons = &minus_tx.genome_builds.get("GRCh38").unwrap().exons;
        assert_eq!(
            minus_exons
                .iter()
                .map(|e| e.alt_start_i)
                .collect::<Vec<_>>(),
            vec![2000, 2300, 2600],
            "minus-strand exons must also be stored in ascending genomic order"
        );
        assert_eq!(
            minus_exons.iter().map(|e| e.ord).collect::<Vec<_>>(),
            vec![2, 1, 0],
            "minus-strand ord must decrease along the (ascending) exon list"
        );
    }
}
