use crate::common::progress::{Progress, open_with_progress};
use crate::db::transcripts::create::cdot_models;
use crate::db::transcripts::create::filter::MITOCHONDRIAL_ACCESSIONS;
use crate::db::transcripts::create::models::{GeneId, TranscriptId, TranscriptLoader};
use anyhow::{Context, Error, anyhow};
use hgvs::data::cdot::json::models::{BioType, Gene, GenomeAlignment, Tag, Transcript};
use indexmap::IndexMap;
use noodles::gff::feature::record::{Phase, Strand};
use noodles::gff::feature::record_buf::attributes::field::tag;
use serde::Deserialize;
use serde::de::IntoDeserializer;
use serde::de::value::StrDeserializer;
use std::collections::{HashMap, HashSet};
use std::io::BufReader;
use std::path::Path;

/// What a transcript row says about its transcript.
struct TranscriptRow {
    /// `ID` of the row, which its exon and CDS rows name as `Parent`.
    raw_id: Option<String>,
    /// `Parent` of the row, usually the `ID` of its gene row.
    raw_parent: Option<String>,
    gene_name: Option<String>,
    contig: String,
    strand: Strand,
    /// GFF3 type of the row, e.g. `mRNA` or `lnc_RNA`.
    ty: String,
    /// Biotype attribute of the row: `biotype` in Ensembl, `transcript_type` in GENCODE.
    biotype: Option<String>,
    tags: Vec<Tag>,
    partial: bool,
}

/// What the exon and CDS rows of a transcript say about it.
#[derive(Default)]
struct TranscriptChildren {
    exons: Vec<(i32, i32)>,
    // Keep each CDS fragment's GFF3 phase (column 8) next to its (start, end). Below, the
    // first fragment in transcript direction is shifted by its phase so the CDS starts on a
    // codon boundary; after that the phase is dropped again.
    cds: Vec<(i32, i32, u8)>,
    protein: Option<String>,
    note: Option<String>,
}

/// An exon and its alignment to the genome, as cdot stores it.
struct AlignedExon {
    /// 0-based, half-open genomic coordinates.
    start: i32,
    end: i32,
    /// 1-based, inclusive transcript coordinates.
    tx_start: i32,
    tx_end: i32,
    /// GFF3 `Gap` operations in transcript direction, empty for a gapless exon.
    gap: Vec<(u8, i32)>,
}

impl AlignedExon {
    /// The CIGAR string in the convention of the hgvs mapper, as `cdot_models::gap_to_cigar`
    /// writes it: GFF3 `I` (transcript-only bases) becomes `D`, and GFF3 `D` becomes `I`.
    fn cigar(&self) -> String {
        if self.gap.is_empty() {
            return format!("{}M", self.end - self.start);
        }
        self.gap
            .iter()
            .map(|&(op, len)| {
                let op = match op {
                    b'I' => 'D',
                    b'D' => 'I',
                    _ => '=',
                };
                format!("{len}{op}")
            })
            .collect()
    }

    /// The 0-based transcript position of the genomic base `pos`, if this exon aligns it.
    fn tx_position(&self, pos: i32, is_reverse: bool) -> Option<i32> {
        if pos < self.start || pos >= self.end {
            return None;
        }
        // Genomic bases between the 5' end of the exon and `pos`.
        let mut offset = if is_reverse {
            self.end - 1 - pos
        } else {
            pos - self.start
        };
        let mut tx_pos = self.tx_start - 1;
        let gapless = [(b'M', self.end - self.start)];
        let ops = if self.gap.is_empty() {
            &gapless[..]
        } else {
            &self.gap[..]
        };
        for &(op, len) in ops {
            match op {
                // Bases in the transcript only.
                b'I' => tx_pos += len,
                // Bases in the genome only.
                b'D' if offset < len => return None,
                b'D' => offset -= len,
                _ if offset < len => return Some(tx_pos + offset),
                _ => {
                    offset -= len;
                    tx_pos += len;
                }
            }
        }
        None
    }
}

/// Parse a GFF3 `Target` value, e.g. `NM_001304717.2 1 1110 +`, into the accession and its
/// 1-based start and end.
fn parse_target(target: &str) -> Result<(String, i32, i32), Error> {
    let mut fields = target.split_whitespace();
    match (fields.next(), fields.next(), fields.next()) {
        (Some(accession), Some(start), Some(end)) => Ok((
            accession.to_string(),
            start
                .parse()
                .with_context(|| format!("invalid Target {target:?}"))?,
            end.parse()
                .with_context(|| format!("invalid Target {target:?}"))?,
        )),
        _ => Err(anyhow!("invalid Target {target:?}")),
    }
}

/// Parse a GFF3 `Gap` value, e.g. `M185 I3 M250`.
fn parse_gap(gap: &str) -> Result<Vec<(u8, i32)>, Error> {
    gap.split_whitespace()
        .map(|op| -> Result<(u8, i32), Error> {
            match op.as_bytes() {
                [code @ (b'M' | b'I' | b'D'), ..] => Ok((
                    *code,
                    op[1..]
                        .parse()
                        .with_context(|| format!("invalid operation {op:?} in Gap {gap:?}"))?,
                )),
                _ => Err(anyhow!("invalid operation {op:?} in Gap {gap:?}")),
            }
        })
        .collect()
}

/// Parse a GFF3 type or biotype, e.g. `lnc_RNA` or `protein_coding`.
fn parse_biotype(value: &str) -> Option<BioType> {
    // Ensembl and GENCODE write `nonsense_mediated_decay`. `BioType` names it by the SO term
    // `NMD_transcript_variant`.
    if value == "nonsense_mediated_decay" {
        return Some(BioType::NmdTranscriptVariant);
    }
    let deserializer: StrDeserializer<'_, serde::de::value::Error> = value.into_deserializer();
    BioType::deserialize(deserializer)
        .inspect_err(|_| tracing::debug!("skipping unknown biotype {value}"))
        .ok()
}

/// Load and extract from standard generic GFF3 using noodles::gff.
///
/// Reads the same transcript fields as cdot does, plus the biotype attribute of transcripts.
/// For RefSeq, these fields include the `cDNA_match` rows, which align some transcripts to the
/// genome, including gaps. Genes keep the gene ID of the annotation, whereas `load_cdot` names
/// them by HGNC ID.
pub fn load_gff3(
    loader: &mut TranscriptLoader,
    path: impl AsRef<Path>,
    progress: &dyn Progress,
) -> Result<(), Error> {
    let file = open_with_progress(path.as_ref(), progress)?;
    let bar = file.progress.clone();
    let reader: Box<dyn std::io::Read> = if path.as_ref().extension().is_some_and(|e| e == "gz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let reader = BufReader::new(reader);
    let mut gff_reader = noodles::gff::io::Reader::new(reader);

    let mut tx_rows: HashMap<String, TranscriptRow> = HashMap::new();
    // The `ID` of every transcript row, including those of skipped locations.
    let mut tx_raw_ids: HashSet<String> = HashSet::new();
    // Keyed by the `Parent` of the exon and CDS rows.
    let mut raw_children: HashMap<String, TranscriptChildren> = HashMap::new();
    // Keyed by the transcript accession in `Target`, with the contig of each row.
    let mut cdna_matches: HashMap<String, Vec<(String, AlignedExon)>> = HashMap::new();
    let mut raw_id_to_gene_id: HashMap<String, GeneId> = HashMap::new();
    let mut gene_symbols: HashMap<GeneId, String> = HashMap::new();

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
        // All values of an attribute, e.g. `Dbxref=GeneID:4535,HGNC:HGNC:7455`.
        let get_values = |key: &str| -> Vec<String> {
            attrs
                .get(key.as_bytes())
                .map(|v| v.iter().map(|s| s.to_string()).collect())
                .unwrap_or_default()
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
            // As cdot does, read only rows without `Parent` as genes. A RefSeq gene segment
            // below its gene names the same NCBI Gene ID and would replace the gene.
            f if f.contains("gene") && raw_parent.is_none() => {
                // Genes keep the ID of the annotation. GFF3 has no shared attribute for it:
                // GENCODE and Ensembl write `gene_id`, RefSeq writes the NCBI Gene ID into
                // `Dbxref`, and `ID` is only unique within the file.
                let ncbi_gene_id = get_values("Dbxref")
                    .iter()
                    .find_map(|x| x.strip_prefix("GeneID:"))
                    .map(String::from);
                let resolved_gene_id = resolve_id(
                    get_attr("gene_id").or(ncbi_gene_id),
                    get_attr("version").or_else(|| get_attr("gene_version")),
                    &["gene:"],
                );
                let gene_id = GeneId::Gene(resolved_gene_id.clone());

                if let Some(rid) = &raw_id {
                    raw_id_to_gene_id.insert(rid.clone(), gene_id.clone());
                }

                if !resolved_gene_id.is_empty() {
                    loader.gene_id_to_gene.insert(
                        gene_id.clone(),
                        Gene {
                            hgnc: Some(gene_id.to_string()),
                            gene_symbol: name.clone(),
                            aliases: None,
                            biotype: get_attr("gene_biotype")
                                .or_else(|| get_attr("biotype"))
                                .and_then(|biotype| parse_biotype(&biotype))
                                .map(|biotype| vec![biotype]),
                            description: get_attr("description"),
                            map_location: None,
                            summary: None,
                            url: String::new(),
                        },
                    );
                    if let Some(n) = name {
                        gene_symbols.insert(gene_id, n);
                    }
                }
            }
            f if f.contains("transcript") || f.contains("mRNA") || f.ends_with("RNA") => {
                // RefSeq has no `transcript_id` on mitochondrial mRNAs. Name them as cdot
                // does, plus the version that `load_cdot` adds, e.g. `fake-rna-ND1.0`.
                let mt_mrna_id = (f == "mRNA"
                    && MITOCHONDRIAL_ACCESSIONS.contains(&contig.as_str()))
                .then(|| raw_id.as_ref().map(|id| format!("fake-{id}.0")))
                .flatten();
                let resolved_tx_id = resolve_id(
                    get_attr("transcript_id").or(mt_mrna_id),
                    get_attr("version").or_else(|| get_attr("transcript_version")),
                    &["transcript:", "rna:", "rna-"],
                );

                if let Some(rid) = &raw_id {
                    tx_raw_ids.insert(rid.clone());
                }
                // A transcript can have more than one location, e.g. in the PARs of chrX and
                // chrY. As cdot does, keep the first one and skip the rows of the others.
                if resolved_tx_id.is_empty() || tx_rows.contains_key(&resolved_tx_id) {
                    continue;
                }
                tx_rows.insert(
                    resolved_tx_id,
                    TranscriptRow {
                        raw_id,
                        raw_parent: raw_parent.and_then(|p| p.split(',').next().map(String::from)),
                        gene_name: name,
                        contig,
                        strand,
                        ty: f.to_string(),
                        biotype: get_attr("biotype").or_else(|| get_attr("transcript_type")),
                        tags: get_values("tag")
                            .iter()
                            .map(|t| cdot_models::str_to_tag(t))
                            .collect(),
                        partial: get_attr("partial").is_some(),
                    },
                );
            }
            "exon" => {
                if let Some(p) = raw_parent {
                    for parent_id in p.split(',') {
                        let children = raw_children.entry(parent_id.to_string()).or_default();
                        children.exons.push((start, end));
                        children.note = get_attr("Note").or(children.note.take());
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
                    // RefSeq and GENCODE write `protein_id=NP_001659.1`, Ensembl writes
                    // `protein_id=ENSP00000477624;protein_version=1`.
                    let protein =
                        get_attr("protein_id").map(|id| match get_attr("protein_version") {
                            Some(version) => format!("{id}.{version}"),
                            None => id,
                        });

                    for parent_id in p.split(',') {
                        let children = raw_children.entry(parent_id.to_string()).or_default();
                        children.cds.push((start, end, phase));
                        children.protein = protein.clone().or(children.protein.take());
                        children.note = get_attr("Note").or(children.note.take());
                    }
                }
            }
            "cDNA_match" => {
                let target =
                    get_attr("Target").ok_or_else(|| anyhow!("cDNA_match row without Target"))?;
                let (accession, tx_start, tx_end) = parse_target(&target)?;
                let gap = get_attr("Gap")
                    .map(|gap| parse_gap(&gap))
                    .transpose()?
                    .unwrap_or_default();
                cdna_matches.entry(accession).or_default().push((
                    contig,
                    AlignedExon {
                        start,
                        end,
                        tx_start,
                        tx_end,
                        gap,
                    },
                ));
            }
            _ => {}
        }
    }

    // Finalize transcripts by resolving genomic-to-transcript coordinates
    for (tx_id, row) in tx_rows {
        // As cdot does, skip transcripts whose parent is a transcript, e.g. a mature miRNA
        // below its primary transcript in RefSeq.
        if row
            .raw_parent
            .as_ref()
            .is_some_and(|parent| tx_raw_ids.contains(parent))
        {
            continue;
        }

        let children = row
            .raw_id
            .as_ref()
            .and_then(|id| raw_children.remove(id))
            .or_else(|| raw_children.remove(&tx_id))
            .unwrap_or_default();
        let is_reverse = matches!(row.strand, Strand::Reverse);

        // RefSeq aligns each transcript to the genome in `cDNA_match` rows. As cdot does,
        // prefer them over the exon rows, which cannot show bases that the genome lacks.
        let mut exons: Vec<AlignedExon> = cdna_matches
            .remove(&tx_id)
            .unwrap_or_default()
            .into_iter()
            .filter(|(contig, _)| *contig == row.contig)
            .map(|(_, exon)| exon)
            .collect();
        let from_exon_rows = exons.is_empty();
        if from_exon_rows {
            exons = children
                .exons
                .iter()
                .map(|&(start, end)| AlignedExon {
                    start,
                    end,
                    tx_start: 0,
                    tx_end: 0,
                    gap: Vec::new(),
                })
                .collect();
        }

        if exons.is_empty() {
            continue;
        }

        // Sort exons into transcript order.
        exons.sort_by_key(|e| e.start);
        if is_reverse {
            exons.reverse();
        }
        if from_exon_rows {
            let mut tx_end = 0;
            for exon in &mut exons {
                exon.tx_start = tx_end + 1;
                tx_end += exon.end - exon.start;
                exon.tx_end = tx_end;
            }
        }

        // Honor the GFF3 CDS `phase`: it counts how many bases of the *previous* codon
        // are already consumed at the first base of a CDS fragment. For a 5'-incomplete
        // transcript (e.g. GENCODE's `cds_start_NF` tag) the first CDS fragment in
        // transcript direction has phase 1 or 2, and ignoring it shifts the translated
        // frame by that many bases. Advance that fragment's genomic start (on `+`) or
        // pull back its genomic end (on `-`) by its phase so the CDS -- and thus
        // translation -- start in the right frame. (3'-incomplete CDS lengths, i.e.
        // `cds_end_NF`, are unrelated to phase and stay flagged as InvalidCdsLength.)
        let mut cds_fragments = children.cds;
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
            // A fragment no longer than its phase holds no base of the CDS.
            .filter(|(start, end)| start < end)
            .collect();

        let tx_strand = if is_reverse {
            cdot_models::Strand::Minus
        } else {
            cdot_models::Strand::Plus
        };

        let cds_start_genomic = cds_fragments.iter().map(|c| c.0).min();
        let cds_end_genomic = cds_fragments.iter().map(|c| c.1).max();

        // Map the first and last CDS base in transcript direction to transcript positions.
        let tx_position = |pos| {
            exons
                .iter()
                .find_map(|exon| exon.tx_position(pos, is_reverse))
        };
        let (tx_cds_start, tx_cds_end) = match (cds_start_genomic, cds_end_genomic) {
            (Some(start), Some(end)) => {
                let (first, last) = if is_reverse {
                    (end - 1, start)
                } else {
                    (start, end - 1)
                };
                (tx_position(first), tx_position(last).map(|pos| pos + 1))
            }
            _ => (None, None),
        };
        if cds_start_genomic.is_some() && (tx_cds_start.is_none() || tx_cds_end.is_none()) {
            tracing::warn!("CDS of {tx_id} is not within its aligned exons");
        }

        let mut final_exons: Vec<_> = exons
            .iter()
            .enumerate()
            .map(|(i, exon)| cdot_models::Exon {
                alt_start_i: exon.start,
                alt_end_i: exon.end,
                ord: i as i32,
                alt_cds_start_i: exon.tx_start,
                alt_cds_end_i: exon.tx_end,
                cigar: exon.cigar(),
            })
            .collect();

        // Store exons in ascending genomic order (as cdot does), regardless of strand;
        // `ord` above already reflects transcript direction and decreases along this list for minus-strand transcripts.
        final_exons.sort_by_key(|e| e.alt_start_i);

        let alignment = GenomeAlignment {
            contig: row.contig,
            strand: tx_strand,
            cds_start: cds_start_genomic,
            cds_end: cds_end_genomic,
            exons: final_exons,
            tag: (!row.tags.is_empty()).then_some(row.tags),
            note: children.note,
        };

        // The biotypes are the GFF3 type of the transcript row, its biotype attribute, and
        // `mRNA` or `ncRNA`. cdot reads the type and `mRNA` or `ncRNA` only. Ensembl and
        // GENCODE write NMD transcripts as `mRNA` and `transcript` rows, so only the
        // attribute marks them.
        let mut biotypes = Vec::new();
        if row.ty != "transcript" {
            biotypes.extend(parse_biotype(&row.ty));
        }
        let coding = if cds_fragments.is_empty() {
            BioType::NcRna
        } else {
            BioType::MRna
        };
        let attribute = row.biotype.as_deref().and_then(parse_biotype);
        for biotype in attribute.into_iter().chain([coding]) {
            if !biotypes.contains(&biotype) {
                biotypes.push(biotype);
            }
        }

        let gene_id = match &row.raw_parent {
            Some(raw) => raw_id_to_gene_id
                .get(raw)
                .cloned()
                .unwrap_or_else(|| GeneId::Gene(raw.clone())),
            None => GeneId::Gene(tx_id.clone()),
        };
        let gene_name = gene_symbols.get(&gene_id).cloned().or(row.gene_name);

        let transcript = Transcript {
            id: tx_id.clone(),
            hgnc: Some(gene_id.to_string()),
            gene_name: gene_name.clone(),
            gene_version: "".to_string(),
            biotype: Some(biotypes.clone()),
            protein: tx_cds_start.map(|_| {
                children
                    .protein
                    .unwrap_or_else(|| "unspecified_protein".to_string())
            }),
            start_codon: tx_cds_start,
            stop_codon: tx_cds_end,
            partial: row.partial.then_some(1),
            genome_builds: IndexMap::from([(loader.genome_release.clone(), alignment)]),
        };

        let t_id = TranscriptId::try_new(tx_id)?;
        loader
            .transcript_id_to_transcript
            .insert(t_id.clone(), transcript);
        loader
            .gene_id_to_transcript_ids
            .entry(gene_id.clone())
            .or_default()
            .push(t_id);

        // Ensure the gene entry exists even if no explicit 'gene' feature was in GFF
        let gene = loader
            .gene_id_to_gene
            .entry(gene_id.clone())
            .or_insert_with(|| Gene {
                hgnc: Some(gene_id.to_string()),
                gene_symbol: gene_name,
                aliases: None,
                biotype: None,
                description: None,
                map_location: None,
                summary: None,
                url: "".into(),
            });
        // As cdot does, a gene also has the biotypes of its transcripts.
        let gene_biotypes = gene.biotype.get_or_insert_with(Vec::new);
        for biotype in biotypes {
            if !gene_biotypes.contains(&biotype) {
                gene_biotypes.push(biotype);
            }
        }
    }

    bar.finish();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::progress::NoProgress;
    use crate::db::transcripts::create::filter::filter_transcripts;
    use crate::db::transcripts::create::models::{Identifier, Reason, TranscriptExt};
    use anyhow::Context;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use hgvs::data::cdot::json::models::{BioType, Tag};
    use hgvs::data::interface::TxExonsRecord;
    use hgvs::mapper::alignment::build_tx_cigar;
    use hgvs::mapper::cigar::CigarMapper;
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

    /// A plus-strand and a minus-strand transcript with two exons each. The CDS runs to the
    /// transcript end and is 190 bases long, as for 3'-incomplete transcripts (GENCODE tag
    /// `cds_end_NF`), so `fix_cds` pads it by 2 bases.
    const GFF3_CDS_END_NF: &str = "\
##gff-version 3
chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=gene:G3P;Name=G3P
chr1\ttest\ttranscript\t1\t400\t.\t+\t.\tID=transcript:T3P;Parent=gene:G3P
chr1\ttest\texon\t1\t100\t.\t+\t.\tID=exon:T3P.1;Parent=transcript:T3P
chr1\ttest\texon\t301\t400\t.\t+\t.\tID=exon:T3P.2;Parent=transcript:T3P
chr1\ttest\tCDS\t11\t100\t.\t+\t0\tID=cds:T3P.1;Parent=transcript:T3P
chr1\ttest\tCDS\t301\t400\t.\t+\t0\tID=cds:T3P.2;Parent=transcript:T3P
chr1\ttest\tgene\t2001\t3000\t.\t-\t.\tID=gene:G3M;Name=G3M
chr1\ttest\ttranscript\t2001\t2400\t.\t-\t.\tID=transcript:T3M;Parent=gene:G3M
chr1\ttest\texon\t2001\t2100\t.\t-\t.\tID=exon:T3M.2;Parent=transcript:T3M
chr1\ttest\texon\t2301\t2400\t.\t-\t.\tID=exon:T3M.1;Parent=transcript:T3M
chr1\ttest\tCDS\t2001\t2100\t.\t-\t0\tID=cds:T3M.2;Parent=transcript:T3M
chr1\ttest\tCDS\t2301\t2390\t.\t-\t0\tID=cds:T3M.1;Parent=transcript:T3M
";

    /// RefSeq-style transcripts. Their `cDNA_match` rows hold the alignment to the genome:
    /// `Target` gives the transcript coordinates, `Gap` the indels.
    /// - NM_000001.1 (`+`): `D2` marks 2 bases that exist in the genome only. The exon rows
    ///   turn them into a fake intron. A PAR copy on chrY follows at the end of the file.
    /// - NM_000002.1 (`-`): `I3` marks 3 bases that exist in the transcript only. The exon
    ///   rows cannot show them.
    /// - NM_000003.1 (`+`): transcript base 1 and bases 102 to 131 align nowhere.
    /// - NR_000004.1: a partial lnc_RNA without `cDNA_match` rows.
    /// - NR_000005.1: a miRNA primary transcript. Its mature miRNA row is not a transcript.
    /// - ND1: a mitochondrial mRNA without `transcript_id`.
    /// - IGKC: a gene with a `C_gene_segment` row below it, which names the same Gene ID.
    const GFF3_REFSEQ: &str = "\
##gff-version 3
NC_000001.11\tBestRefSeq\tgene\t1001\t1300\t.\t+\t.\tID=gene-GA;Dbxref=GeneID:11,HGNC:HGNC:1;Name=GA;gene_biotype=protein_coding
NC_000001.11\tBestRefSeq\tmRNA\t1001\t1300\t.\t+\t.\tID=rna-NM_000001.1;Parent=gene-GA;Dbxref=GeneID:11,HGNC:HGNC:1;Name=NM_000001.1;tag=MANE Select;transcript_id=NM_000001.1
NC_000001.11\tBestRefSeq\texon\t1001\t1100\t.\t+\t.\tID=exon-NM_000001.1-1;Parent=rna-NM_000001.1;transcript_id=NM_000001.1
NC_000001.11\tBestRefSeq\texon\t1201\t1250\t.\t+\t.\tID=exon-NM_000001.1-2;Parent=rna-NM_000001.1;transcript_id=NM_000001.1
NC_000001.11\tBestRefSeq\texon\t1253\t1300\t.\t+\t.\tID=exon-NM_000001.1-3;Parent=rna-NM_000001.1;transcript_id=NM_000001.1
NC_000001.11\tBestRefSeq\tCDS\t1011\t1100\t.\t+\t0\tID=cds-NP_000001.1;Parent=rna-NM_000001.1;Note=UGA stop codon recoded as selenocysteine;protein_id=NP_000001.1
NC_000001.11\tBestRefSeq\tCDS\t1201\t1250\t.\t+\t0\tID=cds-NP_000001.1;Parent=rna-NM_000001.1;Note=UGA stop codon recoded as selenocysteine;protein_id=NP_000001.1
NC_000001.11\tBestRefSeq\tCDS\t1253\t1289\t.\t+\t1\tID=cds-NP_000001.1;Parent=rna-NM_000001.1;Note=UGA stop codon recoded as selenocysteine;protein_id=NP_000001.1
NC_000001.11\tRefSeq\tcDNA_match\t1001\t1100\t.\t+\t.\tID=aln-1;Target=NM_000001.1 1 100 +;gap_count=1
NC_000001.11\tRefSeq\tcDNA_match\t1201\t1300\t.\t+\t.\tID=aln-1;Target=NM_000001.1 101 198 +;gap_count=1;Gap=M50 D2 M48
NC_000001.11\tBestRefSeq\tgene\t2001\t2300\t.\t-\t.\tID=gene-GB;Dbxref=GeneID:22,HGNC:HGNC:2;Name=GB;gene_biotype=protein_coding
NC_000001.11\tBestRefSeq\tmRNA\t2001\t2300\t.\t-\t.\tID=rna-NM_000002.1;Parent=gene-GB;Name=NM_000002.1;transcript_id=NM_000002.1
NC_000001.11\tBestRefSeq\texon\t2201\t2300\t.\t-\t.\tID=exon-NM_000002.1-1;Parent=rna-NM_000002.1;transcript_id=NM_000002.1
NC_000001.11\tBestRefSeq\texon\t2001\t2100\t.\t-\t.\tID=exon-NM_000002.1-2;Parent=rna-NM_000002.1;transcript_id=NM_000002.1
NC_000001.11\tBestRefSeq\tCDS\t2201\t2290\t.\t-\t0\tID=cds-NP_000002.1;Parent=rna-NM_000002.1;protein_id=NP_000002.1
NC_000001.11\tBestRefSeq\tCDS\t2050\t2100\t.\t-\t0\tID=cds-NP_000002.1;Parent=rna-NM_000002.1;protein_id=NP_000002.1
NC_000001.11\tRefSeq\tcDNA_match\t2201\t2300\t.\t-\t.\tID=aln-2;Target=NM_000002.1 1 103 +;gap_count=1;Gap=M40 I3 M60
NC_000001.11\tRefSeq\tcDNA_match\t2001\t2100\t.\t-\t.\tID=aln-2;Target=NM_000002.1 104 203 +;gap_count=1
NC_000001.11\tBestRefSeq\tgene\t3001\t3300\t.\t+\t.\tID=gene-GC;Dbxref=GeneID:33;Name=GC;gene_biotype=protein_coding
NC_000001.11\tBestRefSeq\tmRNA\t3001\t3300\t.\t+\t.\tID=rna-NM_000003.1;Parent=gene-GC;Name=NM_000003.1;tag=RefSeq Select;transcript_id=NM_000003.1
NC_000001.11\tBestRefSeq\texon\t3001\t3100\t.\t+\t.\tID=exon-NM_000003.1-1;Parent=rna-NM_000003.1;transcript_id=NM_000003.1
NC_000001.11\tBestRefSeq\texon\t3201\t3300\t.\t+\t.\tID=exon-NM_000003.1-2;Parent=rna-NM_000003.1;transcript_id=NM_000003.1
NC_000001.11\tBestRefSeq\tCDS\t3011\t3100\t.\t+\t0\tID=cds-NP_000003.1;Parent=rna-NM_000003.1;protein_id=NP_000003.1
NC_000001.11\tBestRefSeq\tCDS\t3201\t3251\t.\t+\t0\tID=cds-NP_000003.1;Parent=rna-NM_000003.1;protein_id=NP_000003.1
NC_000001.11\tRefSeq\tcDNA_match\t3001\t3100\t.\t+\t.\tID=aln-3;Target=NM_000003.1 2 101 +;gap_count=0
NC_000001.11\tRefSeq\tcDNA_match\t3201\t3300\t.\t+\t.\tID=aln-3;Target=NM_000003.1 132 231 +;gap_count=0
NC_000001.11\tBestRefSeq\tgene\t4001\t4200\t.\t+\t.\tID=gene-GD;Dbxref=GeneID:44,HGNC:HGNC:4;Name=GD;gene_biotype=lncRNA
NC_000001.11\tBestRefSeq\tlnc_RNA\t4001\t4200\t.\t+\t.\tID=rna-NR_000004.1;Parent=gene-GD;Name=NR_000004.1;partial=true;transcript_id=NR_000004.1
NC_000001.11\tBestRefSeq\texon\t4001\t4100\t.\t+\t.\tID=exon-NR_000004.1-1;Parent=rna-NR_000004.1;transcript_id=NR_000004.1
NC_000001.11\tBestRefSeq\texon\t4151\t4200\t.\t+\t.\tID=exon-NR_000004.1-2;Parent=rna-NR_000004.1;transcript_id=NR_000004.1
NC_000001.11\tBestRefSeq\tgene\t6001\t6100\t.\t+\t.\tID=gene-MIR1;Dbxref=GeneID:55,HGNC:HGNC:5;Name=MIR1;gene_biotype=miRNA
NC_000001.11\tBestRefSeq\tprimary_transcript\t6001\t6100\t.\t+\t.\tID=rna-NR_000005.1;Parent=gene-MIR1;transcript_id=NR_000005.1
NC_000001.11\tBestRefSeq\texon\t6001\t6100\t.\t+\t.\tID=exon-NR_000005.1-1;Parent=rna-NR_000005.1;transcript_id=NR_000005.1
NC_000001.11\tBestRefSeq\tmiRNA\t6021\t6042\t.\t+\t.\tID=rna-MIR1;Parent=rna-NR_000005.1;gene=MIR1
NC_000001.11\tBestRefSeq\texon\t6021\t6042\t.\t+\t.\tID=exon-MIR1-1;Parent=rna-MIR1;gene=MIR1
NC_000001.11\tCurated Genomic\tgene\t7001\t7300\t.\t-\t.\tID=gene-IGKC;Dbxref=GeneID:3514,HGNC:HGNC:5716;Name=IGKC;description=immunoglobulin kappa constant;gene_biotype=C_region
NC_000001.11\tCurated Genomic\tC_gene_segment\t7001\t7300\t.\t-\t.\tID=id-IGKC;Parent=gene-IGKC;Dbxref=GeneID:3514,HGNC:HGNC:5716;gbkey=C_region;gene=IGKC
NC_000001.11\tCurated Genomic\tCDS\t7001\t7300\t.\t-\t0\tID=cds-IGKC;Parent=id-IGKC;Dbxref=GeneID:3514;gene=IGKC
NC_012920.1\tRefSeq\tgene\t3307\t4262\t.\t+\t.\tID=gene-ND1;Dbxref=GeneID:4535,HGNC:HGNC:7455;Name=ND1;gene_biotype=protein_coding
NC_012920.1\tRefSeq\tmRNA\t3307\t4262\t.\t+\t.\tID=rna-ND1;Parent=gene-ND1;Dbxref=GeneID:4535,HGNC:HGNC:7455;gene=ND1
NC_012920.1\tRefSeq\texon\t3307\t4262\t.\t+\t.\tID=exon-ND1-1;Parent=rna-ND1;gene=ND1
NC_012920.1\tRefSeq\tCDS\t3307\t4262\t.\t+\t0\tID=cds-YP_003024026.1;Parent=rna-ND1;protein_id=YP_003024026.1
NC_000024.10\tBestRefSeq\tgene\t5001\t5300\t.\t+\t.\tID=gene-GA-2;Dbxref=GeneID:11,HGNC:HGNC:1;Name=GA;gene_biotype=protein_coding
NC_000024.10\tBestRefSeq\tmRNA\t5001\t5300\t.\t+\t.\tID=rna-NM_000001.1-2;Parent=gene-GA-2;Name=NM_000001.1;tag=MANE Select;transcript_id=NM_000001.1
NC_000024.10\tBestRefSeq\texon\t5001\t5300\t.\t+\t.\tID=exon-NM_000001.1-2-1;Parent=rna-NM_000001.1-2;transcript_id=NM_000001.1
NC_000024.10\tRefSeq\tcDNA_match\t5001\t5300\t.\t+\t.\tID=aln-4;Target=NM_000001.1 1 300 +;gap_count=0
";

    /// A GENCODE transcript with a PAR_Y copy that has its own `ID` but the same
    /// `transcript_id`, and an Ensembl gene whose `description` names its HGNC ID.
    const GFF3_GENCODE_ENSEMBL: &str = "\
##gff-version 3
chrX\tHAVANA\tgene\t1001\t2000\t.\t+\t.\tID=ENSG00000182378.15;gene_id=ENSG00000182378.15;gene_name=PLCXD1
chrX\tHAVANA\ttranscript\t1001\t2000\t.\t+\t.\tID=ENST00000399012.6;Parent=ENSG00000182378.15;gene_id=ENSG00000182378.15;transcript_id=ENST00000399012.6;tag=basic,Ensembl_canonical,MANE_Select
chrX\tHAVANA\texon\t1001\t2000\t.\t+\t.\tID=exon:ENST00000399012.6:1;Parent=ENST00000399012.6;transcript_id=ENST00000399012.6
chrY\tHAVANA\tgene\t5001\t6000\t.\t+\t.\tID=ENSG00000182378.15_PAR_Y;gene_id=ENSG00000182378.15;gene_name=PLCXD1
chrY\tHAVANA\ttranscript\t5001\t6000\t.\t+\t.\tID=ENST00000399012.6_PAR_Y;Parent=ENSG00000182378.15_PAR_Y;gene_id=ENSG00000182378.15;transcript_id=ENST00000399012.6;tag=basic,Ensembl_canonical,MANE_Select
chrY\tHAVANA\texon\t5001\t6000\t.\t+\t.\tID=exon:ENST00000399012.6_PAR_Y:1;Parent=ENST00000399012.6_PAR_Y;transcript_id=ENST00000399012.6
22\tensembl_havana\tgene\t1001\t2000\t.\t+\t.\tID=gene:ENSG00000100001;Name=GE;biotype=protein_coding;description=test gene [Source:HGNC Symbol%3BAcc:HGNC:8907];gene_id=ENSG00000100001;version=3
22\tensembl_havana\tmRNA\t1001\t2000\t.\t+\t.\tID=transcript:ENST00000100001;Parent=gene:ENSG00000100001;biotype=protein_coding;transcript_id=ENST00000100001;version=2
22\tensembl_havana\texon\t1001\t2000\t.\t+\t.\tParent=transcript:ENST00000100001
";

    /// An Ensembl and a GENCODE NMD transcript. Their GFF3 type is `mRNA` and `transcript`,
    /// so only the biotype attribute marks them as NMD transcripts.
    const GFF3_NMD: &str = "\
##gff-version 3
22\tensembl_havana\tgene\t1001\t2000\t.\t+\t.\tID=gene:ENSG00000100002;Name=GN;biotype=protein_coding;gene_id=ENSG00000100002;version=1
22\thavana\tmRNA\t1001\t2000\t.\t+\t.\tID=transcript:ENST00000100002;Parent=gene:ENSG00000100002;biotype=nonsense_mediated_decay;transcript_id=ENST00000100002;version=1
22\thavana\texon\t1001\t2000\t.\t+\t.\tParent=transcript:ENST00000100002
22\thavana\tCDS\t1101\t1400\t.\t+\t0\tParent=transcript:ENST00000100002;protein_id=ENSP00000100002;protein_version=1
chr22\tHAVANA\tgene\t3001\t4000\t.\t+\t.\tID=ENSG00000100003.1;gene_id=ENSG00000100003.1;gene_type=protein_coding;gene_name=GM
chr22\tHAVANA\ttranscript\t3001\t4000\t.\t+\t.\tID=ENST00000100003.1;Parent=ENSG00000100003.1;gene_id=ENSG00000100003.1;transcript_id=ENST00000100003.1;transcript_type=nonsense_mediated_decay
chr22\tHAVANA\texon\t3001\t4000\t.\t+\t.\tID=exon:ENST00000100003.1:1;Parent=ENST00000100003.1
chr22\tHAVANA\tCDS\t3101\t3400\t.\t+\t0\tID=CDS:ENST00000100003.1;Parent=ENST00000100003.1;protein_id=ENSP00000100003.1
";

    fn load(gff3: &str) -> Result<TranscriptLoader, anyhow::Error> {
        let mut file = tempfile::NamedTempFile::new()?;
        file.write_all(gff3.as_bytes())?;
        let mut loader = TranscriptLoader::new("GRCh38".to_string(), false);
        load_gff3(&mut loader, file.path(), &NoProgress)?;
        Ok(loader)
    }

    #[test]
    fn cds_start_is_advanced_by_phase_on_plus_strand() -> Result<(), anyhow::Error> {
        let loader = load(GFF3)?;

        let tx = loader
            .transcript_id_to_transcript
            .get(&TranscriptId::try_new("T1P")?)
            .context("transcript T1P not loaded")?;
        let alignment = tx
            .genome_builds
            .get("GRCh38")
            .context("T1P has no GRCh38 alignment")?;

        // Phase 1 on the (0-based) fragment (100, 400) moves the genomic CDS start
        // one base to the right; the CDS end is untouched.
        assert_eq!(alignment.cds_start, Some(101));
        assert_eq!(alignment.cds_end, Some(400));

        Ok(())
    }

    #[test]
    fn cds_end_is_pulled_back_by_phase_on_minus_strand() -> Result<(), anyhow::Error> {
        let loader = load(GFF3)?;

        let tx = loader
            .transcript_id_to_transcript
            .get(&TranscriptId::try_new("T2M")?)
            .context("transcript T2M not loaded")?;
        let alignment = tx
            .genome_builds
            .get("GRCh38")
            .context("T2M has no GRCh38 alignment")?;

        // Phase 2 on the (0-based) fragment (2300, 2600) moves the genomic CDS end
        // two bases to the left (the transcript-direction CDS start, since this
        // transcript is on the `-` strand); the CDS start is untouched.
        assert_eq!(alignment.cds_start, Some(2300));
        assert_eq!(alignment.cds_end, Some(2598));

        Ok(())
    }

    /// The bases that `fix_cds` pads exist in the transcript only. The alignment must
    /// therefore keep its genomic length, and the first transcript base must keep its
    /// position.
    #[rstest::rstest]
    #[case("T3P", 1)]
    #[case("T3M", -1)]
    fn fix_cds_pads_the_transcript_only(
        #[case] tx_id: &str,
        #[case] strand: i16,
    ) -> Result<(), anyhow::Error> {
        let mut loader = load(GFF3_CDS_END_NF)?;
        loader.fix_cds();

        let tx = loader
            .transcript_id_to_transcript
            .get(&TranscriptId::try_new(tx_id)?)
            .context("transcript not loaded")?;
        let alignment = tx
            .genome_builds
            .get("GRCh38")
            .context("no GRCh38 alignment")?;
        assert_eq!(tx.cds_length(), Some(192));

        // The alignment as the hgvs mapper sees it, see `Provider::get_tx_exons`.
        let exons = alignment
            .exons
            .iter()
            .map(|exon| TxExonsRecord {
                alt_start_i: exon.alt_start_i,
                alt_end_i: exon.alt_end_i,
                cigar: exon.cigar.clone(),
                ..Default::default()
            })
            .collect::<Vec<_>>();
        let mapper = CigarMapper::new(&build_tx_cigar(&exons, strand)?);

        // Two exons of 100 bases around an intron of 200 bases, plus 2 padding bases.
        assert_eq!(mapper.ref_len, 400);
        assert_eq!(mapper.tgt_len, 202);
        // The first transcript base is the first genomic base on `+` and the last one on
        // `-`, where the mapper counts transcript positions from the genomic start.
        let (ref_pos, tgt_pos) = if strand == 1 { (0, 0) } else { (399, 201) };
        assert_eq!(mapper.map_ref_to_tgt(ref_pos, "start", true)?.pos, tgt_pos);

        Ok(())
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
        load_gff3(&mut loader, &path, &NoProgress)?;

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
    fn exons_are_stored_in_ascending_genomic_order() -> Result<(), anyhow::Error> {
        let loader = load(GFF3_THREE_EXONS)?;

        let plus_tx = loader
            .transcript_id_to_transcript
            .get(&TranscriptId::try_new("T1")?)
            .context("transcript T1 not loaded")?;
        let plus_exons = &plus_tx
            .genome_builds
            .get("GRCh38")
            .context("T1 has no GRCh38 alignment")?
            .exons;
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
            .get(&TranscriptId::try_new("T2")?)
            .context("transcript T2 not loaded")?;
        let minus_exons = &minus_tx
            .genome_builds
            .get("GRCh38")
            .context("T2 has no GRCh38 alignment")?
            .exons;
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

        Ok(())
    }

    fn transcript<'a>(
        loader: &'a TranscriptLoader,
        tx_id: &str,
    ) -> Result<(&'a Transcript, &'a GenomeAlignment), anyhow::Error> {
        let tx = loader
            .transcript_id_to_transcript
            .get(&TranscriptId::try_new(tx_id)?)
            .with_context(|| format!("transcript {tx_id} not loaded"))?;
        let alignment = tx
            .genome_builds
            .get("GRCh38")
            .with_context(|| format!("{tx_id} has no GRCh38 alignment"))?;
        Ok((tx, alignment))
    }

    /// Each exon as `(alt_start_i, alt_end_i, ord, alt_cds_start_i, alt_cds_end_i, cigar)`.
    #[rstest::rstest]
    #[case("NM_000001.1", vec![
        (1000, 1100, 0, 1, 100, "100M"),
        (1200, 1300, 1, 101, 198, "50=2I48="),
    ])]
    #[case("NM_000002.1", vec![
        (2000, 2100, 1, 104, 203, "100M"),
        (2200, 2300, 0, 1, 103, "40=3D60="),
    ])]
    #[case("NM_000003.1", vec![
        (3000, 3100, 0, 2, 101, "100M"),
        (3200, 3300, 1, 132, 231, "100M"),
    ])]
    fn cdna_match_rows_define_the_exons(
        #[case] tx_id: &str,
        #[case] expected: Vec<(i32, i32, i32, i32, i32, &str)>,
    ) -> Result<(), anyhow::Error> {
        let loader = load(GFF3_REFSEQ)?;
        let (_, alignment) = transcript(&loader, tx_id)?;

        let exons = alignment
            .exons
            .iter()
            .map(|e| {
                (
                    e.alt_start_i,
                    e.alt_end_i,
                    e.ord,
                    e.alt_cds_start_i,
                    e.alt_cds_end_i,
                    e.cigar.as_str(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(exons, expected);

        Ok(())
    }

    /// The codon positions count every transcript base, including bases that align nowhere.
    #[rstest::rstest]
    #[case("NM_000001.1", (1010, 1289), (10, 187))]
    #[case("NM_000002.1", (2049, 2290), (10, 154))]
    #[case("NM_000003.1", (3010, 3251), (11, 182))]
    fn codon_positions_follow_the_alignment(
        #[case] tx_id: &str,
        #[case] genomic_cds: (i32, i32),
        #[case] codons: (i32, i32),
    ) -> Result<(), anyhow::Error> {
        let loader = load(GFF3_REFSEQ)?;
        let (tx, alignment) = transcript(&loader, tx_id)?;

        assert_eq!(
            (alignment.cds_start, alignment.cds_end),
            (Some(genomic_cds.0), Some(genomic_cds.1))
        );
        assert_eq!(
            (tx.start_codon, tx.stop_codon),
            (Some(codons.0), Some(codons.1))
        );

        Ok(())
    }

    /// The hgvs mapper must find the 3 transcript-only bases of the minus-strand
    /// NM_000002.1 40 bases after its transcript start, not 40 bases after its genomic start.
    #[test]
    fn gap_is_read_in_transcript_direction_on_minus_strand() -> Result<(), anyhow::Error> {
        let loader = load(GFF3_REFSEQ)?;
        let (_, alignment) = transcript(&loader, "NM_000002.1")?;

        let exons = alignment
            .exons
            .iter()
            .map(|exon| TxExonsRecord {
                alt_start_i: exon.alt_start_i,
                alt_end_i: exon.alt_end_i,
                cigar: exon.cigar.clone(),
                ..Default::default()
            })
            .collect::<Vec<_>>();
        let mapper = CigarMapper::new(&build_tx_cigar(&exons, -1)?);

        assert_eq!(mapper.ref_len, 300);
        assert_eq!(mapper.tgt_len, 203);
        // Genomic base 2259 is transcript base 43, the first one after the insertion. On `-`
        // the mapper counts transcript positions from the genomic start: 202 - 43 = 159.
        assert_eq!(mapper.map_ref_to_tgt(259, "start", true)?.pos, 159);

        Ok(())
    }

    /// A transcript with a second location (e.g. PAR on chrY) keeps its first one, as in cdot.
    #[rstest::rstest]
    #[case(GFF3_REFSEQ, "NM_000001.1", "NC_000001.11", 2)]
    #[case(GFF3_GENCODE_ENSEMBL, "ENST00000399012.6", "chrX", 1)]
    fn transcript_keeps_its_first_location(
        #[case] gff3: &str,
        #[case] tx_id: &str,
        #[case] contig: &str,
        #[case] n_exons: usize,
    ) -> Result<(), anyhow::Error> {
        let loader = load(gff3)?;
        let (_, alignment) = transcript(&loader, tx_id)?;

        assert_eq!(alignment.contig, contig);
        assert_eq!(alignment.exons.len(), n_exons);

        Ok(())
    }

    #[test]
    fn refseq_attributes_are_read_like_cdot() -> Result<(), anyhow::Error> {
        let loader = load(GFF3_REFSEQ)?;

        let (tx, alignment) = transcript(&loader, "NM_000001.1")?;
        assert_eq!(tx.protein.as_deref(), Some("NP_000001.1"));
        assert_eq!(tx.biotype, Some(vec![BioType::MRna]));
        assert_eq!(tx.partial, None);
        assert_eq!(alignment.tag, Some(vec![Tag::ManeSelect]));
        assert_eq!(
            alignment.note.as_deref(),
            Some("UGA stop codon recoded as selenocysteine")
        );

        let (_, alignment) = transcript(&loader, "NM_000003.1")?;
        assert_eq!(alignment.tag, Some(vec![Tag::RefSeqSelect]));

        transcript(&loader, "NR_000005.1")?;
        assert!(
            !loader
                .transcript_id_to_transcript
                .contains_key(&TranscriptId::try_new("MIR1")?)
        );

        // Named as by cdot, plus the version that `load_cdot` adds.
        let (tx, _) = transcript(&loader, "fake-rna-ND1.0")?;
        assert_eq!(tx.protein.as_deref(), Some("YP_003024026.1"));

        let (tx, alignment) = transcript(&loader, "NR_000004.1")?;
        assert_eq!(tx.protein, None);
        assert_eq!(tx.biotype, Some(vec![BioType::LncRna, BioType::NcRna]));
        assert_eq!(tx.partial, Some(1));
        assert_eq!(alignment.tag, None);
        assert_eq!(
            loader.gene_id_to_gene[&GeneId::Gene("44".into())].biotype,
            Some(vec![BioType::LncRna, BioType::NcRna])
        );

        Ok(())
    }

    #[test]
    fn gencode_and_ensembl_attributes_are_read_like_cdot() -> Result<(), anyhow::Error> {
        let loader = load(GFF3_GENCODE_ENSEMBL)?;

        let (_, alignment) = transcript(&loader, "ENST00000399012.6")?;
        assert_eq!(
            alignment.tag,
            Some(vec![Tag::Basic, Tag::EnsemblCanonical, Tag::ManeSelect])
        );

        assert_eq!(
            loader.gene_id_to_gene[&GeneId::Gene("ENSG00000100001.3".into())].biotype,
            Some(vec![BioType::ProteinCoding, BioType::MRna, BioType::NcRna])
        );

        Ok(())
    }

    /// The biotype attribute of a transcript adds to its biotypes, and the filter flags an
    /// NMD transcript for its biotype.
    #[rstest::rstest]
    #[case::ensembl("ENST00000100002.1", vec![BioType::MRna, BioType::NmdTranscriptVariant])]
    #[case::gencode("ENST00000100003.1", vec![BioType::NmdTranscriptVariant, BioType::MRna])]
    fn nmd_biotype_comes_from_the_attributes(
        #[case] tx_id: &str,
        #[case] biotypes: Vec<BioType>,
    ) -> Result<(), anyhow::Error> {
        let mut loader = load(GFF3_NMD)?;
        assert_eq!(transcript(&loader, tx_id)?.0.biotype, Some(biotypes));

        filter_transcripts(&mut loader)?;
        let reason = loader.discards[&Identifier::Transcript(TranscriptId::try_new(tx_id)?)];
        assert!(reason.contains(Reason::Biotype), "{reason:?}");

        Ok(())
    }

    /// Genes keep the ID of the annotation, even if it names an HGNC ID.
    #[rstest::rstest]
    #[case::refseq_ncbi_gene_id(GFF3_REFSEQ, "NM_000001.1", "11")]
    #[case::refseq_without_hgnc(GFF3_REFSEQ, "NM_000003.1", "33")]
    #[case::gencode_gene_id(GFF3_GENCODE_ENSEMBL, "ENST00000399012.6", "ENSG00000182378.15")]
    #[case::ensembl_gene_id_and_version(
        GFF3_GENCODE_ENSEMBL,
        "ENST00000100001.2",
        "ENSG00000100001.3"
    )]
    fn gene_id_comes_from_the_annotation(
        #[case] gff3: &str,
        #[case] tx_id: &str,
        #[case] gene_id: &str,
    ) -> Result<(), anyhow::Error> {
        let loader = load(gff3)?;

        assert!(
            loader.gene_id_to_transcript_ids[&GeneId::Gene(gene_id.into())]
                .contains(&TranscriptId::try_new(tx_id)?)
        );

        Ok(())
    }

    /// A malformed `cDNA_match` row fails the load, and the error names the value.
    #[rstest::rstest]
    #[case("Target=NM_000001.1 1 x +", "invalid Target \"NM_000001.1 1 x +\"")]
    #[case(
        "Target=NM_000001.1 1 100 +;Gap=M50 Dx",
        "invalid operation \"Dx\" in Gap \"M50 Dx\""
    )]
    fn invalid_cdna_match_values_are_named(#[case] attributes: &str, #[case] message: &str) {
        let gff3 = format!(
            "##gff-version 3\nchr1\tRefSeq\tcDNA_match\t1\t100\t.\t+\t.\tID=aln-1;{attributes}\n"
        );
        let error = load(&gff3).map(|_| ()).unwrap_err();
        assert!(format!("{error:#}").contains(message), "{error:#}");
    }

    /// A gene segment names the Gene ID of the gene above it, but must not replace that gene.
    #[test]
    fn gene_segment_does_not_replace_its_gene() -> Result<(), anyhow::Error> {
        let loader = load(GFF3_REFSEQ)?;

        let gene = &loader.gene_id_to_gene[&GeneId::Gene("3514".into())];
        assert_eq!(gene.gene_symbol.as_deref(), Some("IGKC"));
        assert_eq!(
            gene.description.as_deref(),
            Some("immunoglobulin kappa constant")
        );

        Ok(())
    }
}
