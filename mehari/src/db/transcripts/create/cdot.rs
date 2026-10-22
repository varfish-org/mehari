use crate::common::progress::{Progress, open_with_progress};
use crate::db::transcripts::create::cdot_models;
use crate::db::transcripts::create::models::{GeneId, TranscriptId, TranscriptLoader};
use anyhow::Error;
use hgvs::data::cdot::json::models::Gene;
use std::collections::HashMap;
use std::path::Path;

/// Load and extract from cdot JSON.
pub fn load_cdot(
    loader: &mut TranscriptLoader,
    path: impl AsRef<Path>,
    progress: &dyn Progress,
) -> Result<(), Error> {
    let cdot_models::Container {
        genes: cdot_genes,
        transcripts: cdot_transcripts,
        cdot_version,
        ..
    } = read_cdot_json(path.as_ref(), progress)?;
    let cdot_genes = cdot_genes.into_iter().collect::<HashMap<_, _>>();
    let cdot_transcripts = cdot_transcripts
        .into_iter()
        .map(|(txid, mut tx)| {
            let txid = if txid.starts_with("fake-rna-") && !txid.contains('.') {
                format!("{txid}.0")
            } else {
                txid
            };
            tx.id.clone_from(&txid);
            TranscriptId::try_new(txid).map(|t| (t, tx))
        })
        .collect::<Result<HashMap<_, _>, _>>()?;
    loader.annotation_version = cdot_version;

    for (_, gene) in cdot_genes {
        if let Some(gene_id) = gene.hgnc.as_ref() {
            let gene_id: GeneId = gene_id.parse()?;
            loader.gene_id_to_gene.insert(gene_id.clone(), gene);
            loader.gene_id_to_transcript_ids.entry(gene_id).or_default();
        }
    }

    for (tx_id, mut tx) in cdot_transcripts {
        shift_codons_past_unaligned_bases(&mut tx, &loader.genome_release);
        let gene_id = if let Some(hgnc_str) = tx.hgnc.as_ref() {
            Some(hgnc_str.parse()?)
        } else {
            None
        };

        let gene_id = match gene_id {
            Some(id) => id,
            None => {
                // group by gene_name if available, otherwise fall back to transcript id
                let grouping_key = if !tx.gene_version.is_empty() {
                    tx.gene_version.as_str()
                } else {
                    tx.gene_name.as_deref().unwrap_or(tx.id.as_str())
                };
                let fake_id = GeneId::Gene(grouping_key.to_string());

                // only insert the dummy gene if we haven't seen this PseudoId before
                loader
                    .gene_id_to_gene
                    .entry(fake_id.clone())
                    .or_insert_with(|| Gene {
                        hgnc: Some(fake_id.clone().to_string()),
                        gene_symbol: tx.gene_name.clone(),
                        aliases: None,
                        biotype: None,
                        description: None,
                        map_location: None,
                        summary: None,
                        url: "".into(),
                    });
                tx.hgnc = Some(fake_id.clone().to_string());

                fake_id
            }
        };

        loader
            .gene_id_to_transcript_ids
            .entry(gene_id)
            .or_default()
            .push(tx_id.clone());
        loader.transcript_id_to_transcript.insert(tx_id, tx);
    }

    Ok(())
}

/// Count the unaligned bases in front of the codons of a cdot transcript, as the exons do.
///
/// cdot sets `start_codon` and `stop_codon` from the transcript lengths of the exons in front of
/// the codon. It leaves out the transcript bases that align nowhere: those in front of the first
/// exon and those between two exons. `TranscriptLoader::fix_unaligned_bases` adds these bases to
/// the exons, so the codons must count them, too.
///
/// A CDS that starts in unaligned bases begins at the first aligned base in cdot, without its
/// frame. If its length is then no multiple of 3, the frame comes from the stop codon, as GFF3
/// gives it with the phase.
fn shift_codons_past_unaligned_bases(tx: &mut cdot_models::Transcript, genome_release: &str) {
    let Some((_, alignment)) = tx
        .genome_builds
        .iter()
        .find(|(build, _)| build.eq_ignore_ascii_case(genome_release))
    else {
        return;
    };
    let mut exons: Vec<_> = alignment.exons.iter().collect();
    exons.sort_by_key(|exon| exon.ord);
    // The unaligned bases in front of the exon that holds the 0-based transcript position
    // `pos`, which counts the exons as cdot does.
    let unaligned_before = |pos: i32| {
        let mut counted = 0;
        let mut unaligned = 0;
        let mut next_tx_start = 1;
        for exon in &exons {
            unaligned += exon.alt_cds_start_i - next_tx_start;
            counted += exon.alt_cds_end_i - exon.alt_cds_start_i + 1;
            next_tx_start = exon.alt_cds_end_i + 1;
            if pos < counted {
                break;
            }
        }
        unaligned
    };
    tx.start_codon = tx.start_codon.map(|start| start + unaligned_before(start));
    // `stop_codon` is the end of the CDS, so its last base decides.
    tx.stop_codon = tx.stop_codon.map(|stop| stop + unaligned_before(stop - 1));

    let leading = exons.first().map_or(0, |exon| exon.alt_cds_start_i - 1);
    if let (Some(start), Some(stop)) = (tx.start_codon, tx.stop_codon)
        && leading > 0
        && start == leading
    {
        tx.start_codon = Some(start + (stop - start) % 3);
    }
}

pub(crate) fn read_cdot_json(
    path: impl AsRef<Path>,
    progress: &dyn Progress,
) -> Result<cdot_models::Container, Error> {
    let file = open_with_progress(path.as_ref(), progress)?;
    let bar = file.progress.clone();
    let container = if path.as_ref().extension().unwrap_or_default() == "gz" {
        tracing::info!("(from gzip compressed file)");
        serde_json::from_reader(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(
            file,
        )))?
    } else {
        tracing::info!("(from uncompressed file)");
        serde_json::from_reader(std::io::BufReader::new(file))?
    };
    bar.finish();
    Ok(container)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// 5 unaligned bases in front of exon 1 and 5 between exons 1 and 2. cdot counts the start
    /// codon in exon 1 without the first 5 and the end of the CDS in exon 2 without all 10. A CDS
    /// that cdot starts at the first aligned base, at 0, has no frame, so it ends with the stop
    /// codon: it starts at 5 + (160 - 5) % 3 = 7 then.
    #[rstest::rstest]
    #[case::start_in_exon(10, (15, 160))]
    #[case::clipped_start(0, (7, 160))]
    fn codons_count_unaligned_bases(
        #[case] start_codon: i32,
        #[case] expected: (i32, i32),
    ) -> Result<(), Error> {
        let mut tx: cdot_models::Transcript = serde_json::from_value(json!({
            "id": "NM_000001.1",
            "gene_version": "1",
            "start_codon": start_codon,
            "stop_codon": 150,
            "genome_builds": {
                "GRCh38": {
                    "contig": "NC_000001.11",
                    "strand": "+",
                    "cds_start": 1004,
                    "cds_end": 2045,
                    "exons": [
                        [1000, 1100, 0, 6, 105, null],
                        [2000, 2100, 1, 111, 210, null],
                    ],
                },
            },
        }))?;

        shift_codons_past_unaligned_bases(&mut tx, "grch38");

        assert_eq!(
            (tx.start_codon, tx.stop_codon),
            (Some(expected.0), Some(expected.1))
        );
        Ok(())
    }
}
