use crate::db::create::cdot_models;
use crate::db::create::models::{GeneId, TranscriptId, TranscriptLoader};
use anyhow::Error;
use hgvs::data::cdot::json::models::Gene;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

/// Load and extract from cdot JSON.
pub fn load_cdot(loader: &mut TranscriptLoader, path: impl AsRef<Path>) -> Result<(), Error> {
    let cdot_models::Container {
        genes: cdot_genes,
        transcripts: cdot_transcripts,
        cdot_version,
        ..
    } = read_cdot_json(path.as_ref())?;
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

pub(crate) fn read_cdot_json(path: impl AsRef<Path>) -> Result<cdot_models::Container, Error> {
    Ok(if path.as_ref().extension().unwrap_or_default() == "gz" {
        tracing::info!("(from gzip compressed file)");
        serde_json::from_reader(std::io::BufReader::new(flate2::read::GzDecoder::new(
            File::open(path)?,
        )))?
    } else {
        tracing::info!("(from uncompressed file)");
        serde_json::from_reader(std::io::BufReader::new(File::open(path)?))?
    })
}
