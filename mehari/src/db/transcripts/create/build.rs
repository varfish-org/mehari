use crate::db::transcripts::create::cdot_models;
use crate::db::transcripts::create::models::{
    GeneId, Identifier, Reason, TranscriptExt, TranscriptId, TranscriptLoader,
};
use crate::pbs::txs::{GenomeBuild, SourceVersion, TxSeqDatabase};
use anyhow::Error;
use hgvs::data::cdot::json::models::{BioType, Gene, GenomeAlignment, Tag, Transcript};
use hgvs::sequences::aa3_to_aa1;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use std::ops::Not;
use std::time::Instant;

/// Perform protobuf file construction.
///
/// This can be done by simply converting the models from ``hgvs-rs`` to the prost generated data structures.
pub(crate) fn build_protobuf(
    loader: &mut TranscriptLoader,
    sequence_map: &mut HashMap<TranscriptId, String>,
    source_version: SourceVersion,
) -> Result<TxSeqDatabase, Error> {
    tracing::info!("Constructing protobuf data structures …");
    let start = Instant::now();

    let (gene_ids, tx_ids): (Vec<GeneId>, Vec<TranscriptId>) = {
        // ensure we use all gene_ids that are available to us;
        // this includes gene_ids for which we only have transcripts but no genes or vice versa
        let gene_ids = loader
            .gene_id_to_gene
            .keys()
            .cloned()
            .chain(loader.gene_id_to_transcript_ids.keys().cloned())
            .sorted_unstable()
            .dedup()
            .collect_vec();

        // We do, however, discard transcripts that have no sequence or gene id.
        let tx_ids = loader
            .gene_id_to_transcript_ids
            .values()
            .flatten()
            .sorted_unstable()
            .dedup()
            .filter(|&tx_id| {
                loader
                    .discards
                    .get(&Identifier::Transcript(tx_id.clone()))
                    .is_none_or(|reason| !reason.intersects(Reason::hard()))
            })
            .cloned()
            .collect();
        (gene_ids, tx_ids)
    };

    // Construct sequence database.
    tracing::info!("  Constructing sequence database …");
    let seq_db = {
        // filter out transcripts that have no sequence, enumerate, collect; build SequenceDb
        let (aliases, aliases_idx, seqs) = tx_ids
            .iter()
            .filter_map(|tx_id| {
                sequence_map
                    .remove(tx_id)
                    .map(|seq| ((*tx_id).to_string(), seq))
            })
            .enumerate()
            .map(|(idx, (tx_id, seq))| (tx_id, idx as u32, seq))
            .multiunzip();

        crate::pbs::txs::SequenceDb {
            aliases,
            aliases_idx,
            seqs,
        }
    };

    tracing::info!("  Creating transcript records for each gene…");
    let empty_txs = vec![];
    let data_transcripts = gene_ids
        .iter()
        .flat_map(|gene_id| {
            let tx_ids = loader
                .gene_id_to_transcript_ids
                .get(gene_id)
                .unwrap_or(&empty_txs);
            tx_ids
                .iter()
                .sorted_unstable()
                .map(|tx_id| protobuf_transcript(loader, gene_id, tx_id))
        })
        .collect();

    tracing::info!(" … done creating transcripts in {:#?}", start.elapsed());

    // Build mapping of gene HGNC symbol to transcript IDs.
    tracing::info!("  Build gene symbol to transcript ID mapping …");
    let empty = vec![];
    let gene_to_tx = gene_ids
        .iter()
        .map(|gene_id| {
            let tx_ids = loader
                .gene_id_to_transcript_ids
                .get(gene_id)
                .unwrap_or(&empty);
            let gene_id_reason = loader.discards.get(&Identifier::Gene(gene_id.clone()));
            let filtered = gene_id_reason.is_some_and(|reason| reason.intersects(Reason::hard()));
            let hgnc_id = loader.gene_id_to_gene.get(gene_id).and_then(hgnc_id);
            crate::pbs::txs::GeneToTxId {
                #[allow(deprecated)]
                legacy_gene_id: legacy_gene_id(loader, gene_id, tx_ids, hgnc_id.as_deref()),
                tx_ids: tx_ids
                    .iter()
                    .sorted_unstable()
                    .map(|tx_id| tx_id.to_string())
                    .collect(),
                filtered: Some(filtered),
                filter_reason: gene_id_reason.map(|r| r.bits()),
                gene_id: gene_id.to_string(),
                hgnc_id,
            }
        })
        .collect::<Vec<_>>();
    tracing::info!(" … done building gene symbol to transcript ID mapping");

    // Compose transcript database from transcripts and gene to transcript mapping.
    tracing::info!("  Composing transcript seq database …");
    let tx_db = crate::pbs::txs::TranscriptDb {
        transcripts: data_transcripts,
        gene_to_tx,
    };

    // Compose the final transcript and sequence database.
    let tx_seq_db = TxSeqDatabase {
        tx_db: Some(tx_db),
        seq_db: Some(seq_db),
        version: Some(crate::common::version().to_string()),
        source_version: vec![source_version],
    };

    tracing::info!(" … done composing transcript seq database");

    Ok(tx_seq_db)
}

/// The HGNC ID of `gene`, e.g. `HGNC:1100`.
fn hgnc_id(gene: &Gene) -> Option<String> {
    gene.hgnc.as_ref().map(|hgnc| format!("HGNC:{hgnc}"))
}

/// The deprecated `legacy_gene_id` of a gene. It is the HGNC ID of the gene if a transcript
/// of the gene has an HGNC ID as `legacy_gene_id`, as cdot transcripts do. Otherwise it is the
/// greatest value of the transcripts with the highest gene version, so GFF3 genes keep `GENE:`
/// and the gene ID with version. A gene without transcripts gets the value that its loader
/// records, else `GENE:` and the gene ID.
fn legacy_gene_id(
    loader: &TranscriptLoader,
    gene_id: &GeneId,
    tx_ids: &[TranscriptId],
    hgnc_id: Option<&str>,
) -> String {
    let tx_legacy_gene_ids = tx_ids
        .iter()
        .filter_map(|tx_id| loader.legacy_gene_ids.get(tx_id))
        .collect_vec();
    let named_by_hgnc_id = tx_legacy_gene_ids
        .iter()
        .any(|(_, legacy_gene_id)| legacy_gene_id.starts_with("HGNC:"));
    match hgnc_id {
        Some(hgnc_id) if named_by_hgnc_id => hgnc_id.to_string(),
        _ => tx_legacy_gene_ids
            .into_iter()
            .max()
            .map(|(_, legacy_gene_id)| legacy_gene_id)
            .or_else(|| loader.gene_legacy_gene_ids.get(gene_id))
            .cloned()
            .unwrap_or_else(|| format!("GENE:{gene_id}")),
    }
}

fn protobuf_transcript(
    loader: &TranscriptLoader,
    gene_id: &GeneId,
    tx_id: &TranscriptId,
) -> crate::pbs::txs::Transcript {
    let tx = loader
        .transcript_id_to_transcript
        .get(tx_id)
        .unwrap_or_else(|| panic!("No transcript for id {:?}", tx_id));

    // Combine discard reasons for transcript and gene level.
    let tx_reason = loader
        .discards
        .get(&Identifier::Transcript(tx_id.clone()))
        .copied()
        .unwrap_or_default();
    let gene_reason = loader
        .discards
        .get(&Identifier::Gene(gene_id.clone()))
        .copied()
        .unwrap_or_default();
    let combined_reason = tx_reason | gene_reason;

    // Only filter out the transcript if it intersects with a "hard" reason.
    // Soft-filtered transcripts will have `filtered = false`, but their `filter_reason`
    // will still contain the bitmask for downstream inspection.
    let filtered = combined_reason.intersects(Reason::hard());

    // ... build genome alignment for selected:
    let (genome_alignments, tags): (Vec<_>, Vec<_>) = tx
        .genome_builds
        .iter()
        .map(|(genome_build, alignment)| protobuf_genome_alignment(genome_build, alignment))
        .unzip();

    let seleno_tag = tx
        .biotype
        .as_ref()
        .is_some_and(|bt| bt.contains(&BioType::Selenoprotein))
        .then(|| crate::pbs::txs::TranscriptTag::Selenoprotein.into());

    let tags = tags
        .into_iter()
        .flatten()
        .chain(seleno_tag)
        .sorted_unstable()
        .dedup()
        .collect();

    // Now, just obtain the basic properties and create a new `data::Transcript`.
    let gene = loader
        .gene_id_to_gene
        .get(gene_id)
        .cloned()
        .unwrap_or_else(|| Gene {
            aliases: None,
            biotype: None,
            description: None,
            gene_symbol: None,
            hgnc: None,
            map_location: None,
            summary: None,
            url: "".to_string(),
        });
    let hgnc_id = hgnc_id(&gene);

    let Transcript {
        protein,
        start_codon,
        stop_codon,
        ..
    } = tx.clone();
    let biotype = if tx.protein_coding() {
        crate::pbs::txs::TranscriptBiotype::Coding.into()
    } else {
        crate::pbs::txs::TranscriptBiotype::NonCoding.into()
    };

    crate::pbs::txs::Transcript {
        id: (*tx_id).to_string(),
        gene_symbol: gene.gene_symbol.unwrap_or("<MISSING>".to_string()),
        #[allow(deprecated)]
        legacy_gene_id: loader
            .legacy_gene_ids
            .get(tx_id)
            .map(|(_, legacy_gene_id)| legacy_gene_id.clone())
            .unwrap_or_else(|| format!("GENE:{gene_id}")),
        biotype,
        tags,
        protein,
        start_codon,
        stop_codon,
        genome_alignments,
        filtered: Some(filtered),
        filter_reason: combined_reason
            .is_empty()
            .not()
            .then(|| combined_reason.bits()),
        translation_exceptions: translation_exceptions(tx),
        gene_id: gene_id.to_string(),
        hgnc_id,
    }
}

/// The one-letter code of an amino acid in `transl_except`, as hgvs-rs reads cdot: `Other` is
/// an unspecified amino acid (`X`) and `TERM` a stop codon (`*`).
fn transl_except_amino_acid(name: &str) -> Option<char> {
    match name {
        "Other" => Some('X'),
        "TERM" => Some('*'),
        _ if name.len() == 3 => aa3_to_aa1(name).ok()?.chars().next(),
        _ => None,
    }
}

/// The translation exceptions of `tx`, ordered by position. An amino acid without one-letter
/// code is skipped with a warning.
fn translation_exceptions(tx: &Transcript) -> Vec<crate::pbs::txs::TranslationException> {
    let mut exceptions = Vec::new();
    let transl_except = tx
        .translation
        .as_ref()
        .and_then(|translation| translation.transl_except.as_ref());
    for (name, positions) in transl_except.into_iter().flatten() {
        let Some(amino_acid) = transl_except_amino_acid(name) else {
            tracing::warn!(
                "skipping transl_except {name} of {}: unknown amino acid",
                tx.id
            );
            continue;
        };
        exceptions.extend(positions.iter().map(|&position| {
            crate::pbs::txs::TranslationException {
                position,
                amino_acid: amino_acid.to_string(),
            }
        }));
    }
    exceptions.sort_by_key(|exception| exception.position);
    exceptions
}

fn protobuf_genome_alignment(
    genome_build: &str,
    alignment: &GenomeAlignment,
) -> (crate::pbs::txs::GenomeAlignment, HashSet<i32>) {
    let mut tags = HashSet::new();
    let genome_build = genome_build.to_lowercase();

    let GenomeAlignment {
        contig,
        cds_start,
        cds_end,
        ..
    } = alignment.clone();
    let strand = match alignment.strand {
        cdot_models::Strand::Plus => crate::pbs::txs::Strand::Plus,
        cdot_models::Strand::Minus => crate::pbs::txs::Strand::Minus,
    };
    if let Some(tag) = alignment.tag.as_ref() {
        for t in tag {
            let elem = match t {
                Tag::Basic => crate::pbs::txs::TranscriptTag::Basic.into(),
                Tag::EnsemblCanonical => crate::pbs::txs::TranscriptTag::EnsemblCanonical.into(),
                Tag::ManeSelect => crate::pbs::txs::TranscriptTag::ManeSelect.into(),
                Tag::ManePlusClinical => crate::pbs::txs::TranscriptTag::ManePlusClinical.into(),
                Tag::RefSeqSelect => crate::pbs::txs::TranscriptTag::RefSeqSelect.into(),
                Tag::GencodePrimary => crate::pbs::txs::TranscriptTag::GencodePrimary.into(),
                Tag::Other(v) => match v.as_str() {
                    "EnsemblGraft" => crate::pbs::txs::TranscriptTag::EnsemblGraft.into(),
                    "basic-backport" => crate::pbs::txs::TranscriptTag::BasicBackport.into(),
                    "ensembl_canonical-backport" => {
                        crate::pbs::txs::TranscriptTag::EnsemblCanonicalBackport.into()
                    }
                    "mane_select-backport" => {
                        crate::pbs::txs::TranscriptTag::ManeSelectBackport.into()
                    }
                    "mane_plus_clinical-backport" => {
                        crate::pbs::txs::TranscriptTag::ManePlusClinicalBackport.into()
                    }
                    "ref_seq_select-backport" => {
                        crate::pbs::txs::TranscriptTag::RefSeqSelectBackport.into()
                    }
                    "selenoprotein-backport" => {
                        crate::pbs::txs::TranscriptTag::SelenoproteinBackport.into()
                    }
                    "gencode_primary-backport" => {
                        crate::pbs::txs::TranscriptTag::GencodePrimaryBackport.into()
                    }
                    x if x.ends_with("-backport") => {
                        crate::pbs::txs::TranscriptTag::OtherBackport.into()
                    }
                    _ => crate::pbs::txs::TranscriptTag::Other.into(),
                },
            };
            tags.insert(elem);
        }
    }
    // Look into any "note" string for a selenoprotein marker and
    // add this as a tag.
    if let Some(note) = alignment.note.as_ref() {
        let needle = "UGA stop codon recoded as selenocysteine";
        if note.contains(needle) {
            tags.insert(crate::pbs::txs::TranscriptTag::Selenoprotein.into());
        }
    }
    // and construct vector of all exons
    let exons: Vec<_> = alignment
        .exons
        .iter()
        .map(|exon| {
            let cdot_models::Exon {
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

    let genome_build_enum = match genome_build.to_lowercase().as_str() {
        "grch37" | "grch37p10" => GenomeBuild::Grch37,
        "grch38" => GenomeBuild::Grch38,
        _ => GenomeBuild::Unknown,
    };
    let genome_alignment = crate::pbs::txs::GenomeAlignment {
        #[allow(deprecated)]
        genome_build_enum: genome_build_enum.into(),
        genome_build,
        contig,
        cds_start,
        cds_end,
        strand: strand.into(),
        exons,
    };
    (genome_alignment, tags)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use anyhow::Context;

    /// `gene_id`, `hgnc_id` and `legacy_gene_id` of a transcript or a gene in the database.
    pub(crate) type GeneFields = (String, Option<String>, String);

    /// The gene fields of transcript `tx_id` and of its gene, as the database holds them.
    #[allow(deprecated)]
    pub(crate) fn gene_fields(
        loader: &mut TranscriptLoader,
        tx_id: &str,
    ) -> Result<[GeneFields; 2], Error> {
        let db = build_protobuf(loader, &mut HashMap::new(), Default::default())?;
        let tx_db = db.tx_db.context("no transcript database")?;
        let tx = tx_db
            .transcripts
            .iter()
            .find(|tx| tx.id == tx_id)
            .with_context(|| format!("{tx_id} not in the database"))?;
        let gene = tx_db
            .gene_to_tx
            .iter()
            .find(|gene| gene.tx_ids.contains(&tx.id))
            .with_context(|| format!("no gene of {tx_id} in the database"))?;
        Ok([
            (
                tx.gene_id.clone(),
                tx.hgnc_id.clone(),
                tx.legacy_gene_id.clone(),
            ),
            (
                gene.gene_id.clone(),
                gene.hgnc_id.clone(),
                gene.legacy_gene_id.clone(),
            ),
        ])
    }

    /// The gene fields of the gene `gene_id`, as the database holds them.
    #[allow(deprecated)]
    pub(crate) fn gene_to_tx_fields(
        loader: &mut TranscriptLoader,
        gene_id: &str,
    ) -> Result<GeneFields, Error> {
        let db = build_protobuf(loader, &mut HashMap::new(), Default::default())?;
        let gene = db
            .tx_db
            .context("no transcript database")?
            .gene_to_tx
            .into_iter()
            .find(|gene| gene.gene_id == gene_id)
            .with_context(|| format!("{gene_id} not in the database"))?;
        Ok((gene.gene_id, gene.hgnc_id, gene.legacy_gene_id))
    }

    /// A transcript whose gene has no record gets the symbol `<MISSING>` and no HGNC ID.
    #[test]
    fn transcript_without_gene_record_is_kept() -> Result<(), Error> {
        let mut loader = TranscriptLoader::new("GRCh38".to_string(), false);
        let tx_id = TranscriptId::try_new("NM_000001.1")?;
        loader.transcript_id_to_transcript.insert(
            tx_id.clone(),
            serde_json::from_value(serde_json::json!({
                "id": "NM_000001.1",
                "gene_version": "G1",
                "genome_builds": {},
            }))?,
        );
        loader
            .gene_id_to_transcript_ids
            .insert(GeneId("G1".to_string()), vec![tx_id]);

        let db = build_protobuf(&mut loader, &mut HashMap::new(), Default::default())?;
        let tx = &db.tx_db.context("no transcript database")?.transcripts[0];

        assert_eq!(
            (
                tx.gene_symbol.as_str(),
                tx.gene_id.as_str(),
                tx.hgnc_id.as_deref()
            ),
            ("<MISSING>", "G1", None)
        );
        Ok(())
    }
}
