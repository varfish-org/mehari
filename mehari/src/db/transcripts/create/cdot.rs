use crate::common::progress::{Progress, open_with_progress};
use crate::db::transcripts::create::cdot_models;
use crate::db::transcripts::create::models::{
    GeneId, TranscriptId, TranscriptLoader, split_version,
};
use anyhow::Error;
use hgvs::data::cdot::json::models::{Gene, Transcript};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet, HashMap};
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
    // In order of ID, so that the same transcript makes the gene of transcripts without gene
    // record in every run.
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
        .collect::<Result<BTreeMap<_, _>, _>>()?;
    let codons_count_unaligned = codons_count_unaligned_bases(&cdot_version);
    loader.annotation_version = cdot_version;

    // The HGNC IDs that the versions and the transcripts of each gene name.
    let mut hgnc_ids: HashMap<GeneId, BTreeSet<String>> = HashMap::new();

    // All versions of a gene form one gene, with the record of the highest version.
    let mut gene_versions = HashMap::new();
    for (gene_id, gene) in cdot_genes {
        let (gene_id, version) = split_version(&gene_id);
        let gene_id = GeneId(gene_id.to_string());
        if let Some(hgnc) = &gene.hgnc {
            hgnc_ids
                .entry(gene_id.clone())
                .or_default()
                .insert(hgnc.clone());
        }
        if gene_versions
            .get(&gene_id)
            .is_some_and(|&highest| highest > version)
        {
            continue;
        }
        gene_versions.insert(gene_id.clone(), version);
        loader.gene_id_to_gene.insert(gene_id.clone(), gene);
        loader.gene_id_to_transcript_ids.entry(gene_id).or_default();
    }

    for (tx_id, mut tx) in cdot_transcripts {
        shift_codons_past_unaligned_bases(&mut tx, &loader.genome_release, codons_count_unaligned);
        let (gene_id, key) = gene_key(&tx);
        if let Some(hgnc) = &tx.hgnc {
            hgnc_ids
                .entry(gene_id.clone())
                .or_default()
                .insert(hgnc.clone());
        }
        loader
            .transcript_gene_ids
            .insert(tx_id.clone(), key.to_string());
        tx.gene_version = gene_id.to_string();

        // only insert the dummy gene if we haven't seen this gene ID before
        loader
            .gene_id_to_gene
            .entry(gene_id.clone())
            .or_insert_with(|| Gene {
                hgnc: tx.hgnc.clone(),
                gene_symbol: tx.gene_name.clone(),
                aliases: None,
                biotype: None,
                description: None,
                map_location: None,
                summary: None,
                url: "".into(),
            });
        loader
            .gene_id_to_transcript_ids
            .entry(gene_id)
            .or_default()
            .push(tx_id.clone());
        loader.transcript_id_to_transcript.insert(tx_id, tx);
    }

    // The database gives all transcripts of a gene the HGNC ID of the gene record.
    for (gene_id, mut others) in hgnc_ids {
        let hgnc = loader
            .gene_id_to_gene
            .get(&gene_id)
            .and_then(|gene| gene.hgnc.as_ref());
        if let Some(hgnc) = hgnc {
            others.remove(hgnc);
        }
        if !others.is_empty() {
            tracing::warn!(
                "gene {gene_id} has HGNC ID {}, but its transcripts or other versions name {}",
                hgnc.map_or("none".to_string(), |hgnc| format!("HGNC:{hgnc}")),
                others.iter().map(|hgnc| format!("HGNC:{hgnc}")).join(", ")
            );
        }
    }

    Ok(())
}

/// The gene of a cdot transcript: its `GeneId`, and the key that `Transcript.gene_id` keeps. The
/// key is `gene_version`, e.g. `ENSG00000121410.14`, else the gene name, else the transcript ID.
/// Only `gene_version` loses its version in the `GeneId`.
fn gene_key(tx: &Transcript) -> (GeneId, &str) {
    if tx.gene_version.is_empty() {
        let key = tx.gene_name.as_deref().unwrap_or(tx.id.as_str());
        (GeneId(key.to_string()), key)
    } else {
        let (gene_id, _) = split_version(&tx.gene_version);
        (GeneId(gene_id.to_string()), &tx.gene_version)
    }
}

/// Whether the cdot data release `cdot_version`, e.g. `0.2.36`, counts the transcript bases that
/// align nowhere in `start_codon` and `stop_codon`, as in the transcript positions of the exons.
/// It does since 0.2.34 (SACGF/cdot#123).
fn codons_count_unaligned_bases(cdot_version: &str) -> bool {
    let version: Vec<u32> = (cdot_version.split('.'))
        .map_while(|part| part.parse().ok())
        .collect();
    version.as_slice() >= [0, 2, 34].as_slice()
}

/// Count the unaligned bases in front of the codons of a cdot transcript, as the exons do.
///
/// Up to data release 0.2.33, cdot sets `start_codon` and `stop_codon` from the transcript
/// lengths of the exons in front of the codon. It leaves out the transcript bases that align
/// nowhere: those in front of the first exon and those between two exons.
/// `TranscriptLoader::fix_unaligned_bases` adds these bases to the exons, so the codons must count
/// them, too. Since 0.2.34, cdot counts them itself, which `codons_count_unaligned` tells.
///
/// A CDS that starts in unaligned bases begins at the first aligned base in cdot, without its
/// frame. If its length is then no multiple of 3, the frame comes from the stop codon, as GFF3
/// gives it with the phase.
fn shift_codons_past_unaligned_bases(
    tx: &mut cdot_models::Transcript,
    genome_release: &str,
    codons_count_unaligned: bool,
) {
    let Some((_, alignment)) = tx
        .genome_builds
        .iter()
        .find(|(build, _)| build.eq_ignore_ascii_case(genome_release))
    else {
        return;
    };
    let mut exons: Vec<_> = alignment.exons.iter().collect();
    exons.sort_by_key(|exon| exon.ord);
    if !codons_count_unaligned {
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
    }

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
    use crate::common::progress::NoProgress;
    use crate::db::transcripts::create::build::build_protobuf;
    use crate::db::transcripts::create::build::tests::{gene_fields, gene_to_tx_fields};
    use crate::db::transcripts::create::filter::filter_transcripts;
    use crate::db::transcripts::create::models::{Identifier, Reason};
    use anyhow::Context;
    use serde_json::{Value, json};
    use std::io::Write;

    fn load(path: impl AsRef<Path>) -> Result<TranscriptLoader, Error> {
        let mut loader = TranscriptLoader::new("GRCh38".to_string(), false);
        load_cdot(&mut loader, path, &NoProgress)?;
        Ok(loader)
    }

    /// A cdot file with `genes` and `transcripts`.
    fn load_json(genes: Value, transcripts: Value) -> Result<TranscriptLoader, Error> {
        let mut file = tempfile::NamedTempFile::new()?;
        let container = json!({
            "cdot_version": "0.2.36",
            "genome_builds": ["GRCh38"],
            "genes": genes,
            "transcripts": transcripts,
        });
        file.write_all(container.to_string().as_bytes())?;
        load(file.path())
    }

    /// A cdot transcript of the gene `gene_version` with one exon.
    fn cdot_tx(id: &str, gene_version: &str, hgnc: Option<&str>) -> Value {
        json!({
            "id": id,
            "gene_version": gene_version,
            "hgnc": hgnc,
            "genome_builds": {
                "GRCh38": {
                    "contig": "NC_000001.11",
                    "strand": "+",
                    "exons": [[1000, 1100, 0, 1, 100, null]],
                },
            },
        })
    }

    /// The database names the gene of a transcript by `gene_version` without version, and by
    /// its HGNC ID.
    #[rstest::rstest]
    #[case::refseq(
        "tests/data/db/create/txs/cdot-0.2.22.refseq.grch37.brca1_opa1.json",
        "NM_007294.4",
        "672",
        Some("HGNC:1100")
    )]
    #[case::ensembl_without_hgnc(
        "tests/data/db/create/mitochondrial/cdot-0.2.23.ensembl.chrMT.grch37.gff3.json",
        "ENST00000387314.1",
        "ENSG00000210049",
        None
    )]
    fn gene_fields_reach_the_database(
        #[case] path: &str,
        #[case] tx_id: &str,
        #[case] gene_id: &str,
        #[case] hgnc_id: Option<&str>,
    ) -> Result<(), Error> {
        let mut loader = load(path)?;
        let expected = (gene_id.to_string(), hgnc_id.map(String::from));
        assert_eq!(
            gene_fields(&mut loader, tx_id)?,
            [expected.clone(), expected]
        );
        Ok(())
    }

    /// A file with several releases names one gene with several versions. They form one gene
    /// with the record of the highest version and the gene ID without version. Each
    /// transcript keeps its own gene version in `gene_id`.
    #[test]
    fn gene_versions_form_one_gene() -> Result<(), Error> {
        let mut loader = load_json(
            json!({
                "ENSG00000000001.9": {"gene_symbol": "OLD", "url": ""},
                "ENSG00000000001.10": {"gene_symbol": "NEW", "url": ""},
            }),
            json!({
                "ENST00000000001.1": cdot_tx("ENST00000000001.1", "ENSG00000000001.9", None),
                "ENST00000000002.1": cdot_tx("ENST00000000002.1", "ENSG00000000001.10", None),
            }),
        )?;
        let db = build_protobuf(&mut loader, &mut HashMap::new(), Default::default())?;
        let tx_db = db.tx_db.context("no transcript database")?;

        let genes: Vec<_> = tx_db
            .gene_to_tx
            .iter()
            .map(|gene| (gene.gene_id.as_str(), gene.tx_ids.clone()))
            .collect();
        assert_eq!(
            genes,
            [(
                "ENSG00000000001",
                vec![
                    "ENST00000000001.1".to_string(),
                    "ENST00000000002.1".to_string()
                ]
            )]
        );
        let txs: Vec<_> = tx_db
            .transcripts
            .iter()
            .map(|tx| (tx.id.as_str(), tx.gene_id.as_str(), tx.gene_symbol.as_str()))
            .collect();
        assert_eq!(
            txs,
            [
                ("ENST00000000001.1", "ENSG00000000001.9", "NEW"),
                ("ENST00000000002.1", "ENSG00000000001.10", "NEW"),
            ]
        );

        Ok(())
    }

    /// All transcripts of a gene get the HGNC ID of the gene record, also if they name
    /// another one or none. `load_cdot` warns about the other one.
    #[tracing_test::traced_test]
    #[test]
    fn transcripts_get_the_hgnc_id_of_their_gene() -> Result<(), Error> {
        let mut loader = load_json(
            json!({"2": {"gene_symbol": "A2M", "hgnc": "7", "url": ""}}),
            json!({
                "NM_000001.1": cdot_tx("NM_000001.1", "2", Some("7")),
                "NM_000002.1": cdot_tx("NM_000002.1", "2", Some("8")),
                "NM_000003.1": cdot_tx("NM_000003.1", "2", None),
            }),
        )?;

        for tx_id in ["NM_000001.1", "NM_000002.1", "NM_000003.1"] {
            let gene = ("2".to_string(), Some("HGNC:7".to_string()));
            assert_eq!(gene_fields(&mut loader, tx_id)?, [gene.clone(), gene]);
        }
        assert!(logs_contain(
            "gene 2 has HGNC ID HGNC:7, but its transcripts or other versions name HGNC:8"
        ));

        Ok(())
    }

    /// A gene without transcripts keeps its gene ID and HGNC ID.
    #[rstest::rstest]
    #[case::hgnc("2", Some("HGNC:7"))]
    #[case::without_hgnc("3", None)]
    fn gene_without_transcripts_keeps_its_gene_fields(
        #[case] gene_id: &str,
        #[case] hgnc_id: Option<&str>,
    ) -> Result<(), Error> {
        let mut loader = load_json(
            json!({
                "2": {"gene_symbol": "A2M", "hgnc": "7", "url": ""},
                "3": {"gene_symbol": "GX", "url": ""},
            }),
            json!({}),
        )?;
        assert_eq!(
            gene_to_tx_fields(&mut loader, gene_id)?,
            (gene_id.to_string(), hgnc_id.map(String::from))
        );

        Ok(())
    }

    /// A transcript without HGNC ID gets the HGNC ID of its gene record, so the filters do
    /// not flag it with `MissingHgncId`.
    #[test]
    fn transcript_without_hgnc_id_in_gene_with_hgnc_id_is_not_flagged() -> Result<(), Error> {
        let mut loader = load_json(
            json!({"2": {"gene_symbol": "A2M", "hgnc": "7", "url": ""}}),
            json!({"NM_000003.1": cdot_tx("NM_000003.1", "2", None)}),
        )?;
        filter_transcripts(&mut loader)?;

        let reason = loader
            .discards
            .get(&Identifier::Transcript(TranscriptId::try_new(
                "NM_000003.1",
            )?))
            .copied()
            .unwrap_or_default();
        assert!(!reason.contains(Reason::MissingHgncId), "{reason:?}");

        Ok(())
    }

    #[test]
    fn highest_gene_version_wins_in_any_order() -> Result<(), Error> {
        let genes: serde_json::Map<String, Value> = (1..=9)
            .map(|version| {
                (
                    format!("ENSG00000000001.{version}"),
                    json!({"gene_symbol": format!("V{version}"), "url": ""}),
                )
            })
            .collect();
        let loader = load_json(Value::Object(genes), json!({}))?;
        let gene = loader
            .gene_id_to_gene
            .get(&GeneId("ENSG00000000001".to_string()))
            .context("no gene")?;
        assert_eq!(gene.gene_symbol.as_deref(), Some("V9"));
        Ok(())
    }

    /// A transcript whose gene is not in `genes` gets a gene from its own fields.
    #[test]
    fn transcript_without_gene_record_makes_a_gene() -> Result<(), Error> {
        let mut tx = cdot_tx("NM_000001.1", "ENSG00000000001.5", Some("7"));
        tx["gene_name"] = json!("GX");
        let loader = load_json(json!({}), json!({"NM_000001.1": tx}))?;

        let gene = &loader.gene_id_to_gene[&GeneId("ENSG00000000001".to_string())];
        assert_eq!(gene.gene_symbol.as_deref(), Some("GX"));
        assert_eq!(gene.hgnc.as_deref(), Some("7"));
        Ok(())
    }

    /// Of several transcripts without gene record, the one with the lowest ID makes the gene,
    /// whatever the order of the file.
    #[test]
    fn lowest_transcript_id_makes_the_gene() -> Result<(), Error> {
        let transcripts: serde_json::Map<String, Value> = (1..=9)
            .map(|i| {
                let id = format!("NM_00000{i}.1");
                let mut tx = cdot_tx(&id, "ENSG00000000001.5", Some(&i.to_string()));
                tx["gene_name"] = json!(format!("G{i}"));
                (id, tx)
            })
            .collect();
        let loader = load_json(json!({}), Value::Object(transcripts))?;

        let gene = &loader.gene_id_to_gene[&GeneId("ENSG00000000001".to_string())];
        assert_eq!(
            (gene.gene_symbol.as_deref(), gene.hgnc.as_deref()),
            (Some("G1"), Some("1"))
        );
        Ok(())
    }

    /// Without `gene_version`, the gene is the gene name, else the transcript ID.
    #[rstest::rstest]
    #[case::gene_name(Some("GX"), "GX")]
    #[case::transcript_id(None, "NM_000001.1")]
    fn gene_key_falls_back(
        #[case] gene_name: Option<&str>,
        #[case] expected: &str,
    ) -> Result<(), Error> {
        let tx: cdot_models::Transcript = serde_json::from_value(json!({
            "id": "NM_000001.1",
            "gene_version": "",
            "gene_name": gene_name,
            "genome_builds": {},
        }))?;
        assert_eq!(gene_key(&tx), (GeneId(expected.to_string()), expected));
        Ok(())
    }

    /// cdot counts the codons in transcript positions. Since data release 0.2.34, these include
    /// the bases that align nowhere, as the transcript positions of the exons do
    /// (SACGF/cdot#123). Up to 0.2.33, they leave these bases out, and the codons move past them.
    /// Either way, the codons end up at the same positions, e.g. with 48 unaligned bases between
    /// two exons, as in NM_001137667.2 (CASP8AP2), or 119 in front of the first exon, as in
    /// NM_053005.5 (MOB2). A CDS that starts in unaligned bases begins at the first aligned base,
    /// without its frame. Its frame then comes from the stop codon: with 5 unaligned bases, it
    /// starts at 5 + (160 - 5) % 3 = 7.
    #[rstest::rstest]
    #[case::between_exons("0.2.36", [(1, 100), (149, 248)], (9, 160), (9, 160))]
    #[case::between_exons_before_0_2_34("0.2.33", [(1, 100), (149, 248)], (9, 112), (9, 160))]
    #[case::in_front("0.2.34", [(6, 105), (106, 205)], (20, 160), (20, 160))]
    #[case::in_front_before_0_2_34("0.2.22", [(6, 105), (106, 205)], (15, 155), (20, 160))]
    #[case::start_in_unaligned_bases("0.2.36", [(6, 105), (106, 205)], (5, 160), (7, 160))]
    #[case::start_in_unaligned_bases_before_0_2_34(
        "0.2.22",
        [(6, 105), (106, 205)],
        (0, 155),
        (7, 160)
    )]
    fn codons_in_transcript_positions(
        #[case] cdot_version: &str,
        #[case] tx_ranges: [(i32, i32); 2],
        #[case] (start_codon, stop_codon): (i32, i32),
        #[case] expected: (i32, i32),
    ) -> Result<(), Error> {
        let [(start_1, end_1), (start_2, end_2)] = tx_ranges;
        let tx = json!({
            "id": "NM_000001.1",
            "gene_version": "1",
            "start_codon": start_codon,
            "stop_codon": stop_codon,
            "genome_builds": {
                "GRCh38": {
                    "contig": "NC_000001.11",
                    "strand": "+",
                    "cds_start": 1000 + expected.0 + 1 - start_1,
                    "cds_end": 2000 + expected.1 - start_2 + 1,
                    "exons": [
                        [1000, 1100, 0, start_1, end_1, null],
                        [2000, 2100, 1, start_2, end_2, null],
                    ],
                },
            },
        });
        let mut file = tempfile::NamedTempFile::new()?;
        let container = json!({
            "cdot_version": cdot_version,
            "genome_builds": ["GRCh38"],
            "genes": {},
            "transcripts": {"NM_000001.1": tx},
        });
        file.write_all(container.to_string().as_bytes())?;
        let loader = load(file.path())?;

        let tx = &loader.transcript_id_to_transcript[&TranscriptId::try_new("NM_000001.1")?];
        assert_eq!(
            (tx.start_codon, tx.stop_codon),
            (Some(expected.0), Some(expected.1))
        );
        Ok(())
    }
}
