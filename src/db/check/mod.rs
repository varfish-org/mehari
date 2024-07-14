use crate::annotate::seqvars::load_tx_db;
use anyhow::Result;
use clap::Parser;
use enumflags2::{bitflags, BitFlag, BitFlags};
use hgvs::data::cdot::json::models::Container as Cdot;
use serde;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

/// Command line arguments for `db check` sub command.
#[derive(Parser, Debug)]
#[command(about = "Check transcript database", long_about = None)]
pub struct Args {
    /// Path to database file to check
    #[arg(long)]
    pub path_db: PathBuf,

    /// Paths to the cdot JSON transcripts to check against.
    #[arg(long, required = true)]
    pub path_cdot_json: Vec<PathBuf>,

    #[arg(long)]
    pub path_hgnc_json: PathBuf,

    #[arg(long)]
    pub path_disease_gene_tsv: PathBuf,
}

type HgncId = String;
type GeneSymbol = String;
type OtherId = String;

type IdCollection = (HashSet<HgncId>, HashSet<OtherId>);

fn load_cdot_files(paths: &[PathBuf]) -> Result<IdCollection> {
    let cdot_containers = paths
        .iter()
        .map(crate::db::create::read_cdot_json)
        .collect::<Result<Vec<_>>>()?;
    let Cdot {
        transcripts, genes, ..
    } = cdot_containers
        .into_iter()
        .reduce(|mut a, b| {
            assert_eq!(a.cdot_version, b.cdot_version);
            assert_eq!(a.genome_builds, b.genome_builds);
            a.genes.extend(b.genes);
            a.transcripts.extend(b.transcripts);
            a
        })
        .unwrap();

    Ok((
        genes
            .into_values()
            .flat_map(|g| g.hgnc.map(|hgnc| format!("HGNC:{hgnc}")))
            .collect(),
        transcripts
            .into_values()
            .flat_map(|t| t.hgnc.map(|hgnc| format!("HGNC:{hgnc}")))
            .collect(),
    ))
}

#[derive(Debug, Serialize, Deserialize)]
struct HgncEntry {
    hgnc_id: String,
    location: String,
    location_sortable: String,
    locus_group: String,
    locus_type: String,
    name: String,
    status: String,
    symbol: String,

    mane_select: Option<Vec<String>>,
    ensembl_gene_id: Option<String>,
    entrez_id: Option<String>,
    refseq_accession: Option<Vec<String>>,
    uuid: Option<String>,
    // agr: Option<String>,
    // alias_name: Vec<String>,
    // alias_symbol: Vec<String>,
    // bioparadigms_slc: Option<String>,
    // ccds_id: Option<Vec<String>>,
    // cd: Option<String>,
    // cosmic: Option<String>,
    // date_approved_reserved: Option<String>,
    // date_modified: Option<String>,
    // date_name_changed: Option<String>,
    // date_symbol_changed: Option<String>,
    // ena: Option<Vec<String>>,
    // enzyme_id: Option<String>,
    // gencc: Option<String>,
    // gene_group: Option<Vec<String>>,
    // gene_group_id: Vec<String>,
    // // gtrnadb,
    // homeodb: Option<String>,
    // horde_id: Option<String>,
    // imgt: Option<String>,
    // iuphar: Option<String>,
    // lncipedia: Option<String>,
    // lncrnadb: Option<String>,
    // lsdb: Option<String>,
    // // mamit-trnadb,
    // #[serde(rename = "mamit-trnadb")]
    // mamit_trnadb: Option<String>,
    //
    // merops: Option<String>,
    // mgd_id: Option<Vec<String>>,
    // mirbase: Option<String>,
    // name: String,
    // omim_id: Option<Vec<String>>,
    // orphanet: Option<String>,
    // prev_name: Option<Vec<String>>,
    // prev_symbol: Option<Vec<String>>,
    // #[serde(rename = "pseudogene.org")]
    // pseudogene_org: Option<String>,
    // pubmed_id: Option<Vec<String>>,
    // rgd_id: Option<Vec<String>>,
    // // rna_central_id,
    // snornabase: Option<String>,
    // // symbol_report_tag,
    // ucsc_id: Option<String>,
    // uniprot_ids: Option<Vec<String>>,
    // vega_id: Option<String>,
}

#[bitflags]
#[repr(u8)]
#[derive(Debug, Clone, Copy, Serialize, Hash, PartialEq, Eq)]
enum Annotation {
    NoNcbiGeneId,
    NoNcbiTranscriptId,
    NoEnsemblGeneId,
    NoEnsemblTranscriptId,
}
fn load_hgnc_set(
    path: impl AsRef<Path>,
) -> Result<(HashMap<HgncId, BitFlags<Annotation>>, HashSet<OtherId>)> {
    let hgnc: Value = {
        let reader = File::open(path).map(BufReader::new)?;
        serde_json::from_reader::<_, Value>(reader)?["response"]["docs"].to_owned()
    };
    let entries: Vec<HgncEntry> = serde_json::from_value(hgnc)?;
    let mut external_ids = HashSet::new();
    let mut hgnc_ids = HashMap::new();
    entries.into_iter().for_each(|entry| {
        let hgnc_id = entry.hgnc_id;
        let ncbi_gene_id = entry.entrez_id;
        let ensembl_gene_id = entry.ensembl_gene_id;
        let mut annotation = Annotation::empty();
        if let Some(gene_id) = ncbi_gene_id {
            external_ids.insert(gene_id);
        } else {
            annotation |= Annotation::NoNcbiGeneId;
        }
        if let Some(gene_id) = ensembl_gene_id {
            external_ids.insert(gene_id);
        } else {
            annotation |= Annotation::NoEnsemblGeneId;
        }

        if let Some(refseq_accessions) = entry.refseq_accession {
            if refseq_accessions.is_empty() {
                annotation |= Annotation::NoNcbiTranscriptId;
            }
            external_ids.extend(refseq_accessions);
        }

        if let Some(mane_select) = entry.mane_select {
            if !mane_select.iter().any(|r| r.starts_with("ENST")) {
                annotation |= Annotation::NoEnsemblTranscriptId;
            }
            external_ids.extend(mane_select);
        }

        hgnc_ids.insert(hgnc_id, annotation);
    });
    Ok((hgnc_ids, external_ids))
}

fn load_disease_gene_set(path: impl AsRef<Path>) -> Result<(HashSet<HgncId>, HashSet<GeneSymbol>)> {
    let reader = File::open(path).map(BufReader::new)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    let mut hgnc_ids = HashSet::new();
    let mut gene_symbols = HashSet::new();
    let header = rdr.headers().cloned().unwrap();
    let header = header
        .iter()
        .enumerate()
        .map(|(idx, h)| (h, idx))
        .collect::<HashMap<_, _>>();
    for record in rdr.into_records() {
        let record = record?;
        let gene_symbol = record[header["gene_symbol"]].to_string();
        let hgnc_id = record[header["hgnc_id"]].to_string();
        hgnc_ids.insert(hgnc_id.clone());
        gene_symbols.insert(gene_symbol);
    }
    Ok((hgnc_ids, gene_symbols))
}

pub fn run(_common: &crate::common::Args, args: &Args) -> Result<()> {
    let tx_db = load_tx_db(&format!("{}", args.path_db.display()))?
        .tx_db
        .unwrap();
    let tx_db_hgnc_ids: HashSet<HgncId> =
        tx_db.gene_to_tx.iter().map(|g| g.gene_id.clone()).collect();
    let tx_db_tx_ids: HashSet<OtherId> = tx_db.transcripts.iter().map(|tx| tx.id.clone()).collect();
    let (cdot_hgnc_ids, cdot_tx_ids) = load_cdot_files(&args.path_cdot_json)?;
    let (hgnc_hgnc_ids, hgnc_other_ids) = load_hgnc_set(&args.path_hgnc_json)?;
    let (disease_gene_hgnc_ids, disease_gene_gene_symbols) =
        load_disease_gene_set(&args.path_disease_gene_tsv)?;

    for hgnc_id in cdot_hgnc_ids.difference(&tx_db_hgnc_ids) {
        tracing::warn!("HGNC ID {} not found in transcript database", hgnc_id);
        if let Some(annotations) = hgnc_hgnc_ids.get(hgnc_id) {
            if !annotations.is_empty() {
                tracing::warn!("HGNC ID {} has annotations {:?}", hgnc_id, annotations);
            }
        }
    }

    for hgnc_id in disease_gene_hgnc_ids.difference(&cdot_hgnc_ids) {
        tracing::warn!("HGNC ID {} not found in cdot JSON", hgnc_id);
    }

    Ok(())
}
