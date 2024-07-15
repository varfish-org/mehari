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
use strum::Display;

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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, PartialOrd, Ord, Display)]
enum Id {
    Hgnc(HgncId),
    GeneSymbol(GeneSymbol),
    NcbiAccession(String),
    NcbiGene(String),
    NcbiTranscript(String),
    EnsemblAccession(String),
    EnsemblGene(String),
    EnsemblTranscript(String),
}

impl From<String> for Id {
    fn from(value: String) -> Self {
        if value.starts_with("HGNC:") {
            Id::Hgnc(value)
        } else if value.starts_with("ENSG") {
            Id::EnsemblGene(value)
        } else if value.starts_with("ENST") {
            Id::EnsemblTranscript(value)
        } else if value.starts_with("NM_")
            || value.starts_with("NR_")
            || value.starts_with("XM_")
            || value.starts_with("XR_")
        {
            Id::NcbiTranscript(value)
        } else if value.parse::<usize>().is_ok() {
            Id::NcbiGene(value)
        } else {
            Id::GeneSymbol(value)
        }
    }
}

type IdCollection = HashSet<Id>;

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

    Ok(genes
        .iter()
        .flat_map(|(gene_id, g)| {
            let gid = Id::from(gene_id.clone());
            let mut ids = vec![gid];
            if let Some(ref hgnc_id) = g.hgnc {
                ids.push(Id::Hgnc(format!("HGNC:{hgnc_id}")));
            }
            if let Some(ref gene_symbol) = g.gene_symbol {
                ids.push(Id::GeneSymbol(gene_symbol.clone()));
            }
            ids
        })
        .chain(transcripts.iter().flat_map(|(tx_id, t)| {
            let tid = Id::from(tx_id.clone());
            let mut ids = vec![tid];
            if let Some(ref hgnc_id) = t.hgnc {
                ids.push(Id::Hgnc(format!("HGNC:{hgnc_id}")));
            }
            if let Some(ref gene_symbol) = t.gene_name {
                ids.push(Id::GeneSymbol(gene_symbol.clone()));
            }
            ids
        }))
        .collect())
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

fn load_hgnc_set(path: impl AsRef<Path>) -> Result<HashMap<Id, BitFlags<Annotation>>> {
    let hgnc: Value = {
        let reader = File::open(path).map(BufReader::new)?;
        serde_json::from_reader::<_, Value>(reader)?["response"]["docs"].to_owned()
    };
    let entries: Vec<HgncEntry> = serde_json::from_value(hgnc)?;
    let mut ids = HashMap::new();
    entries.into_iter().for_each(|entry| {
        let hgnc_id = Id::Hgnc(entry.hgnc_id);
        let ncbi_gene_id = entry.entrez_id.map(Id::NcbiGene);
        let ensembl_gene_id = entry.ensembl_gene_id.map(Id::EnsemblGene);
        let mut annotation = Annotation::empty();
        if let Some(gene_id) = ncbi_gene_id {
            ids.entry(gene_id).or_default();
        } else {
            annotation |= Annotation::NoNcbiGeneId;
        }
        if let Some(gene_id) = ensembl_gene_id {
            ids.entry(gene_id).or_default();
        } else {
            annotation |= Annotation::NoEnsemblGeneId;
        }

        if let Some(refseq_accessions) = entry.refseq_accession {
            if refseq_accessions.is_empty() {
                annotation |= Annotation::NoNcbiTranscriptId;
            }
            ids.extend(
                refseq_accessions
                    .into_iter()
                    .map(|id| (Id::NcbiAccession(id), BitFlags::empty())),
            );
        }

        if let Some(mane_select) = entry.mane_select {
            if !mane_select.iter().any(|r| r.starts_with("ENST")) {
                annotation |= Annotation::NoEnsemblTranscriptId;
            }
            ids.extend(
                mane_select
                    .into_iter()
                    .map(|id| (Id::from(id), BitFlags::empty())),
            )
        }

        *ids.entry(hgnc_id).or_default() |= annotation;
    });
    Ok(ids)
}

fn load_disease_gene_set(path: impl AsRef<Path>) -> Result<HashSet<Id>> {
    let reader = File::open(path).map(BufReader::new)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    let mut ids = HashSet::new();
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
        ids.insert(Id::Hgnc(hgnc_id.clone()));
        ids.insert(Id::GeneSymbol(gene_symbol.clone()));
    }
    Ok(ids)
}

pub fn run(_common: &crate::common::Args, args: &Args) -> Result<()> {
    let tx_db = load_tx_db(&format!("{}", args.path_db.display()))?
        .tx_db
        .unwrap();
    let tx_db_ids: HashSet<Id> = tx_db
        .gene_to_tx
        .iter()
        .map(|g| Id::from(g.gene_id.clone()))
        .chain(tx_db.transcripts.iter().flat_map(|tx| {
            vec![
                Id::from(tx.id.clone()),
                Id::from(tx.gene_id.clone()),
                Id::from(tx.gene_symbol.clone()),
            ]
        }))
        .collect();
    let cdot_ids = load_cdot_files(&args.path_cdot_json)?;
    let hgnc_ids = load_hgnc_set(&args.path_hgnc_json)?;
    let disease_gene_ids = load_disease_gene_set(&args.path_disease_gene_tsv)?;

    for id in &cdot_ids - &tx_db_ids {
        tracing::warn!("cdot {:?} not found in transcript database", &id);
        if let Some(annotations) = hgnc_ids.get(&id) {
            if !annotations.is_empty() {
                tracing::warn!("HGNC {:?} has annotations {:?}", &id, annotations);
            }
        }
    }

    for id in &disease_gene_ids - &cdot_ids {
        tracing::warn!("disease gene {:?} not found in cdot JSON", &id);
    }

    Ok(())
}
