use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::iter::repeat;
use std::path::{Path, PathBuf};

use anyhow::Result;
use clap::Parser;
use hgvs::data::cdot::json::models::Container as Cdot;
use itertools::Itertools;
use serde;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use strum::Display;

use crate::annotate::seqvars::load_tx_db;

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

    #[arg(long)]
    pub output: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, PartialOrd, Ord, Display)]
enum Id {
    Hgnc(String),
    GeneSymbol(String),
    GeneName(String),
    NcbiAccession(String),
    NcbiGene(String),
    NcbiTranscript(String),
    EnsemblAccession(String),
    EnsemblGene(String),
    EnsemblTranscript(String),
    Empty,
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

type IdCollection = HashMap<Id, Vec<Id>>;

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
            let hgnc_id = g
                .hgnc
                .as_ref()
                .map_or(Id::Empty, |h| Id::Hgnc(format!("HGNC:{h}")));
            let gid = Id::from(gene_id.clone());
            let mut ids = vec![(hgnc_id.clone(), gid)];

            if let Some(ref gene_symbol) = g.gene_symbol {
                ids.push((hgnc_id, Id::GeneSymbol(gene_symbol.clone())));
            }
            ids
        })
        .chain(transcripts.iter().flat_map(|(tx_id, t)| {
            let hgnc_id = t
                .hgnc
                .as_ref()
                .map_or(Id::Empty, |h| Id::Hgnc(format!("HGNC:{h}")));
            let tid = Id::from(tx_id.clone());
            let mut ids = vec![(hgnc_id.clone(), tid)];
            if let Some(ref gene_symbol) = t.gene_name {
                ids.push((hgnc_id, Id::GeneSymbol(gene_symbol.clone())));
            }
            ids
        }))
        .into_group_map())
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

fn load_hgnc_set(path: impl AsRef<Path>) -> Result<IdCollection> {
    let hgnc: Value = {
        let reader = File::open(path).map(BufReader::new)?;
        serde_json::from_reader::<_, Value>(reader)?["response"]["docs"].to_owned()
    };
    let entries: Vec<HgncEntry> = serde_json::from_value(hgnc)?;
    let result = entries
        .into_iter()
        .flat_map(|entry| {
            let mut ids = vec![];
            let hgnc_id = Id::Hgnc(entry.hgnc_id);
            let ncbi_gene_id = entry.entrez_id.map(Id::NcbiGene);
            let ensembl_gene_id = entry.ensembl_gene_id.map(Id::EnsemblGene);

            ids.push(Id::GeneSymbol(entry.symbol));
            ids.push(Id::GeneName(entry.name));

            if let Some(gene_id) = ncbi_gene_id {
                ids.push(gene_id);
            }

            if let Some(gene_id) = ensembl_gene_id {
                ids.push(gene_id);
            }

            if let Some(refseq_accessions) = entry.refseq_accession {
                ids.extend(refseq_accessions.into_iter().map(Id::NcbiAccession));
            }

            if let Some(mane_select) = entry.mane_select {
                ids.extend(mane_select.into_iter().map(Id::from))
            }
            repeat(hgnc_id).zip(ids)
        })
        .into_group_map();
    Ok(result)
}

fn load_disease_gene_set(path: impl AsRef<Path>) -> Result<IdCollection> {
    let reader = File::open(path).map(BufReader::new)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    let mut ids: IdCollection = HashMap::new();
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
        ids.entry(Id::Hgnc(hgnc_id.clone()))
            .or_default()
            .push(Id::GeneSymbol(gene_symbol.clone()));
    }
    Ok(ids)
}

pub fn run(_common: &crate::common::Args, args: &Args) -> Result<()> {
    let tx_db = load_tx_db(&format!("{}", args.path_db.display()))?
        .tx_db
        .unwrap();
    let filtered_ids: HashSet<Id> = tx_db
        .gene_to_tx
        .iter()
        .flat_map(|g| {
            if let Some(true) = g.filtered {
                Some(Id::Hgnc(g.gene_id.clone()))
            } else {
                None
            }
        })
        .chain(tx_db.transcripts.iter().flat_map(|tx| {
            if let Some(true) = tx.filtered {
                Some(Id::from(tx.id.clone()))
            } else {
                None
            }
        }))
        .collect();

    let tx_db_ids: IdCollection = tx_db
        .gene_to_tx
        .iter()
        .flat_map(|g| {
            let hgnc_id = Id::from(g.gene_id.clone());
            repeat(hgnc_id).zip(g.tx_ids.iter().cloned().map(Id::from))
        })
        .chain(tx_db.transcripts.iter().flat_map(|tx| {
            let hgnc_id = Id::from(tx.gene_id.clone());
            vec![
                (hgnc_id.clone(), Id::from(tx.id.clone())),
                (hgnc_id.clone(), Id::from(tx.gene_symbol.clone())),
            ]
        }))
        .into_group_map();
    let cdot_ids = load_cdot_files(&args.path_cdot_json)?;
    let hgnc_ids = load_hgnc_set(&args.path_hgnc_json)?;
    let disease_gene_ids = load_disease_gene_set(&args.path_disease_gene_tsv)?;

    let cdot_keys: HashSet<_> = cdot_ids.keys().collect();
    let all_tx_db_keys: HashSet<_> = tx_db_ids.keys().collect();
    let disease_gene_keys: HashSet<_> = disease_gene_ids.keys().collect();

    let mut out = File::create(&args.output).map(BufWriter::new)?;
    let mut valid = true;
    for id in &cdot_keys - &all_tx_db_keys {
        let info = hgnc_ids.get(&id);
        let has_transcripts = info.map_or(false, |info| {
            info.iter()
                .any(|a| matches!(a, Id::NcbiTranscript(_) | Id::EnsemblTranscript(_)))
        });
        let is_filtered = filtered_ids.contains(&id);
        if has_transcripts {
            writeln!(&mut out,
                "cdot {id:?} (filtered: {is_filtered}) not found in transcript database, even though the id has associated transcripts in the hgnc complete set: {info:?}",
            )?;
            if !is_filtered {
                valid = false;
            }
        }
    }

    for id in &disease_gene_keys - &cdot_keys {
        let info = hgnc_ids.get(&id);
        writeln!(
            &mut out,
            "disease gene {id:?} not found in cdot, hgnc info: {info:?}"
        )?;
    }

    writeln!(&mut out, "Status:\t{}", if valid { "OK" } else { "ERROR" })?;
    out.flush()?;

    Ok(())
}
