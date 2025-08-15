use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use derive_new::new;
use enumflags2::BitFlags;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use serde_with::{serde_as, DisplayFromStr};
use strum::Display;

use crate::annotate::seqvars::load_tx_db;
use crate::db::create::Reason as FilterReason;
use crate::pbs::txs::TranscriptTag;

/// Check a mehari transcript database against information from
/// cdot, HGNC, human-phenotype-ontology and a collection of known issues.
#[derive(Parser, Debug, Clone)]
#[command(about = "Check transcript database", long_about = None)]
pub struct Args {
    /// Path to the transcript database file to check.
    #[arg(long)]
    pub db: PathBuf,

    /// Paths to the cdot JSON files to check against.
    #[arg(long, required = true)]
    pub cdot: Vec<PathBuf>,

    /// Path to the HGNC JSON file to check against.
    #[arg(long)]
    pub hgnc: PathBuf,

    /// Path to the disease gene TSV file to check against.
    #[arg(long)]
    pub disease_genes: PathBuf,

    /// Path to the known issues TSV.
    #[arg(long)]
    pub known_issues: PathBuf,

    /// Path to the output TSV file.
    #[arg(long)]
    pub output: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct KnownIssue {
    id_type: String,
    id: String,
    description: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, PartialOrd, Ord, Display)]
#[serde(tag = "type", content = "id", rename_all = "snake_case")]
enum Id {
    Hgnc(String),
    GeneSymbol(String),
    GeneName(String),
    GeneAlias(String),
    GeneVersion(String),
    NcbiAccession(String),
    NcbiGene(String),
    NcbiTranscript(String),
    EnsemblAccession(String),
    EnsemblGene(String),
    EnsemblTranscript(String),
    Empty,
}

impl Id {
    pub fn to_parts(&self) -> (String, String) {
        match self {
            Id::Hgnc(s) => ("hgnc".into(), s.clone()),
            Id::GeneSymbol(s) => ("gene_symbol".into(), s.clone()),
            Id::GeneName(s) => ("gene_name".into(), s.clone()),
            Id::GeneAlias(s) => ("gene_alias".into(), s.clone()),
            Id::GeneVersion(s) => ("gene_version".into(), s.clone()),
            Id::NcbiAccession(s) => ("ncbi_accession".into(), s.clone()),
            Id::NcbiGene(s) => ("ncbi_gene".into(), s.clone()),
            Id::NcbiTranscript(s) => ("ncbi_transcript".into(), s.clone()),
            Id::EnsemblAccession(s) => ("ensembl_accession".into(), s.clone()),
            Id::EnsemblGene(s) => ("ensembl_gene".into(), s.clone()),
            Id::EnsemblTranscript(s) => ("ensembl_transcript".into(), s.clone()),
            Id::Empty => ("empty".into(), "".into()),
        }
    }
}

impl<T> From<T> for Id
where
    T: AsRef<str>,
{
    fn from(value: T) -> Self {
        let value = value.as_ref().to_string();
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

type IdentifierMap = HashMap<Id, HashSet<Id>>;

#[derive(Debug, new)]
struct ComparisonDataContainer {
    tx_db_data: TxDbData,
    cdot_genes: HashSet<Id>,
    hgnc_map: IdentifierMap,
    disease_genes: HashSet<Id>,
    known_issues: HashMap<Id, String>,
}

impl ComparisonDataContainer {
    pub fn from_args(args: &Args) -> Result<Self> {
        let tx_db_data = TxDbData::from_path(&args.db)?;
        let cdot_map = load_cdot_files(&args.cdot)?;
        let hgnc_map = load_hgnc_set(&args.hgnc)?;
        let disease_gene_map = load_disease_gene_set(&args.disease_genes)?;
        let known_issues = load_known_issues(&args.known_issues)?;

        Ok(ComparisonDataContainer::new(
            tx_db_data,
            cdot_map.into_keys().collect(),
            hgnc_map,
            disease_gene_map.into_keys().collect(),
            known_issues,
        ))
    }
}

#[derive(Debug)]
struct TxDbData {
    /// All gene and transcript IDs present in the database.
    all_ids: HashSet<Id>,
    /// The reason for filtering, for each filtered ID.
    filter_reasons: HashMap<Id, BitFlags<FilterReason>>,
    /// Discarded MANE transcripts and the reason for their exclusion.
    discarded_mane_entries: HashMap<Id, BitFlags<FilterReason>>,
}

impl TxDbData {
    pub fn from_path(path: &Path) -> Result<Self> {
        let tx_db = load_tx_db(format!("{}", path.display()))?.tx_db.unwrap();
        let mut all_ids = HashSet::new();

        for g in &tx_db.gene_to_tx {
            let hgnc_id = Id::Hgnc(g.gene_id.clone());
            all_ids.insert(hgnc_id.clone());

            for tx_id in &g.tx_ids {
                all_ids.insert(Id::from(tx_id));
            }
        }

        let filter_reasons = tx_db
            .transcripts
            .iter()
            .filter(|tx| tx.filtered.unwrap_or(false))
            .map(|tx| {
                (
                    Id::from(&tx.id),
                    BitFlags::<FilterReason>::from_bits(tx.filter_reason.unwrap()).unwrap(),
                )
            })
            .chain(
                tx_db
                    .gene_to_tx
                    .iter()
                    .filter(|g| g.filtered.unwrap_or(false))
                    .map(|g| {
                        (
                            Id::Hgnc(g.gene_id.clone()),
                            BitFlags::<FilterReason>::from_bits(g.filter_reason.unwrap()).unwrap(),
                        )
                    }),
            )
            .collect();

        let discarded_mane_entries = tx_db
            .transcripts
            .iter()
            .filter(|tx| {
                tx.filtered.unwrap_or(false)
                    && (tx.tags.contains(&(TranscriptTag::ManeSelect as i32))
                        || tx.tags.contains(&(TranscriptTag::ManePlusClinical as i32)))
            })
            .map(|tx| {
                (
                    Id::from(&tx.id),
                    BitFlags::<FilterReason>::from_bits(tx.filter_reason.unwrap()).unwrap(),
                )
            })
            .collect();

        Ok(TxDbData {
            all_ids,
            filter_reasons,
            discarded_mane_entries,
        })
    }
}

fn load_cdot_files(paths: &[PathBuf]) -> Result<IdentifierMap> {
    let cdot_container = paths
        .iter()
        .map(crate::db::create::read_cdot_json)
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .reduce(|mut a, b| {
            assert_eq!(a.cdot_version, b.cdot_version);
            assert_eq!(a.genome_builds, b.genome_builds);
            a.genes.extend(b.genes);
            a.transcripts.extend(b.transcripts);
            a
        })
        .context("No cdot files provided")?;

    let mut map = IdentifierMap::new();

    for (gene_id, g) in &cdot_container.genes {
        let hgnc_id = g
            .hgnc
            .as_ref()
            .map_or(Id::Empty, |h| Id::Hgnc(format!("HGNC:{h}")));
        let entry = map.entry(hgnc_id).or_default();
        entry.insert(Id::from(gene_id));
        if let Some(ref gene_symbol) = g.gene_symbol {
            entry.insert(Id::GeneSymbol(gene_symbol.clone()));
        }
        if let Some(ref aliases) = g.aliases {
            entry.extend(aliases.iter().map(|s| Id::GeneSymbol(s.clone())));
        }
    }

    for (tx_id, t) in &cdot_container.transcripts {
        let hgnc_id = t
            .hgnc
            .as_ref()
            .map_or(Id::Empty, |h| Id::Hgnc(format!("HGNC:{h}")));
        let entry = map.entry(hgnc_id).or_default();
        entry.insert(Id::from(tx_id));
        entry.insert(Id::from(&t.id));

        if let Some((accession, _)) = tx_id.rsplit_once('.') {
            if accession.starts_with("ENS") {
                entry.insert(Id::EnsemblAccession(accession.to_string()));
            } else {
                entry.insert(Id::NcbiAccession(accession.to_string()));
            }
        }
        if t.gene_version.starts_with("ENSG") {
            entry.insert(Id::EnsemblGene(t.gene_version.clone()));
        } else {
            entry.insert(Id::NcbiGene(t.gene_version.clone()));
        }
        if let Some(ref gene_name) = t.gene_name {
            entry.insert(Id::GeneName(gene_name.clone()));
        }
    }

    Ok(map)
}

#[derive(Debug, Serialize, Deserialize)]
struct HgncEntry {
    hgnc_id: String,
    location: Option<String>,
    location_sortable: Option<String>,
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

fn load_hgnc_set(path: impl AsRef<Path>) -> Result<IdentifierMap> {
    let reader = File::open(path).map(BufReader::new)?;
    let hgnc_json: Value = serde_json::from_reader(reader)?;
    let entries: Vec<HgncEntry> = serde_json::from_value(hgnc_json["response"]["docs"].clone())?;

    let mut map = IdentifierMap::new();
    for entry in entries {
        let hgnc_id = Id::Hgnc(entry.hgnc_id);
        let ids = map.entry(hgnc_id).or_default();

        ids.insert(Id::GeneSymbol(entry.symbol));
        ids.insert(Id::GeneName(entry.name));

        if let Some(entrez) = entry.entrez_id {
            ids.insert(Id::NcbiGene(entrez));
        }
        if let Some(ensg) = entry.ensembl_gene_id {
            ids.insert(Id::EnsemblGene(ensg));
        }
        if let Some(refseq) = entry.refseq_accession {
            ids.extend(refseq.into_iter().map(Id::NcbiAccession));
        }
        if let Some(mane) = entry.mane_select {
            ids.extend(mane.into_iter().map(|s| Id::from(&s)));
        }
    }
    Ok(map)
}

fn load_disease_gene_set(path: impl AsRef<Path>) -> Result<IdentifierMap> {
    let reader = File::open(path).map(BufReader::new)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);

    let mut map = IdentifierMap::new();
    for result in rdr.deserialize() {
        let record: HashMap<String, String> = result?;
        let gene_symbol = record
            .get("gene_symbol")
            .expect("`gene_symbol` column not found")
            .clone();
        let hgnc_id = record
            .get("hgnc_id")
            .expect("`hgnc_id` column not found")
            .clone();
        map.entry(Id::Hgnc(hgnc_id))
            .or_default()
            .insert(Id::GeneSymbol(gene_symbol));
    }
    Ok(map)
}

fn load_known_issues(path: impl AsRef<Path>) -> Result<HashMap<Id, String>> {
    let reader = File::open(path).map(BufReader::new)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);

    let mut issues = HashMap::new();
    for result in rdr.deserialize() {
        let record: KnownIssue = result?;
        let id = match record.id_type.as_str() {
            "hgnc_id" => Id::Hgnc(record.id),
            "gene_symbol" => Id::GeneSymbol(record.id),
            "refseq_id" => Id::NcbiAccession(record.id),
            "ensembl_id" => Id::EnsemblTranscript(record.id),
            _ => Id::from(&record.id),
        };
        issues.insert(id, record.description);
    }
    Ok(issues)
}

#[derive(Debug, Serialize, Display)]
#[serde(rename_all = "snake_case")]
enum Status {
    Ok,
    Err,
}

#[serde_as]
#[derive(Debug, Serialize, new)]
#[serde(rename_all = "snake_case")]
struct ReportEntry {
    status: Status,
    id_type: String,
    id_value: String,
    check: String,
    is_known_issue: bool,
    #[serde_as(as = "DisplayFromStr")]
    filter_reason: BitFlags<FilterReason>,
    details: Option<String>,
}

struct Reporter {
    entries: Vec<ReportEntry>,
    known_issues: HashMap<Id, String>,
    has_errors: bool,
}

impl Reporter {
    fn new(known_issues: HashMap<Id, String>) -> Self {
        Self {
            entries: Vec::new(),
            known_issues,
            has_errors: false,
        }
    }

    fn add_entry(
        &mut self,
        status: Status,
        id: Id,
        check: String,
        filter_reason: BitFlags<FilterReason>,
        details: Option<String>,
    ) {
        let (is_known_issue, final_status, issue_details) =
            if let Some(desc) = self.known_issues.get(&id) {
                (true, Status::Ok, Some(desc.clone()))
            } else {
                (false, status, details)
            };

        if matches!(final_status, Status::Err) {
            self.has_errors = true;
        }

        let (id_type, id_value) = id.to_parts();

        self.entries.push(ReportEntry::new(
            final_status,
            id_type,
            id_value,
            check,
            is_known_issue,
            filter_reason,
            issue_details,
        ));
    }

    fn write_to_file(mut self, path: &Path) -> Result<()> {
        self.add_entry(
            if self.has_errors {
                Status::Err
            } else {
                Status::Ok
            },
            Id::Empty,
            "Overall Status".to_string(),
            BitFlags::empty(),
            None,
        );

        let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(path)?;
        for entry in self.entries {
            writer.serialize(entry)?;
        }
        writer.flush()?;
        Ok(())
    }
}

pub struct DatabaseChecker {
    args: Args,
}

impl DatabaseChecker {
    pub fn new(args: Args) -> Self {
        Self { args }
    }

    pub fn run(&self) -> Result<()> {
        let data = ComparisonDataContainer::from_args(&self.args)?;
        let mut reporter = Reporter::new(data.known_issues.clone());

        self.check_discarded_mane(&data, &mut reporter);
        self.check_unexpected_filtered_genes(&data, &mut reporter);
        self.check_missing_cdot_entries(&data, &mut reporter);
        self.check_missing_disease_genes(&data, &mut reporter);

        reporter.write_to_file(&self.args.output)?;
        Ok(())
    }

    fn check_discarded_mane(&self, data: &ComparisonDataContainer, reporter: &mut Reporter) {
        for (id, &reason) in &data.tx_db_data.discarded_mane_entries {
            reporter.add_entry(
                Status::Err,
                id.clone(),
                "Discarded MANE transcript".to_string(),
                reason,
                None,
            );
        }
    }

    fn check_unexpected_filtered_genes(
        &self,
        data: &ComparisonDataContainer,
        reporter: &mut Reporter,
    ) {
        let uninteresting_reasons = FilterReason::NoTranscripts
            | FilterReason::PredictedTranscriptsOnly
            | FilterReason::Biotype
            | FilterReason::Pseudogene
            | FilterReason::DeselectedGene
            | FilterReason::UseNmTranscriptInsteadOfNr;

        for (id, &reason) in &data.tx_db_data.filter_reasons {
            if matches!(id, Id::Hgnc(_)) && !reason.intersects(uninteresting_reasons) {
                reporter.add_entry(
                    Status::Err,
                    id.clone(),
                    "Gene filtered for unexpected reason".to_string(),
                    reason,
                    None,
                );
            }
        }
    }

    fn check_missing_cdot_entries(&self, data: &ComparisonDataContainer, reporter: &mut Reporter) {
        for id in data.cdot_genes.difference(&data.tx_db_data.all_ids) {
            let reason = data
                .tx_db_data
                .filter_reasons
                .get(id)
                .cloned()
                .unwrap_or_default();
            if reason.contains(FilterReason::NoTranscripts) {
                continue;
            }

            let info = data.hgnc_map.get(id);

            let has_transcripts = info.is_some_and(|info| {
                info.iter()
                    .any(|a| matches!(a, Id::NcbiTranscript(_) | Id::EnsemblTranscript(_)))
            });

            if !has_transcripts {
                continue;
            }

            let details = info.map(|hgnc_ids| format!("hgnc_ids: {:?}", hgnc_ids));
            reporter.add_entry(
                Status::Err,
                id.clone(),
                "Cdot source missing from tx_db".to_string(),
                reason,
                details,
            );
        }
    }

    fn check_missing_disease_genes(&self, data: &ComparisonDataContainer, reporter: &mut Reporter) {
        for id in data.disease_genes.difference(&data.cdot_genes) {
            let info = data.hgnc_map.get(id);
            let has_transcripts = info.is_some_and(|info| {
                info.iter()
                    .any(|a| matches!(a, Id::NcbiTranscript(_) | Id::EnsemblTranscript(_)))
            });

            let details = if has_transcripts {
                info.map(|hgnc_ids| format!("hgnc_ids: {:?}", hgnc_ids))
            } else {
                Some("No transcripts available in HGNC".to_string())
            };

            reporter.add_entry(
                if has_transcripts {
                    Status::Err
                } else {
                    Status::Ok
                },
                id.clone(),
                "Disease gene missing from cdot".to_string(),
                BitFlags::empty(),
                details,
            );
        }
    }
}

pub fn run(_common: &crate::common::Args, args: &Args) -> Result<()> {
    let checker = DatabaseChecker::new(args.clone());
    checker.run()
}
