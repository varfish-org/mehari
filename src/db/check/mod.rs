use anyhow::{Context, Result};
use clap::Parser;
use derive_new::new;
use enumflags2::BitFlags;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use serde_with::{serde_as, DisplayFromStr};
use std::collections::{HashMap, HashSet};
use std::fmt::Write;
use std::fs::File;
use std::io::BufReader;
use std::ops::Sub;
use std::path::{Path, PathBuf};
use strum::Display;

use crate::annotate::seqvars::load_tx_db;
use crate::db::create::Reason as FilterReason;
use crate::pbs::txs::{TranscriptDb, TranscriptTag};

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

    /// Path to hgncId → count mapping
    #[arg(long)]
    pub clinvar_hgnc_counts: PathBuf,

    /// Path to txAcc → count mapping
    #[arg(long)]
    pub clinvar_tx_acc_counts: PathBuf,

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

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ClinvarTxAccCount {
    #[serde(rename = "txAcc")]
    tx_acc: String,
    count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ClinvarHgncIdCount {
    #[serde(rename = "hgncId")]
    hgnc_id: String,
    count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, PartialOrd, Ord, Display)]
#[serde(tag = "type", content = "id", rename_all = "snake_case")]
enum Id {
    Hgnc(String),
    GeneSymbol(String),
    GeneName(String),
    GeneAlias(String),
    GeneVersion(String),
    Entrez(String),
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
            Id::Entrez(s) => ("entrez_id".into(), s.clone()),
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
            Id::Entrez(value)
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
    clinvar_hgnc_id_counts: HashMap<Id, usize>,
    clinvar_tx_acc_counts: HashMap<Id, usize>,
}

impl ComparisonDataContainer {
    pub fn from_args(args: &Args) -> Result<Self> {
        let tx_db_data = TxDbData::from_path(&args.db)?;
        let cdot_map = load_cdot_files(&args.cdot)?;
        let hgnc_map = load_hgnc_set(&args.hgnc)?;
        let disease_gene_map = load_disease_gene_set(&args.disease_genes)?;
        let known_issues = load_known_issues(&args.known_issues)?;
        let clinvar_hgnc_id_counts = load_clinvar_hgnc_id_counts(&args.clinvar_hgnc_counts)?;
        let clinvar_tx_acc_counts = load_clinvar_tx_acc_counts(&args.clinvar_tx_acc_counts)?;

        Ok(ComparisonDataContainer::new(
            tx_db_data,
            cdot_map.into_keys().collect(),
            hgnc_map,
            disease_gene_map.into_keys().collect(),
            known_issues,
            clinvar_hgnc_id_counts,
            clinvar_tx_acc_counts,
        ))
    }
}

#[derive(Debug)]
struct TxDbData {
    /// The transcript database.
    tx_db: TranscriptDb,

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
            tx_db,
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
            entry.extend(aliases.iter().map(|s| Id::GeneAlias(s.clone())));
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
            entry.insert(Id::from(accession));
        }
        entry.insert(Id::from(t.gene_version.clone()));
        entry.insert(Id::GeneVersion(t.gene_version.clone()));
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
    alias_name: Option<Vec<String>>,
    alias_symbol: Option<Vec<String>>,
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
            ids.insert(Id::Entrez(entrez));
        }
        if let Some(ensg) = entry.ensembl_gene_id {
            ids.insert(Id::EnsemblGene(ensg));
        }
        if let Some(refseq) = entry.refseq_accession {
            ids.extend(refseq.into_iter().map(Id::from));
        }
        if let Some(mane) = entry.mane_select {
            ids.extend(mane.into_iter().map(Id::from));
        }
        if let Some(alias_symbols) = entry.alias_symbol {
            ids.extend(alias_symbols.into_iter().map(Id::GeneAlias));
        }
        if let Some(alias_names) = entry.alias_name {
            ids.extend(alias_names.into_iter().map(Id::GeneAlias));
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
            "refseq_id" => Id::NcbiTranscript(record.id),
            "ensembl_id" => Id::EnsemblTranscript(record.id),
            _ => Id::from(&record.id),
        };
        issues.insert(id, record.description);
    }
    Ok(issues)
}

fn load_clinvar_hgnc_id_counts(path: impl AsRef<Path>) -> Result<HashMap<Id, usize>> {
    let reader = File::open(path).map(BufReader::new)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    let mut counts = HashMap::new();
    for result in rdr.deserialize() {
        let record: ClinvarHgncIdCount = result?;
        let id = Id::Hgnc(record.hgnc_id);
        *counts.entry(id).or_default() += record.count;
    }
    Ok(counts)
}

fn load_clinvar_tx_acc_counts(path: impl AsRef<Path>) -> Result<HashMap<Id, usize>> {
    let reader = File::open(path).map(BufReader::new)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    let mut counts = HashMap::new();
    for result in rdr.deserialize() {
        let record: ClinvarTxAccCount = result?;
        let id = Id::from(record.tx_acc);
        *counts.entry(id).or_default() += record.count;
    }
    Ok(counts)
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
        self.check_disease_genes_missing_from_cdot(&data, &mut reporter);
        self.check_disease_genes_covered_in_tx_db(&data, &mut reporter);
        self.check_clinvar_hgnc_ids_covered(&data, &mut reporter);
        self.check_clinvar_tx_acc_covered(&data, &mut reporter);

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
            | FilterReason::UseNmTranscriptInsteadOfNr
            | FilterReason::CdsStartOrEndNotConfirmed;

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

    /// Check for disease genes that are missing from the cdot file.
    fn check_disease_genes_missing_from_cdot(
        &self,
        data: &ComparisonDataContainer,
        reporter: &mut Reporter,
    ) {
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

    /// Check that every disease gene is represented by at least one unfiltered transcript in the tx_db.
    fn check_disease_genes_covered_in_tx_db(
        &self,
        data: &ComparisonDataContainer,
        reporter: &mut Reporter,
    ) {
        for disease_gene_id in &data.disease_genes {
            if !matches!(disease_gene_id, Id::Hgnc(_)) {
                panic!(
                    "disease_gene_id should be an hgnc id, but is {}",
                    &disease_gene_id
                );
            }

            let gene_in_db = data
                .tx_db_data
                .tx_db
                .gene_to_tx
                .iter()
                .find(|g| g.gene_id == disease_gene_id.to_parts().1);

            let (is_covered, details) = if let Some(gene_entry) = gene_in_db {
                // The gene exists in the database. Now check its transcripts.
                let db_tx_ids = &gene_entry.tx_ids;

                let has_unfiltered_tx = db_tx_ids.iter().any(|tx_id_str| {
                    let tx_id = Id::from(tx_id_str.as_str());
                    !data.tx_db_data.filter_reasons.contains_key(&tx_id)
                });

                if has_unfiltered_tx {
                    (true, None)
                } else {
                    let mut details_str =
                        "Gene is in db, but has no unfiltered transcripts. ".to_string();
                    if db_tx_ids.is_empty() {
                        write!(
                            details_str,
                            "No transcripts are associated with this gene in the db."
                        )
                        .unwrap();
                    } else {
                        let filtered_txs: Vec<String> = db_tx_ids
                            .iter()
                            .filter_map(|tx_id_str| {
                                let tx_id = Id::from(tx_id_str.as_str());
                                data.tx_db_data
                                    .filter_reasons
                                    .get(&tx_id)
                                    .map(|reason| format!("{} (filtered: {})", tx_id, reason))
                            })
                            .collect();
                        write!(
                            details_str,
                            "All associated transcripts were filtered: [{}]",
                            filtered_txs.join(", ")
                        )
                        .unwrap();
                    }
                    (false, Some(details_str))
                }
            } else {
                (
                    false,
                    Some(
                        "Disease gene HGNC ID not found in the transcript database's gene map."
                            .to_string(),
                    ),
                )
            };

            if !is_covered {
                reporter.add_entry(
                    Status::Err,
                    disease_gene_id.clone(),
                    "Disease gene not covered by any transcript in tx_db".to_string(),
                    BitFlags::empty(),
                    details,
                );
            }
        }
    }

    fn check_clinvar_hgnc_ids_covered(
        &self,
        data: &ComparisonDataContainer,
        reporter: &mut Reporter,
    ) {
        let hgnc_ids_in_clinvar: HashSet<Id> =
            data.clinvar_hgnc_id_counts.keys().cloned().collect();
        let hgnc_ids_not_covered = hgnc_ids_in_clinvar.sub(&data.tx_db_data.all_ids);
        for hgnc_id in hgnc_ids_not_covered {
            let count = data
                .clinvar_hgnc_id_counts
                .get(&hgnc_id)
                .expect("HgncId guaranteed to exist.");
            reporter.add_entry(
                Status::Err,
                hgnc_id.clone(),
                "Clinvar HgncId not contained in tx_db".to_string(),
                BitFlags::empty(),
                Some(format!("Count in clinvar: {count}")),
            );
        }
    }

    fn check_clinvar_tx_acc_covered(
        &self,
        data: &ComparisonDataContainer,
        reporter: &mut Reporter,
    ) {
        let db_base_accessions: HashSet<&str> = data
            .tx_db_data
            .all_ids
            .iter()
            .filter_map(|id| {
                let s = match id {
                    Id::NcbiTranscript(s) => s.as_str(),
                    Id::EnsemblTranscript(s) => s.as_str(),
                    _ => return None,
                };
                s.rsplit_once('.').map(|(base, _)| base).or(Some(s))
            })
            .collect();

        for (tx_acc, count) in &data.clinvar_tx_acc_counts {
            if data.tx_db_data.all_ids.contains(tx_acc) {
                continue;
            }

            let (_, id_value) = tx_acc.to_parts();
            let query_base = id_value
                .rsplit_once('.')
                .map(|(base, _)| base)
                .unwrap_or(&id_value);

            if db_base_accessions.contains(query_base) {
                continue;
            }

            reporter.add_entry(
                Status::Err,
                tx_acc.clone(),
                "Clinvar TxAcc not contained in tx_db".to_string(),
                BitFlags::empty(),
                Some(format!("Count in clinvar: {count}")),
            );
        }
    }
}

pub fn run(_common: &crate::common::Args, args: &Args) -> Result<()> {
    let checker = DatabaseChecker::new(args.clone());
    checker.run()
}
