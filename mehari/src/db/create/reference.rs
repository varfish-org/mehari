use crate::db::create::cli::Args;
use anyhow::{Error, anyhow};
use noodles::core::Region;
use seqrepo::{AliasOrSeqId, Interface, SeqRepo};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::ops::RangeFull;
use std::path::{Path, PathBuf};
use std::time::Instant;

pub struct BasicIndexedFasta {
    fasta_path: PathBuf,
    index: noodles::fasta::fai::Index,
    id_lookup: HashMap<String, String>,
}

impl BasicIndexedFasta {
    pub fn new(fasta_path: PathBuf, index_path: PathBuf) -> Result<Self, Error> {
        let index = File::open(index_path)
            .map(BufReader::new)
            .map(noodles::fasta::fai::io::Reader::new)?
            .read_index()?;

        // Build lookup map once
        let mut id_lookup = HashMap::new();
        for record in index.as_ref() {
            let rec_name = String::from_utf8_lossy(record.name()).to_string();
            let rec_name_base = rec_name.split('.').next().unwrap_or(&rec_name).to_string();

            // Store exact name mapping
            id_lookup.insert(rec_name.clone(), rec_name.clone());
            // Store base name mapping (without version)
            id_lookup.entry(rec_name_base).or_insert(rec_name);
        }

        Ok(Self { fasta_path, index, id_lookup })
    }

    pub fn get_sequence(&self, id: &str) -> Result<Option<String>, Error> {
        let clean_id = id.strip_prefix("transcript:").unwrap_or(id);
        let clean_id_base = clean_id.split('.').next().unwrap_or(clean_id);

        // Use precomputed lookup map
        let target_name = self.id_lookup.get(clean_id)
            .or_else(|| self.id_lookup.get(clean_id_base));

        let target_name = match target_name {
            Some(name) => name.clone(),
            None => {
                tracing::debug!(
                    "Sequence not found in FAI index for clean_id: '{}' (original id: '{}')",
                    clean_id,
                    id
                );
                return Ok(None);
            }
        };

        let mut reader = File::open(&self.fasta_path)
            .map(BufReader::new)
            .map(|r| noodles::fasta::io::IndexedReader::new(r, self.index.clone()))?;

        let region = Region::new::<String, RangeFull>(target_name, RangeFull);
        match reader.query(&region) {
            Ok(record) => {
                let seq = String::from_utf8(record.sequence().as_ref().to_vec())?;
                Ok(Some(seq))
            }
            Err(e) => {
                tracing::warn!("Failed to query region from indexed FASTA: {}", e);
                Ok(None)
            }
        }
    }
}

pub enum SequenceProvider {
    SeqRepo(SeqRepo),
    IndexedFasta(BasicIndexedFasta),
    FastaMap(HashMap<String, String>),
}

impl SequenceProvider {
    pub(crate) fn fetch_sequence(&self, alias: &AliasOrSeqId) -> Result<String, Error> {
        let id = match alias {
            AliasOrSeqId::Alias { value, .. } => value,
            AliasOrSeqId::SeqId(id) => id,
        };

        match self {
            SequenceProvider::SeqRepo(repo) => repo.fetch_sequence(alias).map_err(Into::into),
            SequenceProvider::FastaMap(map) => {
                let clean_id = id.strip_prefix("transcript:").unwrap_or(id);

                map.get(clean_id)
                    .or_else(|| map.get(clean_id.split('.').next().unwrap_or("")))
                    .cloned()
                    .ok_or_else(|| {
                        tracing::debug!(
                            "Sequence not found in FASTA map for clean_id: '{}' (original id: '{}')",
                            clean_id, id
                        );
                        anyhow!("Sequence not found in FASTA map for {}", id)
                    })
            }
            SequenceProvider::IndexedFasta(reader) => {
                if let Some(seq) = reader.get_sequence(id)? {
                    return Ok(seq);
                }

                Err(anyhow!("Sequence not found in indexed FASTA for {}", id))
            }
        }
    }
}

/// Create file-backed `SeqRepo`.
pub fn open_seqrepo(path: impl AsRef<Path>) -> Result<SeqRepo, Error> {
    tracing::info!("Opening seqrepo…");
    let start = Instant::now();
    let seqrepo = PathBuf::from(path.as_ref());
    let p = path.as_ref().to_str();
    let path = seqrepo
        .parent()
        .ok_or(anyhow::anyhow!("Could not get parent from {:?}", &p))?
        .to_str()
        .unwrap()
        .to_string();
    let instance = seqrepo
        .file_name()
        .ok_or(anyhow::anyhow!("Could not get basename from {:?}", &p))?
        .to_str()
        .unwrap()
        .to_string();
    let seqrepo = SeqRepo::new(path, &instance)?;
    tracing::info!("… seqrepo opened in {:?}", start.elapsed());
    Ok(seqrepo)
}

pub fn open_sequence_provider(args: &Args) -> Result<SequenceProvider, Error> {
    if let Some(seqrepo_path) = &args.seqrepo {
        return Ok(SequenceProvider::SeqRepo(open_seqrepo(seqrepo_path)?));
    }

    if let Some(fasta_path) = &args.transcript_sequences {
        let mut fai_path = fasta_path.clone();
        fai_path.as_mut_os_string().push(".fai");

        if fai_path.exists() {
            tracing::info!(
                "Opening indexed FASTA: {} (using index {})",
                fasta_path.display(),
                fai_path.display()
            );
            return Ok(SequenceProvider::IndexedFasta(BasicIndexedFasta::new(
                fasta_path.clone(),
                fai_path,
            )?));
        }

        tracing::info!(
            "Loading non-indexed FASTA into memory: {}",
            fasta_path.display()
        );
        let mut map = HashMap::new();
        let file = File::open(fasta_path)?;

        let is_gz = fasta_path
            .extension()
            .map(|ext| ext == "gz" || ext == "bgz")
            .unwrap_or(false);

        let buf_reader: Box<dyn std::io::BufRead> = if is_gz {
            Box::new(noodles_bgzf::io::Reader::new(file))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut reader = noodles::fasta::io::Reader::new(buf_reader);
        for result in reader.records() {
            let record = result?;
            let raw_name = String::from_utf8_lossy(record.name()).to_string();

            let id = raw_name
                .split_whitespace()
                .next()
                .unwrap_or(&raw_name)
                .to_string();
            let seq = String::from_utf8_lossy(record.sequence().as_ref()).to_string();

            if let Some((base, _version)) = id.split_once('.') {
                map.insert(base.to_string(), seq.clone());
            }

            map.insert(id, seq);
        }
        return Ok(SequenceProvider::FastaMap(map));
    }

    // clap guarantees this is unreachable, but we need to keep the compiler happy
    unreachable!("Either seqrepo or transcript_sequences must be set");
}