use crate::annotate::seqvars::reference::{ReferenceReader, UnbufferedIndexedFastaAccess};
use crate::common::contig::ContigManager;
use crate::db::create::cli::Args;
use anyhow::{Error, anyhow};
use seqrepo::{AliasOrSeqId, Interface, SeqRepo};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::time::Instant;

pub enum SequenceProvider {
    SeqRepo(SeqRepo),
    IndexedFasta(UnbufferedIndexedFastaAccess),
    FastaMap(HashMap<String, String>),
}

impl SequenceProvider {
    pub(crate) fn fetch_sequence(&self, alias: &AliasOrSeqId) -> Result<String, Error> {
        let id = match alias {
            AliasOrSeqId::Alias { value, .. } => value,
            AliasOrSeqId::SeqId(id) => id,
        };

        let clean_id = id.strip_prefix("transcript:").unwrap_or(id);
        let versionless_clean_id = clean_id.split('.').next().unwrap_or("");

        match self {
            SequenceProvider::SeqRepo(repo) => repo.fetch_sequence(alias).map_err(Into::into),
            SequenceProvider::FastaMap(map) => map
                .get(clean_id)
                .or_else(|| map.get(versionless_clean_id))
                .cloned()
                .ok_or_else(|| {
                    tracing::debug!(
                        "Sequence not found in FASTA map for clean_id: '{}' (original id: '{}')",
                        clean_id,
                        id
                    );
                    anyhow!("Sequence not found in FASTA map for {}", id)
                }),
            SequenceProvider::IndexedFasta(reader) => {
                if let Some(seq) = reader.get(id, None, None)? {
                    return Ok(String::from_utf8(seq)?);
                }

                if clean_id != id
                    && let Some(seq) = reader.get(clean_id, None, None)?
                {
                    return Ok(String::from_utf8(seq)?);
                }
                if versionless_clean_id != clean_id
                    && versionless_clean_id != id
                    && let Some(seq) = reader.get(versionless_clean_id, None, None)?
                {
                    return Ok(String::from_utf8(seq)?);
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
    let p = path.as_ref();
    let path = seqrepo
        .parent()
        .ok_or(anyhow::anyhow!("Could not get parent from {:?}", &p))?
        .to_str()
        .ok_or_else(|| anyhow!("Could not convert parent path to UTF-8 from {:?}", &p))?
        .to_string();
    let instance = seqrepo
        .file_name()
        .ok_or(anyhow::anyhow!("Could not get basename from {:?}", &p))?
        .to_str()
        .ok_or_else(|| anyhow!("Could not convert basename to UTF-8 from {:?}", &p))?
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
            return Ok(SequenceProvider::IndexedFasta(
                UnbufferedIndexedFastaAccess::from_path(
                    fasta_path.clone(),
                    std::sync::Arc::new(ContigManager::new(&args.assembly)),
                )?,
            ));
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
