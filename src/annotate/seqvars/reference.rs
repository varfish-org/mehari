use anyhow::{Context, Error};
use biocommons_bioutils::assemblies::{Assembly, ASSEMBLY_INFOS};
use noodles::fasta::repository::Repository;
use std::collections::HashSet;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Contig {
    name: String,
    length: u64,
}

impl Contig {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn length(&self) -> u64 {
        self.length
    }
}

#[derive(Clone)]
pub struct Reference {
    repository: Repository,
    contigs: Vec<Contig>,
}

impl Reference {
    pub fn contigs(&self) -> &[Contig] {
        &self.contigs
    }

    pub fn repository(&self) -> &Repository {
        &self.repository
    }

    pub fn guess_assembly(&self) -> Option<Assembly> {
        let these_contigs = self
            .contigs
            .iter()
            .map(|c| (c.name().to_string(), c.length()))
            .collect::<HashSet<_>>();
        for (assembly, info) in ASSEMBLY_INFOS.iter() {
            let assembly_contigs = info
                .sequences
                .iter()
                .flat_map(|s| {
                    s.aliases
                        .iter()
                        .cloned()
                        .chain([s.name.clone(), s.refseq_ac.clone()])
                        .map(|name| (name, s.length as u64))
                })
                .collect::<HashSet<_>>();
            if assembly_contigs.intersection(&these_contigs).count() == these_contigs.len() {
                return Some(assembly);
            }
        }
        None
    }

    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, Error> {
        genome_reference(path)
    }
}

pub fn genome_reference(path: impl AsRef<Path>) -> Result<Reference, Error> {
    tracing::info!("Loading reference from {}", path.as_ref().display());
    let reader = noodles::fasta::io::indexed_reader::Builder::default()
        .build_from_path(path.as_ref().canonicalize()?)?;
    let contigs = reader
        .index()
        .iter()
        .map(|record| Contig {
            name: std::str::from_utf8(record.name())
                .with_context(|| format!("Failed parsing contig name '{:?}'", record.name()))
                .unwrap()
                .to_string(),
            length: record.length(),
        })
        .collect();
    let adapter = noodles::fasta::repository::adapters::IndexedReader::new(reader);
    Ok(Reference {
        repository: Repository::new(adapter),
        contigs,
    })
}
