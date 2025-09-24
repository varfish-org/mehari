//! Contig name harmonization.

use biocommons_bioutils::assemblies::{Assembly, Sequence, ASSEMBLY_INFOS};
use std::collections::HashMap;

/// A manager for contig name harmonization.
#[derive(Debug, Clone)]
pub struct ContigManager {
    /// Mapping from any known alias (e.g., "1", "chr1", "NC_000001.10") to the RefSeq accession.
    alias_to_accession: HashMap<String, String>,
    /// Mapping from the RefSeq accession back to the primary sequence info.
    accession_to_info: HashMap<String, Sequence>,
    /// Mapping from the primary name (e.g., "1", "X", "MT") to the chromosome number.
    name_to_chrom_no: HashMap<String, u32>,
}

pub struct ContigInfo {
    /// The name without a "chr" prefix (e.g., "1", "X", "MT").
    pub name_without_chr: String,

    /// The name with a "chr" prefix (e.g., "chr1", "chrX", "chrMT").
    pub name_with_chr: String,

    /// The RefSeq accession.
    pub accession: String,

    /// The chromosome number (1-25) only.
    pub chrom_no: u32,
}

const CHR_X: u32 = 23;
const CHR_Y: u32 = 24;
const CHR_M: u32 = 25;

impl ContigManager {
    /// Create a new manager for a given assembly.
    pub fn new(assembly: Assembly) -> Self {
        let mut alias_to_accession = HashMap::new();
        let mut accession_to_info = HashMap::new();
        let mut name_to_chrom_no = HashMap::new();

        for seq in &ASSEMBLY_INFOS[assembly].sequences {
            // Skip non-primary sequences, but keep chrMT.
            if !["Primary Assembly", "non-nuclear"].contains(&&*seq.assembly_unit)
                || seq.sequence_role != "assembled-molecule"
            {
                tracing::debug!("Skipping non-primary sequence: {:?}", seq);
                continue;
            }
            // Store mapping from accession to the full sequence info.
            accession_to_info.insert(seq.refseq_ac.clone(), seq.clone());

            // Map all known identifiers to the RefSeq accession.
            alias_to_accession.insert(seq.name.clone(), seq.refseq_ac.clone());
            alias_to_accession.insert(seq.refseq_ac.clone(), seq.refseq_ac.clone());
            for alias in &seq.aliases {
                alias_to_accession.insert(alias.clone(), seq.refseq_ac.clone());
            }
        }

        // Build chrom_no map based on the primary sequence names.
        for i in 1..=22 {
            name_to_chrom_no.insert(format!("{}", i), i);
            name_to_chrom_no.insert(format!("chr{}", i), i);
        }
        name_to_chrom_no.insert("X".to_string(), CHR_X);
        name_to_chrom_no.insert("chrX".to_string(), CHR_X);
        name_to_chrom_no.insert("Y".to_string(), CHR_Y);
        name_to_chrom_no.insert("chrY".to_string(), CHR_Y);
        name_to_chrom_no.insert("MT".to_string(), CHR_M);
        name_to_chrom_no.insert("chrMT".to_string(), CHR_M);

        // Add "M" and "chrM" as aliases for mitochondrial
        name_to_chrom_no.insert("M".to_string(), CHR_M);
        name_to_chrom_no.insert("chrM".to_string(), CHR_M);

        let mut additional_aliases = HashMap::new();
        for (accession, info) in &accession_to_info {
            let name = &info.name;

            // Case 1: "1", "X", etc. Add "chr1", "chrX" as an alias.
            if !name.starts_with("chr") {
                let chr_name = format!("chr{}", name);
                if !alias_to_accession.contains_key(&chr_name) {
                    additional_aliases.insert(chr_name, accession.clone());
                }
            }
            // Case 2: "chr1", "chrX". Add "1", "X" as an alias.
            else if let Some(stripped_name) = name.strip_prefix("chr") {
                if !alias_to_accession.contains_key(stripped_name) {
                    additional_aliases.insert(stripped_name.to_string(), accession.clone());
                }
            }
        }

        alias_to_accession.extend(additional_aliases);

        let mt_acc = alias_to_accession.get("chrMT").cloned();
        if let Some(ref mt_acc) = mt_acc {
            alias_to_accession.insert("M".to_string(), mt_acc.clone());
            alias_to_accession.insert("chrM".to_string(), mt_acc.clone());
        }

        Self {
            alias_to_accession,
            accession_to_info,
            name_to_chrom_no,
        }
    }

    /// Get the RefSeq accession for any given contig name/alias.
    #[inline]
    pub fn get_accession(&self, alias: &str) -> Option<&String> {
        self.alias_to_accession.get(alias)
    }

    /// Get the primary display name (e.g., "1", "X", "MT") for any given alias.
    #[inline]
    pub(crate) fn get_primary_name(&self, alias: &str) -> Option<&String> {
        self.get_accession(alias)
            .and_then(|ac| self.accession_to_info.get(ac))
            .map(|info| &info.name)
    }

    /// Get the chromosome number (1-22, 23=X, 24=Y, 25=MT) for any given alias.
    #[inline]
    pub fn get_chrom_no(&self, alias: &str) -> Option<u32> {
        self.get_primary_name(alias)
            .and_then(|name| self.name_to_chrom_no.get(name).copied())
            .or_else(|| self.name_to_chrom_no.get(alias).copied())
    }

    /// Check if the contig is an autosome (chr1-22).
    #[inline]
    pub fn is_autosomal(chrom_no: u32) -> bool {
        (1..=22).contains(&chrom_no)
    }

    /// Check if the contig is a gonosome (chrX or chrY).
    #[inline]
    pub fn is_gonosomal(chrom_no: u32) -> bool {
        chrom_no == CHR_X || chrom_no == CHR_Y
    }

    /// Check if the contig is chromosome X.
    #[inline]
    pub fn is_chr_x(chrom_no: u32) -> bool {
        chrom_no == CHR_X
    }

    /// Check if the contig is chromosome Y.
    #[inline]
    pub fn is_chr_y(chrom_no: u32) -> bool {
        chrom_no == CHR_Y
    }

    /// Check if the contig is mitochondrial DNA.
    #[inline]
    pub fn is_mitochondrial(chrom_no: u32) -> bool {
        chrom_no == CHR_M
    }

    pub fn is_canonical(chrom_no: u32) -> bool {
        (1..=25).contains(&chrom_no)
    }

    /// Check if the contig is an autosome (chr1-22).
    #[inline]
    pub fn is_autosomal_alias(&self, alias: &str) -> bool {
        self.get_chrom_no(alias).is_some_and(Self::is_autosomal)
    }

    /// Check if the contig is a gonosome (chrX or chrY).
    #[inline]
    pub fn is_gonosomal_alias(&self, alias: &str) -> bool {
        self.get_chrom_no(alias).is_some_and(Self::is_gonosomal)
    }

    /// Check if the contig is chromosome X.
    #[inline]
    pub fn is_chr_x_alias(&self, alias: &str) -> bool {
        self.get_chrom_no(alias).is_some_and(Self::is_chr_x)
    }

    /// Check if the contig is chromosome Y.
    #[inline]
    pub fn is_chr_y_alias(&self, alias: &str) -> bool {
        self.get_chrom_no(alias).is_some_and(Self::is_chr_y)
    }

    /// Check if the contig is mitochondrial DNA.
    #[inline]
    pub fn is_mitochondrial_alias(&self, alias: &str) -> bool {
        self.get_chrom_no(alias).is_some_and(Self::is_mitochondrial)
    }

    /// Check if the contig is a canonical chromosome (chr1-22, chrX, chrY, or chrMT).
    #[inline]
    pub fn is_canonical_alias(&self, alias: &str) -> bool {
        self.get_chrom_no(alias).is_some_and(Self::is_canonical)
    }

    pub fn get_contig_info(&self, alias: &str) -> Option<ContigInfo> {
        let accession = self.get_accession(alias)?;
        let seq_info = self.accession_to_info.get(accession)?;
        let chrom_no = self.get_chrom_no(alias)?;

        let primary_name = &seq_info.name;
        let (name_with_chr, name_without_chr) = if primary_name.starts_with("chr") {
            (
                primary_name.clone(),
                primary_name.strip_prefix("chr").unwrap().to_string(),
            )
        } else {
            (format!("chr{}", primary_name), primary_name.clone())
        };

        Some(ContigInfo {
            name_without_chr,
            name_with_chr,
            accession: accession.clone(),
            chrom_no,
        })
    }
}
