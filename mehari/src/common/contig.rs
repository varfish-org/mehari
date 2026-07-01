//! Contig name harmonization.

use biocommons_bioutils::assemblies::{ASSEMBLY_INFOS, Assembly, Sequence};
use indexmap::IndexMap;

/// A manager for contig name harmonization.
#[derive(Debug, Clone)]
pub struct ContigManager {
    /// Mapping from any known alias (e.g., "1", "chr1", "NC_000001.10") to the RefSeq accession.
    alias_to_accession: IndexMap<String, String>,
    /// Mapping from the RefSeq accession back to the primary sequence info.
    accession_to_info: IndexMap<String, Sequence>,
    /// Mapping from the primary name (e.g., "1", "X", "MT") to the chromosome number.
    name_to_chrom_no: IndexMap<String, u32>,
    /// Name of the assembly.
    assembly: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
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
    /// Create a new manager for a given assembly name.
    pub fn new(assembly_name: &str) -> Self {
        let mut alias_to_accession = IndexMap::new();
        let mut accession_to_info = IndexMap::new();
        let mut name_to_chrom_no = IndexMap::new();

        let assembly = match assembly_name.to_lowercase().as_str() {
            "grch37" | "grch37p10" => Some(Assembly::Grch37p10),
            "grch38" => Some(Assembly::Grch38),
            _ => None,
        };

        if let Some(assembly) = assembly {
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
        } else {
            // Fallback for unknown assemblies: leave assembly-derived alias/accession mappings empty.
            tracing::debug!(
                "Unknown assembly '{}', no assembly-derived contig mappings available",
                assembly_name
            );
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

        let mut additional_aliases = IndexMap::new();
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
            else if let Some(stripped_name) = name.strip_prefix("chr")
                && !alias_to_accession.contains_key(stripped_name)
            {
                additional_aliases.insert(stripped_name.to_string(), accession.clone());
            }
        }

        alias_to_accession.extend(additional_aliases);

        let mt_acc = alias_to_accession
            .get("chrMT")
            .or_else(|| alias_to_accession.get("MT"))
            .cloned();
        if let Some(ref mt_acc) = mt_acc {
            alias_to_accession.insert("M".to_string(), mt_acc.clone());
            alias_to_accession.insert("chrM".to_string(), mt_acc.clone());
        }

        Self {
            assembly: assembly_name.to_string(),
            alias_to_accession,
            accession_to_info,
            name_to_chrom_no,
        }
    }

    pub fn assembly(&self) -> &str {
        &self.assembly
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

    pub fn is_homo_sapiens_canonical(chrom_no: u32) -> bool {
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
    pub fn is_homo_sapiens_canonical_alias(&self, alias: &str) -> bool {
        self.get_chrom_no(alias)
            .is_some_and(Self::is_homo_sapiens_canonical)
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

    pub fn sequences(&self) -> impl Iterator<Item = &Sequence> {
        self.accession_to_info.values()
    }

    /// Canonicalize an alias to be used as a key for annonars' DB accesses.
    /// Basically strips leading "chr" and converts "M" to "MT".
    #[inline]
    pub fn canonicalize(alias: &str) -> String {
        annonars::common::cli::canonicalize(alias)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn get_manager() -> ContigManager {
        ContigManager::new("grch37")
    }

    #[test]
    fn test_canonicalize() {
        assert_eq!(ContigManager::canonicalize("chr1"), "1");
        assert_eq!(ContigManager::canonicalize("1"), "1");
        assert_eq!(ContigManager::canonicalize("chrX"), "X");
        assert_eq!(ContigManager::canonicalize("chrM"), "MT");
        assert_eq!(ContigManager::canonicalize("M"), "MT");
        assert_eq!(ContigManager::canonicalize("MT"), "MT");
    }

    #[test]
    fn test_is_canonical_alias() {
        let cm = get_manager();

        assert!(cm.is_homo_sapiens_canonical_alias("chr1"));
        assert!(cm.is_homo_sapiens_canonical_alias("1"));
        assert!(cm.is_homo_sapiens_canonical_alias("NC_000001.10")); // GRCh37 chr1
        assert!(cm.is_homo_sapiens_canonical_alias("chrX"));
        assert!(cm.is_homo_sapiens_canonical_alias("chrY"));
        assert!(cm.is_homo_sapiens_canonical_alias("chrM"));
        assert!(cm.is_homo_sapiens_canonical_alias("MT"));

        assert!(!cm.is_homo_sapiens_canonical_alias("chrUn_gl000211"));
        assert!(!cm.is_homo_sapiens_canonical_alias("random_string"));
    }

    #[test]
    fn test_is_autosomal_alias() {
        let cm = get_manager();

        assert!(cm.is_autosomal_alias("chr1"));
        assert!(cm.is_autosomal_alias("22"));
        assert!(cm.is_autosomal_alias("NC_000001.10"));

        assert!(!cm.is_autosomal_alias("chrX"));
        assert!(!cm.is_autosomal_alias("MT"));
    }

    #[test]
    fn test_is_gonosomal_alias() {
        let cm = get_manager();

        assert!(cm.is_gonosomal_alias("chrX"));
        assert!(cm.is_gonosomal_alias("Y"));
        assert!(cm.is_gonosomal_alias("NC_000023.10")); // GRCh37 chrX

        assert!(!cm.is_gonosomal_alias("chr1"));
        assert!(!cm.is_gonosomal_alias("MT"));
    }

    #[test]
    fn test_is_mitochondrial_alias() {
        let cm = get_manager();

        assert!(cm.is_mitochondrial_alias("chrM"));
        assert!(cm.is_mitochondrial_alias("MT"));
        assert!(cm.is_mitochondrial_alias("M"));
        assert!(cm.is_mitochondrial_alias("NC_012920.1")); // GRCh37 MT

        assert!(!cm.is_mitochondrial_alias("chr1"));
        assert!(!cm.is_mitochondrial_alias("chrX"));
    }

    #[test]
    fn test_get_chrom_no() {
        let cm = get_manager();

        assert_eq!(cm.get_chrom_no("chr1"), Some(1));
        assert_eq!(cm.get_chrom_no("NC_000001.10"), Some(1));
        assert_eq!(cm.get_chrom_no("chrX"), Some(23));
        assert_eq!(cm.get_chrom_no("Y"), Some(24));
        assert_eq!(cm.get_chrom_no("chrM"), Some(25));
        assert_eq!(cm.get_chrom_no("MT"), Some(25));
        assert_eq!(cm.get_chrom_no("unknown_contig"), None);
    }

    #[test]
    fn test_get_primary_name() {
        let cm = get_manager();

        assert_eq!(cm.get_primary_name("chr1").map(|s| s.as_str()), Some("1"));
        assert_eq!(
            cm.get_primary_name("NC_000001.10").map(|s| s.as_str()),
            Some("1")
        );
        assert_eq!(cm.get_primary_name("chrM").map(|s| s.as_str()), Some("MT"));
        assert_eq!(cm.get_primary_name("M").map(|s| s.as_str()), Some("MT"));
        assert_eq!(cm.get_primary_name("chrX").map(|s| s.as_str()), Some("X"));
    }

    #[test]
    fn test_get_accession() {
        let cm = get_manager();

        assert_eq!(
            cm.get_accession("chr1").map(|s| s.as_str()),
            Some("NC_000001.10")
        );
        assert_eq!(
            cm.get_accession("1").map(|s| s.as_str()),
            Some("NC_000001.10")
        );
        assert_eq!(
            cm.get_accession("chrM").map(|s| s.as_str()),
            Some("NC_012920.1")
        );
    }
}
