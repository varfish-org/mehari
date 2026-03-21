//! Parsing of PLINK PED files.
//!
//! PED files are TSV files with the following entries
//!
//! 1. Family name (string)
//! 2. Individual name (string)
//! 3. Paternal ID (string)
//! 4. Maternal ID (string)
//! 5. Sex (integer; M/1=male; F/2=female; other=unknown)
//! 6. Affected (integer; 1=unaffected; 2=affected; other=unknown)

use std::{fmt::Display, path::Path, str::FromStr};

use serde::{Deserialize, Serialize};

/// Encode the sex of an individual in a PED file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Sex {
    Male,
    Female,
    #[default]
    Unknown,
}

impl FromStr for Sex {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "1" | "m" | "M" => Ok(Sex::Male),
            "2" | "f" | "F" => Ok(Sex::Female),
            _ => Ok(Sex::Unknown),
        }
    }
}

impl Display for Sex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Sex::Male => write!(f, "1"),
            Sex::Female => write!(f, "2"),
            Sex::Unknown => write!(f, "0"),
        }
    }
}

/// Encode the disease status of an individual in a PED file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Disease {
    Affected,
    Unaffected,
    #[default]
    Unknown,
}

impl FromStr for Disease {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "1" => Ok(Disease::Unaffected),
            "2" => Ok(Disease::Affected),
            _ => Ok(Disease::Unknown),
        }
    }
}

impl Display for Disease {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Disease::Unaffected => write!(f, "1"),
            Disease::Affected => write!(f, "2"),
            Disease::Unknown => write!(f, "0"),
        }
    }
}

/// Encode an individual in a PED file.
#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct Individual {
    /// Family of the individual.
    pub family: String,
    /// Individual's ID.
    pub name: String,
    /// ID of the individual's father, `"0"` for founder.
    #[serde(with = "string_option")]
    pub father: Option<String>,
    /// ID of the individual's mother, `"0"` for founder.
    #[serde(with = "string_option")]
    pub mother: Option<String>,
    /// Sex of the individual.
    #[serde(with = "string")]
    pub sex: Sex,
    /// Disease of the individual.
    #[serde(with = "string")]
    pub disease: Disease,
}

/// Individuals in a PED file.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Pedigree {
    pub individuals: Vec<Individual>,
}

impl Pedigree {
    /// Load `Pedigree` from the given path.
    pub fn from_path<P>(path: P) -> Result<Self, anyhow::Error>
    where
        P: AsRef<Path>,
    {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .has_headers(false)
            .from_path(path)?;
        let mut individuals = Vec::new();
        for result in rdr.deserialize() {
            let individual: Individual = result?;
            individuals.push(individual);
        }
        Ok(Pedigree { individuals })
    }
}

/// Individuals in a PED file, by name.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PedigreeByName {
    pub individuals: indexmap::IndexMap<String, Individual>,
}

impl PedigreeByName {
    /// Load `Pedigree` from the given path.
    pub fn from_path<P>(path: P) -> Result<Self, anyhow::Error>
    where
        P: AsRef<Path>,
    {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .has_headers(false)
            .from_path(path)?;
        let mut individuals = indexmap::IndexMap::new();
        for result in rdr.deserialize() {
            let individual: Individual = result?;
            individuals.insert(individual.name.clone(), individual);
        }
        Ok(PedigreeByName { individuals })
    }
}

/// Helper for serialize using Display/FromStr.
///
/// cf. https://github.com/serde-rs/serde/issues/1316#issue-332908452
mod string {
    use std::fmt::Display;
    use std::str::FromStr;

    use serde::{de, Deserialize, Deserializer, Serializer};

    pub fn serialize<T, S>(value: &T, serializer: S) -> Result<S::Ok, S::Error>
    where
        T: Display,
        S: Serializer,
    {
        serializer.collect_str(value)
    }

    pub fn deserialize<'de, T, D>(deserializer: D) -> Result<T, D::Error>
    where
        T: FromStr,
        T::Err: Display,
        D: Deserializer<'de>,
    {
        String::deserialize(deserializer)?
            .parse()
            .map_err(de::Error::custom)
    }
}

/// Helper for serializing `Option<String>` where `"0"` encodes `None`.
///
/// cf. https://github.com/serde-rs/serde/issues/1316#issue-332908452
mod string_option {
    use serde::{Deserialize, Deserializer, Serializer};

    pub fn serialize<S>(value: &Option<String>, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match value {
            Some(s) => serializer.collect_str(s),
            None => serializer.collect_str("0"),
        }
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Option<String>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        if s == "0" {
            Ok(None)
        } else {
            Ok(Some(s))
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn sex_display() -> Result<(), anyhow::Error> {
        assert_eq!(format!("{}", Sex::Unknown), "0");
        assert_eq!(format!("{}", Sex::Male), "1");
        assert_eq!(format!("{}", Sex::Female), "2");

        Ok(())
    }

    #[test]
    fn sex_from_str() -> Result<(), anyhow::Error> {
        assert_eq!(Sex::from_str("0")?, Sex::Unknown);
        assert_eq!(Sex::from_str("1")?, Sex::Male);
        assert_eq!(Sex::from_str("2")?, Sex::Female);

        Ok(())
    }

    #[test]
    fn disease_display() -> Result<(), anyhow::Error> {
        assert_eq!(format!("{}", Disease::Unknown), "0");
        assert_eq!(format!("{}", Disease::Unaffected), "1");
        assert_eq!(format!("{}", Disease::Affected), "2");

        Ok(())
    }

    #[test]
    fn disease_from_str() -> Result<(), anyhow::Error> {
        assert_eq!(Disease::from_str("0")?, Disease::Unknown);
        assert_eq!(Disease::from_str("1")?, Disease::Unaffected);
        assert_eq!(Disease::from_str("2")?, Disease::Affected);

        Ok(())
    }

    #[test]
    fn pedigree_from_path() -> Result<(), anyhow::Error> {
        let pedigree = Pedigree::from_path("tests/data/ped/family.ped")?;
        assert_eq!(pedigree.individuals.len(), 3);

        Ok(())
    }

    #[test]
    fn pedigree_by_name_from_path() -> Result<(), anyhow::Error> {
        let pedigree = PedigreeByName::from_path("tests/data/ped/family.ped")?;
        assert_eq!(pedigree.individuals.len(), 3);

        Ok(())
    }
}
