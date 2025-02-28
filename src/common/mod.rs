//! Commonly used code.

use annonars::common::cli::CANONICAL;
use annonars::freqs::cli::import::reading::ContigMap;
use std::collections::HashMap;
use std::ops::Range;

use crate::pbs::txs::GenomeBuild;
use biocommons_bioutils::assemblies::{Assembly, ASSEMBLY_INFOS};
use byte_unit::{Byte, UnitType};
use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};

pub mod io;
pub mod noodles;

/// Commonly used command line arguments.
#[derive(Parser, Debug, Default)]
pub struct Args {
    /// Verbosity of the program
    #[clap(flatten)]
    pub verbose: Verbosity<InfoLevel>,
}

/// Helper to print the current memory resident set size via `tracing`.
pub fn trace_rss_now() {
    let me = procfs::process::Process::myself().unwrap();
    let page_size = procfs::page_size();
    tracing::debug!(
        "RSS now: {}",
        Byte::from_u128((me.stat().unwrap().rss * page_size) as u128)
            .expect("RSS memory computation failed")
            .get_appropriate_unit(UnitType::Binary)
    );
}

/// Select the genome release to use.
#[derive(
    clap::ValueEnum,
    serde::Serialize,
    serde::Deserialize,
    Clone,
    Copy,
    Debug,
    PartialEq,
    Eq,
    Hash,
    Default,
    utoipa::ToSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum GenomeRelease {
    #[default]
    Grch37,
    Grch38,
}

impl GenomeRelease {
    pub fn name(&self) -> String {
        match self {
            GenomeRelease::Grch37 => String::from("GRCh37"),
            GenomeRelease::Grch38 => String::from("GRCh38"),
        }
    }
}

impl From<GenomeRelease> for Assembly {
    fn from(val: GenomeRelease) -> Self {
        match val {
            GenomeRelease::Grch37 => Assembly::Grch37p10,
            GenomeRelease::Grch38 => Assembly::Grch38,
        }
    }
}

impl From<Assembly> for GenomeRelease {
    fn from(assembly: Assembly) -> Self {
        match assembly {
            Assembly::Grch37 | Assembly::Grch37p10 => GenomeRelease::Grch37,
            Assembly::Grch38 => GenomeRelease::Grch38,
        }
    }
}

impl From<GenomeRelease> for GenomeBuild {
    fn from(val: GenomeRelease) -> Self {
        match val {
            GenomeRelease::Grch37 => GenomeBuild::Grch37,
            GenomeRelease::Grch38 => GenomeBuild::Grch38,
        }
    }
}

impl TryFrom<GenomeBuild> for GenomeRelease {
    type Error = anyhow::Error;

    fn try_from(value: GenomeBuild) -> Result<Self, Self::Error> {
        match value {
            GenomeBuild::Grch37 => Ok(GenomeRelease::Grch37),
            GenomeBuild::Grch38 => Ok(GenomeRelease::Grch38),
            _ => anyhow::bail!("Unknown genome build"),
        }
    }
}

// Compute reciprocal overlap between two ranges.
pub fn reciprocal_overlap(lhs: Range<i32>, rhs: Range<i32>) -> f32 {
    let lhs_b = lhs.start;
    let lhs_e = lhs.end;
    let rhs_b = rhs.start;
    let rhs_e = rhs.end;
    let ovl_b = std::cmp::max(lhs_b, rhs_b);
    let ovl_e = std::cmp::min(lhs_e, rhs_e);
    if ovl_b >= ovl_e {
        0f32
    } else {
        let ovl_len = (ovl_e - ovl_b) as f32;
        let x1 = ovl_len / (lhs_e - lhs_b) as f32;
        let x2 = ovl_len / (rhs_e - rhs_b) as f32;
        x1.min(x2)
    }
}

/// The version of `mehari` package.
#[cfg(not(test))]
const VERSION: &str = env!("CARGO_PKG_VERSION");

/// This allows us to override the version to `0.0.0` in tests.
pub fn version() -> &'static str {
    #[cfg(test)]
    return "0.0.0";
    #[cfg(not(test))]
    return VERSION;
}

/// Version information that is returned by the HTTP server.
#[derive(serde::Serialize, serde::Deserialize, Default, Debug, Clone)]
#[serde_with::skip_serializing_none]
#[serde(rename_all = "snake_case")]
pub struct Version {
    /// Version of the transcript database data.
    pub tx_db: Option<String>,
    /// Version of the `mehari` package.
    pub mehari: String,
}

impl Version {
    /// Construct a new version.
    ///
    /// The mehari version is filled automatically.
    pub fn new(tx_db: Option<String>) -> Self {
        Self {
            tx_db,
            mehari: version().to_string(),
        }
    }
}

#[macro_export]
macro_rules! set_snapshot_suffix {
    ($($expr:expr_2021),*) => {
        let mut settings = insta::Settings::clone_current();
        settings.set_snapshot_suffix(format!($($expr,)*));
        let _guard = settings.bind_to_scope();
    }
}

pub use set_snapshot_suffix;

/// Guess the assembly from the given header.
///
/// If the header only contains chrM, for example, the result may be ambiguous. Use `ambiguous_ok`
/// to allow or disallow this.  You can specify an initial value for the assembly to overcome
/// issues.  If the result is incompatible with the `initial_assembly` then an error will
/// be returned.
pub fn guess_assembly(
    vcf_header: &::noodles::vcf::Header,
    ambiguous_ok: bool,
    initial_assembly: Option<Assembly>,
) -> Result<Assembly, anyhow::Error> {
    let mut result = initial_assembly;

    let assembly_infos = [
        (Assembly::Grch37p10, &ASSEMBLY_INFOS[Assembly::Grch37p10]),
        (Assembly::Grch38, &ASSEMBLY_INFOS[Assembly::Grch38]),
    ];

    // Check each assembly.
    for (assembly, info) in assembly_infos.iter() {
        // Collect contig name / length pairs for the assembly.
        let contig_map = ContigMap::new(*assembly);
        let mut lengths = HashMap::new();
        for seq in &info.sequences {
            if CANONICAL.contains(&seq.name.as_str()) {
                lengths.insert(
                    contig_map.name_map.get(seq.name.as_str()).unwrap(),
                    seq.length,
                );
            }
        }

        // Count compatible and incompatible contigs.
        let mut incompatible = 0;
        let mut compatible = 0;
        for (name, data) in vcf_header.contigs() {
            if let Some(length) = data.length() {
                let idx = contig_map.name_map.get(name);
                if let Some(idx) = idx {
                    let name = &info.sequences[*idx].name;
                    if CANONICAL.contains(&name.as_ref()) {
                        if *lengths.get(idx).unwrap() == length {
                            compatible += 1;
                        } else {
                            incompatible += 1;
                        }
                    }
                }
            } else {
                tracing::warn!(
                    "Cannot guess assembly because no length for contig {}",
                    &name
                );
                compatible = 0;
                break;
            }
        }

        if compatible > 0 && incompatible == 0 {
            // Found a compatible assembly.  Check if we already have one and bail out if
            // ambiguity is not allowed.  Anyway, we only keep the first found compatible
            // assembly.
            if let Some(result) = result {
                if result != *assembly && !ambiguous_ok {
                    return Err(anyhow::anyhow!(
                        "Found ambiguity;  initial={:?}, previous={:?}, current={:?}",
                        initial_assembly,
                        result,
                        assembly,
                    ));
                }
                // else: do not re-assign
            } else {
                result = Some(*assembly);
            }
        } else {
            // Found incompatible assembly, bail out if is the initial assembly.
            if let Some(initial_assembly) = initial_assembly {
                if initial_assembly == *assembly {
                    return Err(anyhow::anyhow!(
                        "Incompatible with initial assembly {:?}",
                        result.unwrap()
                    ));
                }
            }
        }
    }

    if let Some(result) = result {
        Ok(result)
    } else {
        Err(anyhow::anyhow!("No matching assembly found"))
    }
}
