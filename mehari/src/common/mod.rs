//! Commonly used code.

use std::ops::Range;

use byte_unit::{Byte, UnitType};
use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};

pub mod contig;
pub mod io;
pub mod noodles;

/// Commonly used command line arguments.
#[derive(Parser, Debug, Default)]
pub struct Args {
    /// Verbosity of the program
    #[clap(flatten)]
    pub verbose: Verbosity<InfoLevel>,
}

#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TsvContigStyle {
    /// Use contig name from input VCF as is.
    Passthrough,
    /// Enforce "chr" prefix, using canonical name from assembly info.
    WithChr,
    /// Enforce no "chr" prefix, using canonical name from assembly info.
    WithoutChr,
    /// Use canonical name, with "chr" prefix for GRCh38 and without for GRCh37.
    #[default]
    Auto,
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
