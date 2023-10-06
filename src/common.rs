//! Commonly used code.

use std::{
    fs::File,
    io::{BufRead, BufReader},
    ops::Range,
    path::Path,
};

use byte_unit::Byte;
use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};
use flate2::bufread::MultiGzDecoder;
use hgvs::static_data::Assembly;

/// Commonly used command line arguments.
#[derive(Parser, Debug)]
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
        Byte::from_bytes((me.stat().unwrap().rss * page_size) as u128).get_appropriate_unit(true)
    );
}

/// Select the genome release to use.
#[derive(
    clap::ValueEnum, serde::Serialize, serde::Deserialize, Clone, Copy, Debug, PartialEq, Eq, Hash,
)]
#[serde(rename_all = "snake_case")]
pub enum GenomeRelease {
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

/// Transparently open a file with gzip decoder.
pub fn open_read_maybe_gz<P>(path: P) -> Result<Box<dyn BufRead>, anyhow::Error>
where
    P: AsRef<Path>,
{
    if path.as_ref().extension().map(|s| s.to_str()) == Some(Some("gz")) {
        tracing::trace!("Opening {:?} as gzip for reading", path.as_ref());
        let file = File::open(path)?;
        let bufreader = BufReader::new(file);
        let decoder = MultiGzDecoder::new(bufreader);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        tracing::trace!("Opening {:?} as plain text for reading", path.as_ref());
        let file = File::open(path).map(BufReader::new)?;
        Ok(Box::new(BufReader::new(file)))
    }
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
    ($($expr:expr),*) => {
        let mut settings = insta::Settings::clone_current();
        settings.set_snapshot_suffix(format!($($expr,)*));
        let _guard = settings.bind_to_scope();
    }
}

pub use set_snapshot_suffix;

#[cfg(test)]
mod test {
    #[rstest::rstest]
    #[case(true)]
    #[case(false)]
    fn open_read_maybe_gz(#[case] is_gzip: bool) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{:?}", is_gzip);

        let mut f = super::open_read_maybe_gz(if is_gzip {
            "tests/common/test.txt.gz"
        } else {
            "tests/common/test.txt"
        })?;

        let mut buf = String::new();
        f.read_to_string(&mut buf)?;

        insta::assert_snapshot!(&buf);

        Ok(())
    }
}
