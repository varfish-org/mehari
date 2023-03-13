//! Commonly used code.

use byte_unit::Byte;
use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};

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
