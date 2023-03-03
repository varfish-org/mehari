//! Commonly used code.

use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};

/// Commonly used command line arguments.
#[derive(Parser, Debug)]
pub struct Args {
    /// Verbosity of the program
    #[clap(flatten)]
    pub verbose: Verbosity<InfoLevel>,
}
