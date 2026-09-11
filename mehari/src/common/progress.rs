//! Progress reporting for long-running operations.

use std::fs::File;
use std::path::Path;

pub use indicatif::ProgressBar;
use indicatif::{ProgressBarIter, ProgressDrawTarget};

/// The kind of work unit that a progress bar counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Unit {
    Bytes,
    Transcripts,
}

/// Creates one progress bar per step of a long-running operation.
///
/// The operation advances the bar and finishes it when the step is done. If the operation
/// fails, a bar may stay unfinished.
pub trait Progress: Sync {
    fn bar(&self, step: &str, total: u64, unit: Unit) -> ProgressBar;
}

/// Creates hidden progress bars.
pub struct NoProgress;

impl Progress for NoProgress {
    fn bar(&self, _step: &str, total: u64, _unit: Unit) -> ProgressBar {
        hidden_bar(total)
    }
}

/// A progress bar of length `total` that does not draw.
pub fn hidden_bar(total: u64) -> ProgressBar {
    ProgressBar::with_draw_target(Some(total), ProgressDrawTarget::hidden())
}

/// Opens `path` with a bar "Loading <file name>" that counts the bytes read.
///
/// The caller finishes the bar, which is in the returned reader's `progress` field.
pub fn open_with_progress(
    path: &Path,
    progress: &dyn Progress,
) -> std::io::Result<ProgressBarIter<File>> {
    let file = File::open(path)?;
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    let bar = progress.bar(
        &format!("Loading {name}"),
        file.metadata()?.len(),
        Unit::Bytes,
    );
    Ok(bar.wrap_read(file))
}
