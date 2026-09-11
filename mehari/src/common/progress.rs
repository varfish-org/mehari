//! Progress reporting for long-running operations.

use std::fs::File;
use std::io::Read;
use std::path::Path;

/// The kind of work unit that a progress step counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Unit {
    Bytes,
    Transcripts,
}

/// Receives progress reports from a long-running operation.
///
/// The operation reports one step at a time: `start`, then any number of `advance` calls,
/// then `finish`. `advance` may be called from several threads at once. If the operation
/// fails, it may return without calling `finish`.
pub trait Progress: Sync {
    /// Starts a step with `total` units of work.
    fn start(&self, step: &str, total: u64, unit: Unit);

    /// Reports `n` more units of the current step as done.
    fn advance(&self, n: u64);

    /// Finishes the current step.
    fn finish(&self);
}

/// Ignores all progress reports.
pub struct NoProgress;

impl Progress for NoProgress {
    fn start(&self, _step: &str, _total: u64, _unit: Unit) {}
    fn advance(&self, _n: u64) {}
    fn finish(&self) {}
}

/// Reports the bytes read through it as progress.
pub struct ProgressReader<'a, R> {
    inner: R,
    progress: &'a dyn Progress,
}

impl<R: Read> Read for ProgressReader<'_, R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.progress.advance(n as u64);
        Ok(n)
    }
}

/// Opens `path` and starts a step "Loading <file name>" over the file size.
///
/// Reading the file to its end advances the step to its total. The caller finishes the step.
pub fn open_with_progress<'a>(
    path: &Path,
    progress: &'a dyn Progress,
) -> std::io::Result<ProgressReader<'a, File>> {
    let file = File::open(path)?;
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    progress.start(
        &format!("Loading {name}"),
        file.metadata()?.len(),
        Unit::Bytes,
    );
    Ok(ProgressReader {
        inner: file,
        progress,
    })
}
