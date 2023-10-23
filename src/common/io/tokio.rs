//! Tokio-based async common I/O code.

use async_compression::tokio::bufread::GzipDecoder;
use std::path::Path;
use std::pin::Pin;
use tokio::fs::File;
use tokio::io::{AsyncBufRead, BufReader};

use crate::common::io::std::is_gz;

/// Transparently open a file with gzip decoder.
pub async fn open_read_maybe_gz<P>(path: P) -> Result<Pin<Box<dyn AsyncBufRead>>, anyhow::Error>
where
    P: AsRef<Path>,
{
    tracing::trace!(
        "Opening {} as {} reading",
        path.as_ref().display(),
        "palin text"
    );
    let file = File::open(path.as_ref())
        .await
        .map_err(|e| anyhow::anyhow!("could not open file {}: {}", path.as_ref().display(), e))?;

    if is_gz(path.as_ref()) {
        let bufreader = BufReader::new(file);
        let decoder = {
            let mut decoder = GzipDecoder::new(bufreader);
            decoder.multiple_members(true);
            decoder
        };
        Ok(Box::pin(BufReader::new(decoder)))
    } else {
        Ok(Box::pin(BufReader::new(file)))
    }
}
