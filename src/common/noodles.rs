//! Helper code for using noodles.

use crate::common::io::tokio::open_read_maybe_bgzf;
use futures::future::{join_all, BoxFuture};
use futures::Stream;
use noodles::vcf::variant::{Record, RecordBuf};
use noodles::vcf::Header;
use std::path::Path;
use tokio::io;
use tokio::io::{AsyncBufRead, AsyncWrite};

use super::io::{tokio::open_read_maybe_gz, tokio::open_write_maybe_bgzf};

/// Alias for the async vcf reader type that we will use.
pub type AsyncVcfReader = noodles::vcf::AsyncReader<std::pin::Pin<Box<dyn AsyncBufRead>>>;

/// A variant format reader.
pub(crate) trait AsyncNoodlesReader<R> {
    /// Reads a VCF header.
    async fn read_vcf_header(&mut self) -> io::Result<Header>;

    /// Returns an iterator over records.
    async fn vcf_records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h Header,
    ) -> impl Stream<Item = io::Result<impl Record>> + 'r;
}

impl<R> AsyncNoodlesReader<R> for AsyncVcfReader {
    async fn read_vcf_header(&mut self) -> io::Result<Header> {
        self.read_header().await
    }

    async fn vcf_records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h Header,
    ) -> impl Stream<Item = io::Result<impl noodles::vcf::variant::Record>> + 'r {
        self.record_bufs(header)
    }
}

/// Helper function that opens one VCF reader at the given path.
pub async fn open_vcf_reader(path_in: &str) -> Result<AsyncVcfReader, anyhow::Error> {
    Ok(noodles::vcf::AsyncReader::new(
        open_read_maybe_bgzf(path_in)
            .await
            .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))?,
    ))
}

/// Helper function that opens a list of paths as VCF readers.
pub async fn open_vcf_readers(paths: &[String]) -> Result<Vec<AsyncVcfReader>, anyhow::Error> {
    Ok(join_all(paths.iter().map(open_read_maybe_gz))
        .await
        .into_iter()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))?
        .into_iter()
        .map(noodles::vcf::AsyncReader::new)
        .collect::<Vec<_>>())
}

/// Alias for the async vcf writer that we use.
pub type AsyncVcfWriter = noodles::vcf::AsyncWriter<std::pin::Pin<Box<dyn AsyncWrite>>>;

/// Helper function that opens one VCF write at the given path.
pub async fn open_vcf_writer(path_out: &str) -> Result<AsyncVcfWriter, anyhow::Error> {
    Ok(noodles::vcf::AsyncWriter::new(
        open_write_maybe_bgzf(path_out)
            .await
            .map_err(|e| anyhow::anyhow!("could not build VCF writer: {}", e))?,
    ))
}
