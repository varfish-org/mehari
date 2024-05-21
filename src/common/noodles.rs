//! Helper code for using noodles.

use futures::future::join_all;
use noodles_vcf as vcf;
use tokio::io::{AsyncBufRead, AsyncWrite};

use super::io::{tokio::open_read_maybe_gz, tokio::open_write_maybe_bgzf};

/// Alias for the async vcf reader type that we will use.
pub type AsyncVcfReader = vcf::AsyncReader<std::pin::Pin<Box<dyn AsyncBufRead>>>;

/// Alias for the async vcf reader type that we will use.
pub type AsyncBcfReader = bcf::AsyncReader<noodles::bgzf::AsyncReader<tokio::fs::File>>;

/// A variant format reader.
pub(crate) trait AsyncNoodlesReader<R> {
    /// Reads a VCF header.
    async fn read_variant_header(&mut self) -> io::Result<Header>;

    /// Returns an iterator over records.
    async fn variant_records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h Header,
    ) -> impl Stream<Item = io::Result<Record>> + 'r;
}

impl<R> AsyncNoodlesReader<R> for AsyncVcfReader {
    async fn read_variant_header(&mut self) -> io::Result<Header> {
        self.read_header().await
    }

    async fn variant_records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h Header,
    ) -> impl Stream<Item = io::Result<Record>> + 'r {
        self.records(header)
    }
}

impl<R> AsyncNoodlesReader<R> for AsyncBcfReader {
    async fn read_variant_header(&mut self) -> io::Result<Header> {
        self.read_header().await
    }

    async fn variant_records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h Header,
    ) -> impl Stream<Item = io::Result<Record>> + 'r {
        self.records()
    }
}

/// Helper function that opens one VCF reader at the given path.
pub async fn open_vcf_reader(path_in: &str) -> Result<AsyncVcfReader, anyhow::Error> {
    Ok(vcf::AsyncReader::new(
        open_read_maybe_gz(path_in)
            .await
            .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))?,
    ))
}

pub async fn open_bcf_reader(path_in: impl AsRef<Path>) -> Result<AsyncBcfReader, anyhow::Error> {
    Ok(tokio::fs::File::open(path_in)
        .await
        .map(bcf::AsyncReader::new)?)
}

pub async fn open_reader<R>(
    path_in: impl AsRef<Path>,
) -> Result<std::pin::Pin<Box<dyn AsyncNoodlesReader<R>>>, anyhow::Error> {
    todo!()
}

/// Helper function that opens a list of paths as VCF readers.
pub async fn open_vcf_readers(paths: &[String]) -> Result<Vec<AsyncVcfReader>, anyhow::Error> {
    Ok(join_all(paths.iter().map(open_read_maybe_gz))
        .await
        .into_iter()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))?
        .into_iter()
        .map(vcf::AsyncReader::new)
        .collect::<Vec<_>>())
}

/// Alias for the async vcf writer that we use.
pub type AsyncVcfWriter = vcf::AsyncWriter<std::pin::Pin<Box<dyn AsyncWrite>>>;

/// Helper function that opens one VCF write at the given path.
pub async fn open_vcf_writer(path_out: &str) -> Result<AsyncVcfWriter, anyhow::Error> {
    Ok(vcf::AsyncWriter::new(
        open_write_maybe_bgzf(path_out)
            .await
            .map_err(|e| anyhow::anyhow!("could not build VCF writer: {}", e))?,
    ))
}
