//! Helper code for using noodles.

use std::path::Path;

use futures::future::join_all;
use futures::stream::LocalBoxStream;
use futures::StreamExt;
use noodles::vcf;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::Header;
use tokio::io::{AsyncBufRead, AsyncRead, AsyncWrite};

use crate::common::io::tokio::open_read_maybe_bgzf;

use super::io::{tokio::open_read_maybe_gz, tokio::open_write_maybe_bgzf};

/// Alias for the async vcf reader type that we will use.
pub type AsyncVcfReader = noodles::vcf::AsyncReader<std::pin::Pin<Box<dyn AsyncBufRead>>>;
pub type AsyncBcfReader =
    noodles::bcf::AsyncReader<noodles::bgzf::AsyncReader<std::pin::Pin<Box<dyn AsyncRead>>>>;

/// Helper function that opens one VCF reader at the given path.
pub async fn open_vcf_reader(path: impl AsRef<Path>) -> Result<AsyncVcfReader, anyhow::Error> {
    Ok(vcf::AsyncReader::new(
        open_read_maybe_bgzf(path)
            .await
            .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))?,
    ))
}

pub async fn open_bcf_reader(path: impl AsRef<Path>) -> Result<AsyncBcfReader, anyhow::Error> {
    Ok(noodles::bcf::AsyncReader::new(
        tokio::fs::File::open(path)
            .await
            .map(Box::pin)
            .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))?,
    ))
}

pub async fn open_variant_reader(path: impl AsRef<Path>) -> anyhow::Result<VariantReader> {
    match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("bcf") => open_bcf_reader(path).await.map(VariantReader::Bcf),
        Some("gz") | Some("vcf") | _ => open_vcf_reader(path).await.map(VariantReader::Vcf),
    }
}

pub(crate) enum VariantReader {
    Vcf(AsyncVcfReader),
    Bcf(AsyncBcfReader),
}

pub(crate) trait NoodlesVariantReader {
    async fn read_header(&mut self) -> tokio::io::Result<Header>;
    async fn records<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> LocalBoxStream<std::io::Result<RecordBuf>>;
}

impl NoodlesVariantReader for VariantReader {
    async fn read_header(&mut self) -> std::io::Result<Header> {
        match self {
            VariantReader::Vcf(r) => r.read_header().await,
            VariantReader::Bcf(r) => r.read_header().await,
        }
    }

    async fn records<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> LocalBoxStream<std::io::Result<RecordBuf>> {
        match self {
            VariantReader::Vcf(r) => r.record_bufs(header).boxed_local(),
            VariantReader::Bcf(r) => r
                .records()
                .map(|r| r.and_then(|r| RecordBuf::try_from_variant_record(header, &r)))
                .boxed_local(),
        }
    }
}

/// Helper function that opens a list of paths as VCF readers.
#[allow(dead_code)]
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
