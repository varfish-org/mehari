//! Helper code for using noodles.

use std::path::Path;

use anyhow::Error;
use futures::future::join_all;
use futures::stream::LocalBoxStream;
use futures::StreamExt;
use noodles::bcf;
use noodles::vcf;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::Header;
use tokio::io::{AsyncBufRead, AsyncRead, AsyncWrite};

use crate::annotate::seqvars::AsyncAnnotatedVariantWriter;
use crate::common::io::tokio::open_read_maybe_bgzf;

use super::io::{tokio::open_read_maybe_gz, tokio::open_write_maybe_bgzf};

/// Alias for the async vcf reader type that we will use.
pub type AsyncVcfReader = vcf::AsyncReader<std::pin::Pin<Box<dyn AsyncBufRead>>>;
pub type AsyncBcfReader =
    bcf::AsyncReader<noodles::bgzf::AsyncReader<std::pin::Pin<Box<dyn AsyncRead>>>>;

/// Helper function that opens one VCF reader at the given path.
pub async fn open_vcf_reader(path: impl AsRef<Path>) -> Result<AsyncVcfReader, Error> {
    Ok(vcf::AsyncReader::new(
        open_read_maybe_bgzf(path)
            .await
            .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))?,
    ))
}

pub async fn open_bcf_reader(path: impl AsRef<Path>) -> Result<AsyncBcfReader, Error> {
    Ok(bcf::AsyncReader::new(
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

pub enum VariantReader {
    Vcf(AsyncVcfReader),
    Bcf(AsyncBcfReader),
}

pub trait NoodlesVariantReader {
    #[allow(async_fn_in_trait)]
    async fn read_header(&mut self) -> tokio::io::Result<Header>;
    #[allow(async_fn_in_trait)]
    async fn records<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> LocalBoxStream<std::io::Result<RecordBuf>>;
}

impl NoodlesVariantReader for VariantReader {
    #[allow(async_fn_in_trait)]
    async fn read_header(&mut self) -> std::io::Result<Header> {
        match self {
            VariantReader::Vcf(r) => r.read_header().await,
            VariantReader::Bcf(r) => r.read_header().await,
        }
    }

    #[allow(async_fn_in_trait)]
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
pub async fn open_vcf_readers(paths: &[String]) -> Result<Vec<AsyncVcfReader>, Error> {
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
pub async fn open_vcf_writer(path: impl AsRef<Path>) -> Result<AsyncVcfWriter, Error> {
    Ok(vcf::AsyncWriter::new(
        open_write_maybe_bgzf(path)
            .await
            .map_err(|e| anyhow::anyhow!("could not build VCF writer: {}", e))?,
    ))
}

pub type AsyncBcfWriter =
    bcf::AsyncWriter<noodles::bgzf::AsyncWriter<std::pin::Pin<Box<dyn AsyncWrite>>>>;
/// Helper function that opens one VCF write at the given path.
pub async fn open_bcf_writer(path: impl AsRef<Path>) -> Result<AsyncBcfWriter, Error> {
    Ok(bcf::AsyncWriter::new(
        tokio::fs::File::create(path)
            .await
            .map(Box::pin)
            .map_err(|e| anyhow::anyhow!("could not build VCF writer: {}", e))?,
    ))
}

pub async fn open_variant_writer(path: impl AsRef<Path>) -> anyhow::Result<VariantWriter> {
    match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("bcf") => open_bcf_writer(path).await.map(VariantWriter::Bcf),
        Some("gz") | Some("vcf") | _ => open_vcf_writer(path).await.map(VariantWriter::Vcf),
    }
}

trait NoodlesVariantWriter {
    async fn write_header(&mut self, header: &Header) -> std::io::Result<()>;
    async fn write_record(
        &mut self,
        header: &Header,
        record: &impl vcf::variant::Record,
    ) -> std::io::Result<()>;
}

pub enum VariantWriter {
    Vcf(AsyncVcfWriter),
    Bcf(AsyncBcfWriter),
}

impl NoodlesVariantWriter for VariantWriter {
    async fn write_header(&mut self, header: &Header) -> std::io::Result<()> {
        match self {
            VariantWriter::Vcf(w) => w.write_header(header).await,
            VariantWriter::Bcf(w) => w.write_header(header).await,
        }
    }

    async fn write_record(
        &mut self,
        header: &Header,
        record: &impl vcf::variant::Record,
    ) -> std::io::Result<()> {
        match self {
            VariantWriter::Vcf(w) => w.write_variant_record(header, record).await,
            VariantWriter::Bcf(w) => w.write_variant_record(header, record).await,
        }
    }
}

impl AsyncAnnotatedVariantWriter for VariantWriter {
    async fn write_noodles_header(&mut self, header: &Header) -> Result<(), Error> {
        self.write_header(header).await.map_err(Into::into)
    }

    async fn write_noodles_record(
        &mut self,
        header: &Header,
        record: &RecordBuf,
    ) -> Result<(), Error> {
        self.write_record(header, record).await.map_err(Into::into)
    }

    async fn flush(&mut self) -> Result<(), Error> {
        match self {
            VariantWriter::Vcf(r) => r.flush().await,
            VariantWriter::Bcf(r) => r.flush().await,
        }
    }
}
