//! Helper code for using noodles.

use std::path::Path;

use futures::StreamExt;
use futures::stream::BoxStream;
use noodles::{bcf, vcf, vcf::Header, vcf::variant::RecordBuf};
use tokio::io::{AsyncBufRead, AsyncWrite, AsyncWriteExt};

use crate::annotate::seqvars::{AnnotatedVariant, AsyncAnnotatedVariantWriter};
use crate::common::io::tokio::{open_read_maybe_bgzf, open_write_maybe_bgzf};

/// Alias for the async vcf reader type.
pub type AsyncVcfReader = vcf::r#async::io::Reader<std::pin::Pin<Box<dyn AsyncBufRead + Send>>>;

/// Alias for the async bcf reader type.
pub type AsyncBcfReader = bcf::r#async::io::Reader<std::pin::Pin<Box<dyn AsyncBufRead + Send>>>;

/// Alias for the async vcf writer type.
pub type AsyncVcfWriter = vcf::r#async::io::Writer<std::pin::Pin<Box<dyn AsyncWrite + Send>>>;

/// Alias for the async bcf writer type.
pub type AsyncBcfWriter = bcf::r#async::io::Writer<std::pin::Pin<Box<dyn AsyncWrite + Send>>>;

/// Helper function that opens one VCF reader at the given path.
pub async fn open_vcf_reader(path: impl AsRef<Path>) -> anyhow::Result<AsyncVcfReader> {
    let stream = open_read_maybe_bgzf(path).await?;
    Ok(vcf::r#async::io::Reader::new(stream))
}

/// Helper function that opens one BCF reader at the given path.
pub async fn open_bcf_reader(path: impl AsRef<Path>) -> anyhow::Result<AsyncBcfReader> {
    let stream = open_read_maybe_bgzf(path).await?;
    Ok(AsyncBcfReader::from(stream))
}

/// Helper function that opens one VCF writer at the given path.
pub async fn open_vcf_writer(path: impl AsRef<Path>) -> anyhow::Result<AsyncVcfWriter> {
    let stream = open_write_maybe_bgzf(path).await?;
    Ok(vcf::r#async::io::Writer::new(stream))
}

/// Helper function that opens one BCF writer at the given path.
pub async fn open_bcf_writer(path: impl AsRef<Path>) -> anyhow::Result<AsyncBcfWriter> {
    let stream = open_write_maybe_bgzf(path).await?;
    Ok(AsyncBcfWriter::from(stream))
}

/// Enum for the different types of variant readers.
pub enum VariantReader {
    Vcf(AsyncVcfReader),
    Bcf(AsyncBcfReader),
}

/// Trait for VCF/BCF readers to allow generic usage and mocking.
pub trait NoodlesVariantReader {
    #[allow(async_fn_in_trait)]
    async fn read_header(&mut self) -> anyhow::Result<Header>;

    #[allow(async_fn_in_trait)]
    async fn records<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> BoxStream<'a, anyhow::Result<RecordBuf>>;
}

impl NoodlesVariantReader for VariantReader {
    async fn read_header(&mut self) -> anyhow::Result<Header> {
        match self {
            VariantReader::Vcf(r) => r.read_header().await.map_err(Into::into),
            VariantReader::Bcf(r) => r.read_header().await.map_err(Into::into),
        }
    }

    async fn records<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> BoxStream<'a, anyhow::Result<RecordBuf>> {
        match self {
            VariantReader::Vcf(r) => r
                .record_bufs(header)
                .map(|res| res.map_err(anyhow::Error::from))
                .boxed(),
            VariantReader::Bcf(r) => r
                .records()
                .map(move |res| {
                    res.map_err(anyhow::Error::from).and_then(|rec| {
                        RecordBuf::try_from_variant_record(header, &rec)
                            .map_err(anyhow::Error::from)
                    })
                })
                .boxed(),
        }
    }
}

/// Opens a variant reader from the given path, detecting format by extension.
pub async fn open_variant_reader(path: impl AsRef<Path>) -> anyhow::Result<VariantReader> {
    match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("bcf") => open_bcf_reader(path).await.map(VariantReader::Bcf),
        _ => open_vcf_reader(path).await.map(VariantReader::Vcf),
    }
}

/// Enum for the different types of variant writers.
pub enum VariantWriter {
    Vcf(AsyncVcfWriter),
    Bcf(AsyncBcfWriter),
}

/// Opens a variant writer from the given path, detecting format by extension.
pub async fn open_variant_writer(path: impl AsRef<Path>) -> anyhow::Result<VariantWriter> {
    match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("bcf") => open_bcf_writer(path).await.map(VariantWriter::Bcf),
        _ => open_vcf_writer(path).await.map(VariantWriter::Vcf),
    }
}

impl AsyncAnnotatedVariantWriter for VariantWriter {
    async fn write_noodles_header(&mut self, header: &Header) -> Result<(), anyhow::Error> {
        match self {
            VariantWriter::Vcf(w) => w.write_header(header).await.map_err(Into::into),
            VariantWriter::Bcf(w) => w.write_header(header).await.map_err(Into::into),
        }
    }

    async fn write_annotated_record(
        &mut self,
        header: &Header,
        record: AnnotatedVariant,
    ) -> Result<(), anyhow::Error> {
        match self {
            VariantWriter::Vcf(w) => w.write_annotated_record(header, record).await,
            VariantWriter::Bcf(w) => w.write_annotated_record(header, record).await,
        }
    }

    async fn shutdown(&mut self) -> Result<(), anyhow::Error> {
        match self {
            VariantWriter::Vcf(w) => w.get_mut().shutdown().await.map_err(Into::into),
            VariantWriter::Bcf(w) => w.get_mut().shutdown().await.map_err(Into::into),
        }
    }
}
