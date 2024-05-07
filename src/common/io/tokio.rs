//! Tokio-based async common I/O code.

use async_compression::tokio::bufread::GzipDecoder;
use noodles::bgzf;
use std::path::Path;
use std::pin::Pin;
use tokio::fs::File;
use tokio::io::{AsyncBufRead, AsyncWrite, BufReader, BufWriter};

use crate::common::io::std::is_gz;

/// Transparently open a file with gzip decoder for reading.
///
/// Note that decoding of multi-member gzip files is automatically supported, as is needed for
/// `bgzip`` files.
///
/// # Arguments
///
/// * `path` - A path to the file to open.
pub async fn open_read_maybe_gz<P>(path: P) -> Result<Pin<Box<dyn AsyncBufRead>>, anyhow::Error>
where
    P: AsRef<Path>,
{
    let path_is_gzip = is_gz(path.as_ref());
    tracing::trace!(
        "Opening {} as {} for reading (async)",
        path.as_ref().display(),
        if path_is_gzip {
            "gzip (allow multi-member)"
        } else {
            "plain text"
        }
    );
    let file = File::open(path.as_ref())
        .await
        .map_err(|e| anyhow::anyhow!("could not open file {}: {}", path.as_ref().display(), e))?;

    if path_is_gzip {
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

/// Transparently open a file with bgzip encoder for writing.
///
/// Note that decoding of multi-member gzip files is automatically supported, as is needed for
/// `bgzip`` files.
///
/// # Arguments
///
/// * `path` - A path to the file to open.
pub async fn open_write_maybe_bgzf<P>(path: P) -> Result<Pin<Box<dyn AsyncWrite>>, anyhow::Error>
where
    P: AsRef<Path>,
{
    let path_is_gzip = is_gz(path.as_ref());
    tracing::trace!(
        "Opening {} as {} for writing (async)",
        path.as_ref().display(),
        if path_is_gzip {
            "bgzip (block gzip)"
        } else {
            "plain text"
        }
    );
    let file = File::create(path.as_ref())
        .await
        .map_err(|e| anyhow::anyhow!("could not open file {}: {}", path.as_ref().display(), e))?;

    if path_is_gzip {
        Ok(Box::pin(BufWriter::new(
            bgzf::r#async::writer::Writer::new(file),
        )))
    } else {
        Ok(Box::pin(BufWriter::new(file)))
    }
}

#[cfg(test)]
mod test {
    use tokio::io::{AsyncReadExt, AsyncWriteExt};

    #[rstest::rstest]
    #[case("14kb.txt")]
    #[case("14kb.txt.gz")]
    #[case("14kb.txt.bgz")]
    #[tokio::test]
    async fn open_read_maybe_gz(#[case] path: &str) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", path);
        // Note that the 14kb.txt file contains about 14 KB of data so bgz will have multiple 4KB
        // blocks.

        let mut reader = super::open_read_maybe_gz(&format!("tests/common/io/{}", path)).await?;
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).await?;

        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }

    #[rstest::rstest]
    #[case("14kb.txt")]
    #[case("14kb.txt.gz")]
    #[case("14kb.txt.bgz")]
    #[tokio::test]
    async fn open_write_maybe_bgzf(#[case] filename: &str) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", filename);
        // Note that the 14kb.txt file contains about 14 KB of data so bgz will have multiple 4KB
        // blocks.

        let tmp_dir = temp_testdir::TempDir::default();
        let tmp_file_path = tmp_dir.join(filename);

        {
            let mut writer = super::open_write_maybe_bgzf(&tmp_file_path).await?;
            for i in 1..3000 {
                writer.write_all(format!("{}\n", i).as_bytes()).await?;
            }

            writer.flush().await?;
            writer.shutdown().await?;
        }
        std::thread::sleep(std::time::Duration::from_millis(100));

        let mut buffer: Vec<u8> = Vec::new();
        hxdmp::hexdump(
            &crate::common::io::std::read_to_bytes(&tmp_file_path)?,
            &mut buffer,
        )?;
        insta::assert_snapshot!(String::from_utf8_lossy(&buffer));

        Ok(())
    }
}
