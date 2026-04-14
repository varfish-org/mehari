//! Tokio-based async common I/O code.

use async_compression::tokio::bufread::GzipDecoder;
use noodles_bgzf;
use std::path::Path;
use std::pin::Pin;
use tokio::fs::File;
use tokio::io::{AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncWrite, BufReader};

/// Returns whether the path looks like a gzip or bgzip file.
pub fn is_gz<P>(path: P) -> bool
where
    P: AsRef<Path>,
{
    [Some(Some("gz")), Some(Some("bgz"))].contains(&path.as_ref().extension().map(|s| s.to_str()))
}

/// Transparently open a file with gzip decoder for reading.
///
/// Note that decoding of multi-member gzip files is automatically supported, as is needed for
/// `bgzip`` files.
///
/// # Arguments
///
/// * `path` - A path to the file to open.
pub async fn open_read_maybe_gz<P>(
    path: P,
) -> Result<Pin<Box<dyn AsyncBufRead + Send>>, anyhow::Error>
where
    P: AsRef<Path>,
{
    if path.as_ref().to_str() == Some("-") {
        let mut stdin = BufReader::new(tokio::io::stdin());

        // Peek at the first 2 bytes to check for GZIP magic number
        let buf = stdin.fill_buf().await?;
        let is_gzipped = buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b;

        if is_gzipped {
            tracing::info!("Opening stdin as gzip for reading (async)");
            let decoder = async_compression::tokio::bufread::GzipDecoder::new(stdin);
            return Ok(Box::pin(BufReader::new(decoder)));
        } else {
            tracing::info!("Opening stdin as plain text for reading (async)");
            return Ok(Box::pin(stdin));
        }
    }

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

/// Transparently open a file or stdin, sniffing for BGZF/GZIP magic bytes.
/// Returns a BufRead stream that might be a multithreaded BGZF decoder.
pub async fn open_read_maybe_bgzf<P>(
    path: P,
) -> Result<Pin<Box<dyn AsyncBufRead + Send>>, anyhow::Error>
where
    P: AsRef<Path>,
{
    let mut bufreader = if path.as_ref().to_str() == Some("-") {
        BufReader::new(Box::pin(tokio::io::stdin()) as Pin<Box<dyn AsyncRead + Send>>)
    } else {
        let file = tokio::fs::File::open(path.as_ref()).await?;
        BufReader::new(Box::pin(file) as Pin<Box<dyn AsyncRead + Send>>)
    };

    // Peek at the first 2 bytes to check for GZIP magic number
    let buf = bufreader.fill_buf().await?;
    let is_gzipped = buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b;

    if is_gzipped {
        tracing::debug!("Detected bgzip/gzip stream. Using noodles-bgzf default worker count.");
        let bgzf_reader =
            noodles_bgzf::r#async::io::reader::Builder::default().build_from_reader(bufreader);
        Ok(Box::pin(bgzf_reader))
    } else {
        Ok(Box::pin(bufreader))
    }
}

/// Transparently open a file or stdout, compressing if the extension suggests it.
pub async fn open_write_maybe_bgzf<P>(
    path: P,
) -> Result<Pin<Box<dyn AsyncWrite + Send>>, anyhow::Error>
where
    P: AsRef<Path>,
{
    if path.as_ref().to_str() == Some("-") {
        Ok(Box::pin(tokio::io::stdout()))
    } else {
        let file = File::create(path.as_ref()).await?;
        let ext = path
            .as_ref()
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("");

        let is_gzipped = ["gz", "bgz", "bcf"].contains(&ext);
        if is_gzipped {
            let bgzf_writer =
                noodles_bgzf::r#async::io::writer::Builder::default().build_from_writer(file);
            Ok(Box::pin(bgzf_writer))
        } else {
            Ok(Box::pin(file))
        }
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
        use std::io::Read;
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
        let mut file_buffer: Vec<u8> = Vec::new();
        std::fs::File::open(&tmp_file_path)?.read_to_end(&mut file_buffer)?;
        hxdmp::hexdump(&file_buffer, &mut buffer)?;
        insta::assert_snapshot!(String::from_utf8_lossy(&buffer));

        Ok(())
    }
}
