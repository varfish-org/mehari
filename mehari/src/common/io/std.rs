//! Common I/O code using sync I/O.

use flate2::{Compression, bufread::MultiGzDecoder, write::GzEncoder};
use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
    path::Path,
};

/// Returns whether the path looks like a gzip, bgzip, or BCF file based on extension.
pub fn is_gz<P>(path: P) -> bool
where
    P: AsRef<Path>,
{
    path.as_ref()
        .extension()
        .and_then(|s| s.to_str())
        .map(|ext| ["gz", "bgz", "bcf"].contains(&ext))
        .unwrap_or(false)
}

/// Transparently open a file or stdin for reading.
/// Uses magic-byte sniffing (1f 8b) to detect compression, making it extension-agnostic.
pub fn open_read_maybe_bgzf<P>(path: P) -> Result<Box<dyn BufRead + Send>, anyhow::Error>
where
    P: AsRef<Path>,
{
    let mut bufreader: Box<dyn BufRead + Send> = if path.as_ref().to_str() == Some("-") {
        Box::new(BufReader::new(std::io::stdin()))
    } else {
        Box::new(BufReader::new(File::open(path.as_ref())?))
    };

    // Sniff for GZIP magic bytes (1f 8b)
    let buf = bufreader.fill_buf()?;
    let is_gzipped = buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b;

    if is_gzipped {
        // MultiGzDecoder handles multi-member gzip files (like BGZF)
        let decoder = MultiGzDecoder::new(bufreader);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(bufreader)
    }
}

/// Transparently open a file or stdout for writing.
/// Uses the file extension to decide whether to apply GZIP compression.
pub fn open_write_maybe_bgzf<P>(path: P) -> Result<Box<dyn Write + Send>, anyhow::Error>
where
    P: AsRef<Path>,
{
    if path.as_ref().to_str() == Some("-") {
        Ok(Box::new(std::io::stdout()))
    } else {
        let file = File::create(path.as_ref())?;
        if is_gz(path.as_ref()) {
            Ok(Box::new(GzEncoder::new(file, Compression::default())))
        } else {
            Ok(Box::new(file))
        }
    }
}

/// Return `std::io::Lines<>` for buffered reading from a file.
pub fn read_lines<P: AsRef<Path>>(filename: P) -> std::io::Result<std::io::Lines<BufReader<File>>> {
    let file = File::open(filename)?;
    Ok(BufReader::new(file).lines())
}

/// Helper function to read a file to a byte vector.
pub fn read_to_bytes<P>(path: P) -> Result<Vec<u8>, anyhow::Error>
where
    P: AsRef<std::path::Path>,
{
    std::fs::read(path).map_err(|e| anyhow::anyhow!("Failed to read file: {}", e))
}

/// Given a `BufWriter<File>`, flush buffers and sync the file.
#[macro_export]
macro_rules! finalize_buf_writer {
    ($a:expr_2021) => {
        $a.flush()
            .map_err(|e| anyhow::anyhow!("problem flushing buffers: {}", e))?;
        let file = $a
            .into_inner()
            .map_err(|e| anyhow::anyhow!("problem getting inner file: {}", e))?;
        file.sync_all()
            .map_err(|e| anyhow::anyhow!("problem syncing file: {}", e))?;
    };
}

#[cfg(test)]
mod test {
    use std::io::Read;

    #[rstest::rstest]
    #[case("14kb.txt")]
    #[case("14kb.txt.gz")]
    #[case("14kb.txt.bgz")]
    fn open_read_maybe_gz(#[case] path: &str) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", path);
        // Note that the 14kb.txt file contains about 14 KB of data so bgz will have multiple 4KB
        // blocks.

        let mut reader = super::open_read_maybe_bgzf(format!("tests/common/io/{}", path))?;
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;

        insta::assert_snapshot!(String::from_utf8(buf)?);

        Ok(())
    }

    #[rstest::rstest]
    #[case("14kb.txt")]
    #[case("14kb.txt.gz")]
    #[case("14kb.txt.bgz")]
    fn open_write_maybe_bgzf(#[case] filename: &str) -> Result<(), anyhow::Error> {
        crate::common::set_snapshot_suffix!("{}", filename);

        let tmp_dir = temp_testdir::TempDir::default();
        let tmp_file_path = tmp_dir.join(filename);

        {
            let mut writer = super::open_write_maybe_bgzf(&tmp_file_path)?;
            for i in 1..3000 {
                writer.write_all(format!("{}\n", i).as_bytes())?;
            }

            writer.flush()?;
        }
        std::thread::sleep(std::time::Duration::from_millis(100));

        let mut buffer: Vec<u8> = Vec::new();
        hxdmp::hexdump(&super::read_to_bytes(&tmp_file_path)?, &mut buffer)?;
        insta::assert_snapshot!(String::from_utf8_lossy(&buffer));

        Ok(())
    }

    #[test]
    fn read_lines() -> Result<(), anyhow::Error> {
        let lines =
            super::read_lines("tests/common/io/lines.txt")?.collect::<Result<Vec<_>, _>>()?;

        insta::assert_yaml_snapshot!(lines);

        Ok(())
    }
}
