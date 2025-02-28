//! Common I/O code using sync I/O.

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use flate2::{bufread::MultiGzDecoder, write::GzEncoder, Compression};

/// Returns whether the path looks like a gzip or bgzip file.
pub fn is_gz<P>(path: P) -> bool
where
    P: AsRef<Path>,
{
    [Some(Some("gz")), Some(Some("bgz"))].contains(&path.as_ref().extension().map(|s| s.to_str()))
}

/// Transparently open a file with bgzip encoder for writing.
///
/// Note that decoding of multi-member gzip files is automatically supported, as is needed for
/// `bgzip`` files.
///
/// # Arguments
///
/// * `path` - A path to the file to open.
pub fn open_read_maybe_gz<P>(path: P) -> Result<Box<dyn BufRead>, anyhow::Error>
where
    P: AsRef<Path>,
{
    if is_gz(path.as_ref()) {
        tracing::trace!("Opening {:?} as gzip for reading", path.as_ref());
        let file = File::open(path)?;
        let bufreader = BufReader::new(file);
        let decoder = MultiGzDecoder::new(bufreader);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        tracing::trace!("Opening {:?} as plain text for reading", path.as_ref());
        let file = File::open(path).map(BufReader::new)?;
        Ok(Box::new(BufReader::new(file)))
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
pub fn open_write_maybe_bgzf<P>(path: P) -> Result<Box<dyn Write>, anyhow::Error>
where
    P: AsRef<Path>,
{
    if path.as_ref().extension().map(|s| s.to_str()) == Some(Some("gz")) {
        tracing::trace!("Opening {:?} as gzip for writing", path.as_ref());
        let file = File::create(path)?;
        let bufwriter = BufWriter::new(file);
        let encoder = GzEncoder::new(bufwriter, Compression::default());
        Ok(Box::new(encoder))
    } else {
        tracing::trace!("Opening {:?} as plain text for writing", path.as_ref());
        let file = File::create(path)?;
        Ok(Box::new(file))
    }
}

/// Return `std::io::Lines<>` for buffered reading from a file.
///
/// The output is wrapped in a Result to allow matching on errors
/// Returns an Iterator to the Reader of the lines of the file.
pub fn read_lines<P: AsRef<Path>>(filename: P) -> std::io::Result<std::io::Lines<BufReader<File>>> {
    let file = File::open(filename)?;
    Ok(BufReader::new(file).lines())
}

/// Helper function to read a file to a byte vector.
pub fn read_to_bytes<P>(path: P) -> Result<Vec<u8>, anyhow::Error>
where
    P: AsRef<std::path::Path>,
{
    use std::io::Read;

    let mut f = std::fs::File::open(&path).expect("no file found");
    let metadata = std::fs::metadata(&path).expect("unable to read metadata");
    let mut buffer = vec![0; metadata.len() as usize];
    f.read_exact(&mut buffer).expect("buffer overflow");
    Ok(buffer)
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

        let mut reader = super::open_read_maybe_gz(format!("tests/common/io/{}", path))?;
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
