//! Code for reading maelstrom coverage and mapping quality VCF files.

use std::{
    ops::Range,
    path::{Path, PathBuf},
};

use noodles::core::{Position, Region};
use noodles::vcf;
use noodles::vcf::variant::record::info::field::Value::Integer;
use noodles::vcf::variant::record::samples::series;

/// Read access for maelstrom coverage/mapping quality VCF files.
///
/// This only works robustly for WGS data.  Also, it is assumed that each VCF
/// file contains exactly one sample.
pub struct Reader {
    /// Path to the VCF file.
    pub path: PathBuf,
    /// Name of the single sample in the VCF file.
    pub sample_name: String,
    /// The internal reader.
    pub reader: vcf::io::IndexedReader<noodles::bgzf::Reader<std::fs::File>>,
    /// The header from the VCF file.
    pub header: vcf::header::Header,
}

impl Reader {
    /// Construct reader with given path and sample name.
    ///
    /// The VCF file at `p` must have an index.
    ///
    /// # Arguments
    ///
    /// * `p` - Path to the tabix indexed and bgzip-compressed VCF file.
    pub fn from_path<P>(p: P) -> Result<Self, anyhow::Error>
    where
        P: AsRef<Path> + Clone,
    {
        let path = p.as_ref().to_path_buf();
        let mut reader = vcf::io::indexed_reader::Builder::default().build_from_path(&path)?;
        let header = reader.read_header()?;
        let sample = if header.sample_names().len() == 1 {
            header
                .sample_names()
                .iter()
                .next()
                .expect("just checked for ==1 sample")
                .clone()
        } else {
            anyhow::bail!(
                "VCF file {} must contain exactly one sample",
                p.as_ref().display()
            );
        };

        Ok(Self {
            path,
            sample_name: sample,
            reader,
            header,
        })
    }

    /// Read coverage and mapping quality for given range and compute means.
    ///
    /// # Arguments
    ///
    /// * `chrom` - Chromosome name.
    /// * `range` - 1-based range to read on `chrom`.
    ///
    /// # Returns
    ///
    /// The mean coverage and mean mapping quality bundled into one `Record`.
    pub fn read(&mut self, chrom: &str, range: Range<i32>) -> Result<Record, anyhow::Error> {
        // Jump to the given region.
        let start: usize = range.start.try_into()?;
        let pos_start = Position::try_from(start)?;
        let end: usize = range.end.try_into()?;
        let pos_end = Position::try_from(end)?;
        let region = Region::new(chrom, pos_start..=pos_end);
        let query = self.reader.query(&self.header, &region)?;
        let header = &self.header;

        // Read in the records and compute the mean coverage and mapping quality.  The records
        // describe the mean statistics over a window.  We go all the way to compute fractions
        // of windows.
        let mut cov_sum = 0f64;
        let mut mq_sum = 0f64;
        let mut count = 0f64;

        for result in query {
            let record = result?;

            let window_end = match record.info().get(header, "END").transpose()? {
                Some(Some(Integer(window_end))) => window_end as usize,
                _ => {
                    anyhow::bail!("missing INFO/END in record");
                }
            };
            let window_start: usize = record
                .variant_start()
                .expect("Telomeric breakends not supported")?
                .get();
            let window_size = window_end - window_start + 1;

            // Use first and last window values only in fractions.
            let covered = if window_start < start {
                window_size - (start - window_start)
            } else if window_end > end {
                window_size - (window_end - end)
            } else {
                window_size
            };
            let covered: i32 = covered.try_into()?;
            let covered: f64 = covered.into();
            let window_size: i32 = window_size.try_into()?;
            let window_size: f64 = window_size.into();
            let factor = covered / window_size;
            count += factor;

            let sample = record
                .samples()
                .iter()
                .next()
                .expect("just checked for ==1 sample");

            // The simplest way to obtain the genotype keys is to iterate and call `as_ref()` on the
            // key.
            for (key, value) in sample.iter(header).flatten() {
                match (key, value) {
                    ("CV", Some(series::Value::Float(cov))) => {
                        cov_sum += factor * cov as f64;
                    }
                    ("MQ", Some(series::Value::Float(mq))) => {
                        mq_sum += factor * mq as f64;
                    }
                    // Ignore all other keys.
                    _ => (),
                }
            }
        }

        Ok(if count > 0f64 {
            Record {
                mean_coverage: cov_sum / count,
                mean_mapq: mq_sum / count,
            }
        } else {
            Default::default()
        })
    }
}

/// Result for `Reader::read`.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Record {
    /// Mean quality.
    pub mean_coverage: f64,
    /// Mean mapping quality.
    pub mean_mapq: f64,
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use super::Reader;

    #[test]
    fn reader_from_path() -> Result<(), anyhow::Error> {
        let path = "tests/data/annotate/strucvars/maelstrom/example.SAMPLE.cov.vcf.gz";
        let reader = Reader::from_path(path)?;
        assert_eq!(reader.sample_name, "SAMPLE");
        assert_eq!(format!("{}", reader.path.display()), path);

        Ok(())
    }

    #[test]
    fn reader_read() -> Result<(), anyhow::Error> {
        let path = "tests/data/annotate/strucvars/maelstrom/example.SAMPLE.cov.vcf.gz";
        let mut reader = Reader::from_path(path)?;

        {
            let record = reader.read("1", 10_001..20_000)?;
            assert_eq!(record.mean_coverage, 1.0);
            assert_eq!(record.mean_mapq, 40.0);
        }

        {
            let record = reader.read("1", 10_501..20_000)?;
            assert_eq!(record.mean_coverage, 1.0);
            assert_eq!(record.mean_mapq, 40.0);
        }

        {
            let record = reader.read("1", 10_001..19_500)?;
            assert_eq!(record.mean_coverage, 1.0);
            assert_eq!(record.mean_mapq, 40.0);
        }

        {
            let record = reader.read("1", 10_501..19_500)?;
            assert_eq!(record.mean_coverage, 1.0);
            assert_eq!(record.mean_mapq, 40.0);
        }

        Ok(())
    }
}
