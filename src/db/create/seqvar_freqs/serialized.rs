//! Serialized population frequencies.

use std::{ops::Deref, str::FromStr};

use byteorder::{ByteOrder, LittleEndian};
use noodles::vcf::{
    self,
    header::info::{key::Other as InfoOther, key::Standard as InfoStandard, Key as InfoKey},
    record::Chromosome,
    Record as VcfRecord,
};

/// Record type for storing AN, AC_hom, AC_het counts for chrMT.
#[derive(Default, Debug, PartialEq, Eq, Clone)]
pub struct MtCounts {
    /// Total number of alleles.
    pub an: u32,
    /// Number of homoplasmic alleles.
    pub ac_hom: u32,
    /// Number of heteroplasmic alleles.
    pub ac_het: u32,
}

impl MtCounts {
    /// Create from the given VCF record.
    pub fn from_vcf(value: &VcfRecord) -> Self {
        let ac_hom = match value
            .info()
            .get(&InfoKey::Other(
                InfoOther::from_str("AC_hom").expect("Invalid key: AC_hom?"),
            ))
            .unwrap()
            .unwrap()
        {
            vcf::record::info::field::Value::Integer(ac_hom) => *ac_hom as u32,
            _ => panic!("invalid type for AC_hom"),
        };
        let ac_het = match value
            .info()
            .get(&InfoKey::Other(
                InfoOther::from_str("AC_het").expect("Invalid key: AC_het?"),
            ))
            .unwrap()
            .unwrap()
        {
            vcf::record::info::field::Value::Integer(ac_het) => *ac_het as u32,
            _ => panic!("invalid type for AC_het"),
        };
        let an = match value
            .info()
            .get(&InfoKey::Standard(InfoStandard::TotalAlleleCount))
            .unwrap()
            .unwrap()
        {
            vcf::record::info::field::Value::Integer(an) => *an as u32,
            _ => panic!("invalid type for AN"),
        };

        MtCounts { ac_hom, ac_het, an }
    }

    /// Read from buffer.
    pub fn from_buf(buf: &[u8]) -> Self {
        Self {
            an: LittleEndian::read_u32(&buf[0..4]),
            ac_hom: LittleEndian::read_u32(&buf[4..8]),
            ac_het: LittleEndian::read_u32(&buf[8..12]),
        }
    }

    /// Write to buffer.
    pub fn to_buf(&self, buf: &mut [u8]) {
        LittleEndian::write_u32(&mut buf[0..4], self.an);
        LittleEndian::write_u32(&mut buf[4..8], self.an);
        LittleEndian::write_u32(&mut buf[8..12], self.an);
    }
}

/// Record type for the "mitochondrial" column family.
pub struct MtRecord {
    /// Counts from gnomAD mtDNA.
    pub gnomad_mtdna: MtCounts,
    /// Counts from HelixMtDb.
    pub helix_mtdb: MtCounts,
}

impl MtRecord {
    /// Read from buffer.
    pub fn from_buf(buf: &[u8]) -> Self {
        Self {
            gnomad_mtdna: MtCounts::from_buf(&buf[0..12]),
            helix_mtdb: MtCounts::from_buf(&buf[12..24]),
        }
    }

    /// Write to buffer.
    pub fn to_buf(&self, buf: &mut [u8]) {
        self.gnomad_mtdna.to_buf(&mut buf[0..12]);
        self.helix_mtdb.to_buf(&mut buf[12..24]);
    }
}

/// A chromosomal change `CHROM-POS-REF-ALT`.
#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct VcfVar {
    pub chrom: String,
    pub pos: u32,
    pub reference: String,
    pub alternative: String,
}

impl VcfVar {
    /// Create from the given VCF record.
    pub fn from_vcf(value: &VcfRecord) -> Self {
        let chrom = match value.chromosome() {
            Chromosome::Name(name) | Chromosome::Symbol(name) => name.to_owned(),
        };
        let pos: usize = value.position().into();
        let pos = pos as u32;
        let reference = value.reference_bases().to_string();
        let alternative = value.alternate_bases().deref()[0].to_string();

        VcfVar {
            chrom,
            pos,
            reference,
            alternative,
        }
    }
}

impl Into<Vec<u8>> for VcfVar {
    fn into(self) -> Vec<u8> {
        let mut result = Vec::new();

        result.extend_from_slice(chrom_name_to_key(&self.chrom).as_bytes());
        result.extend_from_slice(&self.pos.to_be_bytes());
        result.extend_from_slice(self.reference.as_bytes());
        result.push(b'>');
        result.extend_from_slice(self.alternative.as_bytes());

        result
    }
}

/// Convert chromosome to key in RocksDB.
pub fn chrom_name_to_key(name: &str) -> String {
    let chrom = if name.starts_with("chr") {
        &name[3..]
    } else {
        &name[..]
    };
    let chrom = if chrom == "M" {
        String::from("MT")
    } else if "XY".contains(chrom) {
        format!(" {}", chrom)
    } else {
        String::from(chrom)
    };
    assert!(chrom.len() <= 2);
    assert!(chrom.len() >= 1);
    if chrom.len() == 1 {
        format!("0{}", chrom)
    } else {
        chrom.to_string()
    }
}

/// Convert from RocksDB chromosome key part to chromosome name.
pub fn chrom_key_to_name(key: &str) -> String {
    assert!(key.len() == 2);
    if key.starts_with("0") || key.starts_with(" ") {
        key[1..].to_string()
    } else {
        key.to_string()
    }
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use noodles_util::variant::reader::Builder as VariantReaderBuilder;

    use super::*;

    #[test]
    fn test_chrom_name_to_key() {
        assert_eq!(chrom_name_to_key("chr1"), "01");
        assert_eq!(chrom_name_to_key("chr21"), "21");
        assert_eq!(chrom_name_to_key("chrX"), " X");
        assert_eq!(chrom_name_to_key("chrY"), " Y");
        assert_eq!(chrom_name_to_key("chrM"), "MT");
        assert_eq!(chrom_name_to_key("chrMT"), "MT");

        assert_eq!(chrom_name_to_key("1"), "01");
        assert_eq!(chrom_name_to_key("21"), "21");
        assert_eq!(chrom_name_to_key("X"), " X");
        assert_eq!(chrom_name_to_key("Y"), " Y");
        assert_eq!(chrom_name_to_key("M"), "MT");
        assert_eq!(chrom_name_to_key("MT"), "MT");
    }

    #[test]
    fn test_chrom_key_to_name() {
        assert_eq!(chrom_key_to_name("01"), "1");
        assert_eq!(chrom_key_to_name("21"), "21");
        assert_eq!(chrom_key_to_name(" X"), "X");
        assert_eq!(chrom_key_to_name(" Y"), "Y");
        assert_eq!(chrom_key_to_name("MT"), "MT");
    }

    #[test]
    fn test_vcf_var_from_vcf() -> Result<(), anyhow::Error> {
        let path = "tests/data/db/create/seqvar_freqs/helix.chrM.vcf";
        let mut reader = VariantReaderBuilder::default().build_from_path(path)?;
        let header = reader.read_header()?;
        let mut records = reader.records(&header);
        let vcf_record = records.next().transpose()?.unwrap();

        let vcf_var = VcfVar::from_vcf(&vcf_record);
        assert_eq!(
            vcf_var,
            VcfVar {
                chrom: String::from("chrM"),
                pos: 5,
                reference: String::from("A"),
                alternative: String::from("C")
            }
        );

        Ok(())
    }

    #[test]
    fn test_mtcounts_from_vcf() -> Result<(), anyhow::Error> {
        let path = "tests/data/db/create/seqvar_freqs/helix.chrM.vcf";
        let mut reader = VariantReaderBuilder::default().build_from_path(path)?;
        let header = reader.read_header()?;
        let mut records = reader.records(&header);
        let vcf_record = records.next().transpose()?.unwrap();

        let mt_counts = MtCounts::from_vcf(&vcf_record);
        assert_eq!(
            mt_counts,
            MtCounts {
                an: 196554,
                ac_hom: 1,
                ac_het: 0
            }
        );

        Ok(())
    }

    #[test]
    fn test_vcf_var_lexicographic_sorting() {
        let var1 = VcfVar {
            chrom: String::from("11"),
            pos: 10,
            reference: String::from("T"),
            alternative: String::from("TC"),
        };
        let var2 = VcfVar {
            chrom: String::from("11"),
            pos: 10,
            reference: String::from("T"),
            alternative: String::from("TG"),
        };
        let var3 = VcfVar {
            chrom: String::from("11"),
            pos: 10,
            reference: String::from("TC"),
            alternative: String::from("T"),
        };
        let var4 = VcfVar {
            chrom: String::from("11"),
            pos: 10,
            reference: String::from("TG"),
            alternative: String::from("T"),
        };

        assert!(var1 < var2);
        assert!(var2 < var3);
        assert!(var3 < var4);
    }
}
