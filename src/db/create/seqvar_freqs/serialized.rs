//! Serialization of population frequencies.

/// Mitochondrial counts.
pub mod mt {
    use std::str::FromStr;

    use byteorder::{ByteOrder, LittleEndian};
    use noodles::vcf::{
        self,
        header::info::{key::Other as InfoOther, key::Standard as InfoStandard, Key as InfoKey},
        Record as VcfRecord,
    };

    /// Record type for storing AN, AC_hom, AC_het counts for chrMT.
    #[derive(Default, Debug, PartialEq, Eq, Clone)]
    pub struct Counts {
        /// Total number of alleles.
        pub an: u32,
        /// Number of homoplasmic alleles.
        pub ac_hom: u32,
        /// Number of heteroplasmic alleles.
        pub ac_het: u32,
    }

    impl Counts {
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

            Counts { ac_hom, ac_het, an }
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
            LittleEndian::write_u32(&mut buf[4..8], self.ac_hom);
            LittleEndian::write_u32(&mut buf[8..12], self.ac_het);
        }
    }

    /// Record type for the "mitochondrial" column family.
    pub struct Record {
        /// Counts from gnomAD mtDNA.
        pub gnomad_mtdna: Counts,
        /// Counts from HelixMtDb.
        pub helix_mtdb: Counts,
    }

    impl Record {
        /// Read from buffer.
        pub fn from_buf(buf: &[u8]) -> Self {
            Self {
                gnomad_mtdna: Counts::from_buf(&buf[0..12]),
                helix_mtdb: Counts::from_buf(&buf[12..24]),
            }
        }

        /// Write to buffer.
        pub fn to_buf(&self, buf: &mut [u8]) {
            self.gnomad_mtdna.to_buf(&mut buf[0..12]);
            self.helix_mtdb.to_buf(&mut buf[12..24]);
        }
    }
}

/// Autosomal counts.
pub mod auto {
    use std::str::FromStr;

    use byteorder::{ByteOrder, LittleEndian};
    use noodles::vcf::{
        self,
        header::info::{key::Other as InfoOther, key::Standard as InfoStandard, Key as InfoKey},
        Record as VcfRecord,
    };

    /// Record type for storing AN, AC_hom, AC_het counts for autosomal chromosomes.
    #[derive(Default, Debug, PartialEq, Eq, Clone)]
    pub struct Counts {
        /// Total number of alleles.
        pub an: u32,
        /// Number of hom. alt. alleles.
        pub ac_hom: u32,
        /// Number of het. alt. alleles.
        pub ac_het: u32,
    }

    impl Counts {
        /// Create from the given VCF record.
        pub fn from_vcf(value: &VcfRecord) -> Self {
            tracing::trace!("@ {:?}", &value);

            let ac = value
                .info()
                .get(&InfoKey::Other(
                    InfoOther::from_str("AC").expect("Invalid key: AC?"),
                ))
                .unwrap_or_default();
            let ac = if let Some(ac) = ac {
                match ac {
                    vcf::record::info::field::Value::Integer(ac) => *ac as u32,
                    _ => panic!("invalid type for AC"),
                }
            } else {
                0
            };

            let ac_hom = value
                .info()
                .get(&InfoKey::Other(
                    InfoOther::from_str("nhomalt").expect("Invalid key: nhomalt?"),
                ))
                .unwrap_or_default();
            let ac_hom = if let Some(ac_hom) = ac_hom {
                match ac_hom {
                    noodles::vcf::record::info::field::Value::IntegerArray(nhomalt) => {
                        nhomalt[0].unwrap() as u32
                    }
                    _ => panic!("invalid type for nhomalt"),
                }
            } else {
                0
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

            Counts {
                ac_hom,
                ac_het: an - 2 * ac_hom,
                an,
            }
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
            LittleEndian::write_u32(&mut buf[4..8], self.ac_hom);
            LittleEndian::write_u32(&mut buf[8..12], self.ac_het);
        }
    }

    /// Record type for the "autosomal" column family.
    pub struct Record {
        /// Counts from gnomAD exomes.
        pub gnomad_exomes: Counts,
        /// Counts from gnomAD genomes.
        pub gnomad_genomes: Counts,
    }

    impl Record {
        /// Read from buffer.
        pub fn from_buf(buf: &[u8]) -> Self {
            Self {
                gnomad_exomes: Counts::from_buf(&buf[0..16]),
                gnomad_genomes: Counts::from_buf(&buf[16..32]),
            }
        }

        /// Write to buffer.
        pub fn to_buf(&self, buf: &mut [u8]) {
            self.gnomad_exomes.to_buf(&mut buf[0..16]);
            self.gnomad_genomes.to_buf(&mut buf[16..32]);
        }
    }
}

/// Gonomosomal counts.
pub mod xy {
    use std::str::FromStr;

    use byteorder::{ByteOrder, LittleEndian};
    use noodles::vcf::{
        self,
        header::info::{key::Other as InfoOther, key::Standard as InfoStandard, Key as InfoKey},
        Record as VcfRecord,
    };

    /// Record type for storing AN, AC_hom, AC_het, AC_hemi counts for chrX/chrY.
    #[derive(Default, Debug, PartialEq, Eq, Clone)]
    pub struct Counts {
        /// Total number of alleles.
        pub an: u32,
        /// Number of hom. alt. alleles.
        pub ac_hom: u32,
        /// Number of het. alt. alleles.
        pub ac_het: u32,
        /// Number of hemi. alt. alleles.
        pub ac_hemi: u32,
    }

    impl Counts {
        /// Create from the given VCF record.
        pub fn from_vcf(value: &VcfRecord) -> Self {
            tracing::trace!("@ {:?}", &value);

            let ac_hom_xx = value
                .info()
                .get(&InfoKey::Other(
                    InfoOther::from_str("nhomalt_female").expect("Invalid key: nhomalt_female?"),
                ))
                .unwrap_or_else(|| {
                    value
                        .info()
                        .get(&InfoKey::Other(
                            InfoOther::from_str("nhomalt_XX").expect("Invalid key: nhomalt_XX?"),
                        ))
                        .unwrap_or_default()
                });
            let ac_hom_xx = if let Some(ac_hom_xx) = ac_hom_xx {
                match ac_hom_xx {
                    noodles::vcf::record::info::field::Value::IntegerArray(ac_hom_xx) => {
                        ac_hom_xx[0].unwrap() as u32
                    }
                    _ => panic!("invalid type for nhomalt_female/nhomalt_XX"),
                }
            } else {
                0
            };

            let ac_xx = value
                .info()
                .get(&InfoKey::Other(
                    InfoOther::from_str("AC_female").expect("Invalid key: AC_female?"),
                ))
                .unwrap_or_else(|| {
                    value
                        .info()
                        .get(&InfoKey::Other(
                            InfoOther::from_str("AC_XX").expect("Invalid key: AC_XX?"),
                        ))
                        .unwrap_or_default()
                });
            let ac_xx = if let Some(ac_xx) = ac_xx {
                match ac_xx {
                    noodles::vcf::record::info::field::Value::IntegerArray(ac_het) => {
                        ac_het[0].unwrap() as u32
                    }
                    _ => panic!("invalid type for AC_female/AC_XX"),
                }
            } else {
                0
            };

            let ac_hom_xy = value
                .info()
                .get(&InfoKey::Other(
                    InfoOther::from_str("nhomalt_male").expect("Invalid key: nhomalt_male?"),
                ))
                .unwrap_or_else(|| {
                    value
                        .info()
                        .get(&InfoKey::Other(
                            InfoOther::from_str("nhomalt_XY").expect("Invalid key: nhomalt_XY?"),
                        ))
                        .unwrap_or_default()
                });
            let ac_hom_xy = if let Some(ac_hom_xy) = ac_hom_xy {
                match ac_hom_xy {
                    vcf::record::info::field::Value::IntegerArray(ac_hom_xx) => {
                        ac_hom_xx[0].unwrap() as u32
                    }
                    _ => panic!("invalid type for nhomalt_male/nhomalt_XY"),
                }
            } else {
                0
            };

            let ac_xy = value
                .info()
                .get(&InfoKey::Other(
                    InfoOther::from_str("AC_male").expect("Invalid key: AC_male/AC_XY?"),
                ))
                .unwrap_or_else(|| {
                    value
                        .info()
                        .get(&InfoKey::Other(
                            InfoOther::from_str("AC_XY").expect("Invalid key: AC_XY?"),
                        ))
                        .unwrap_or_default()
                });
            let ac_xy = if let Some(ac_xy) = ac_xy {
                match ac_xy {
                    vcf::record::info::field::Value::IntegerArray(ac_het) => {
                        ac_het[0].unwrap() as u32
                    }
                    _ => panic!("invalid type for AC_male/AC_XY"),
                }
            } else {
                0
            };

            let nonpar = value
                .info()
                .get(&InfoKey::Other(
                    InfoOther::from_str("nonpar").expect("Invalid key: nonpar?"),
                ))
                .is_some();

            let an = match value
                .info()
                .get(&InfoKey::Standard(InfoStandard::TotalAlleleCount))
                .unwrap()
                .unwrap()
            {
                vcf::record::info::field::Value::Integer(an) => *an as u32,
                _ => panic!("invalid type for AN"),
            };

            if nonpar {
                // not in PAR
                Counts {
                    ac_hom: ac_hom_xx,
                    ac_het: ac_xx - 2 * ac_hom_xx,
                    ac_hemi: ac_xy,
                    an,
                }
            } else {
                // is in PAR
                Counts {
                    ac_hom: ac_hom_xx + ac_hom_xy,
                    ac_het: ac_xx.saturating_sub(2 * ac_hom_xx + 2 * ac_hom_xy),
                    ac_hemi: 0,
                    an,
                }
            }
        }

        /// Read from buffer.
        pub fn from_buf(buf: &[u8]) -> Self {
            Self {
                an: LittleEndian::read_u32(&buf[0..4]),
                ac_hom: LittleEndian::read_u32(&buf[4..8]),
                ac_het: LittleEndian::read_u32(&buf[8..12]),
                ac_hemi: LittleEndian::read_u32(&buf[12..16]),
            }
        }

        /// Write to buffer.
        pub fn to_buf(&self, buf: &mut [u8]) {
            LittleEndian::write_u32(&mut buf[0..4], self.an);
            LittleEndian::write_u32(&mut buf[4..8], self.ac_hom);
            LittleEndian::write_u32(&mut buf[8..12], self.ac_het);
            LittleEndian::write_u32(&mut buf[12..16], self.ac_hemi);
        }
    }

    /// Record type for the "mitochondrial" column family.
    pub struct Record {
        /// Counts from gnomAD exomes.
        pub gnomad_exomes: Counts,
        /// Counts from gnomAD genomes.
        pub gnomad_genomes: Counts,
    }

    impl Record {
        /// Read from buffer.
        pub fn from_buf(buf: &[u8]) -> Self {
            Self {
                gnomad_exomes: Counts::from_buf(&buf[0..16]),
                gnomad_genomes: Counts::from_buf(&buf[16..32]),
            }
        }

        /// Write to buffer.
        pub fn to_buf(&self, buf: &mut [u8]) {
            self.gnomad_exomes.to_buf(&mut buf[0..16]);
            self.gnomad_genomes.to_buf(&mut buf[16..32]);
        }
    }
}

/// Code for reading VCF.
pub mod vcf {
    use std::ops::Deref;

    use noodles::vcf::{record::Chromosome, Record as VcfRecord};

    use super::chrom_name_to_key;

    /// A chromosomal change `CHROM-POS-REF-ALT`.
    #[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
    pub struct Var {
        pub chrom: String,
        pub pos: u32,
        pub reference: String,
        pub alternative: String,
    }

    impl Var {
        /// Create from the given VCF record.
        pub fn from_vcf(value: &VcfRecord) -> Self {
            let chrom = match value.chromosome() {
                Chromosome::Name(name) | Chromosome::Symbol(name) => name.to_owned(),
            };
            let pos: usize = value.position().into();
            let pos = pos as u32;
            let reference = value.reference_bases().to_string();
            let alternative = value.alternate_bases().deref()[0].to_string();

            Var {
                chrom,
                pos,
                reference,
                alternative,
            }
        }
    }

    impl From<Var> for Vec<u8> {
        fn from(val: Var) -> Self {
            let mut result = Vec::new();

            result.extend_from_slice(chrom_name_to_key(&val.chrom).as_bytes());
            result.extend_from_slice(&val.pos.to_be_bytes());
            result.extend_from_slice(val.reference.as_bytes());
            result.push(b'>');
            result.extend_from_slice(val.alternative.as_bytes());

            result
        }
    }
}

/// Convert chromosome to key in RocksDB.
pub fn chrom_name_to_key(name: &str) -> String {
    let chrom = if let Some(stripped) = name.strip_prefix("chr") {
        stripped
    } else {
        name
    };
    let chrom = if chrom == "M" {
        String::from("MT")
    } else if "XY".contains(chrom) {
        format!(" {chrom}")
    } else {
        String::from(chrom)
    };
    assert!(chrom.len() <= 2);
    assert!(!chrom.is_empty());
    if chrom.len() == 1 {
        format!("0{chrom}")
    } else {
        chrom
    }
}

/// Convert from RocksDB chromosome key part to chromosome name.
pub fn chrom_key_to_name(key: &str) -> String {
    assert!(key.len() == 2);
    if key.starts_with('0') || key.starts_with(' ') {
        key[1..].to_string()
    } else {
        key.to_string()
    }
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use noodles_util::variant::reader::Builder as VariantReaderBuilder;

    use super::mt::Counts as MtCounts;
    use super::vcf::Var as VcfVar;
    use super::{chrom_key_to_name, chrom_name_to_key};

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
        let path = "tests/data/db/create/seqvar_freqs/mt/helix.chrM.vcf";
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
        let path = "tests/data/db/create/seqvar_freqs/mt/helix.chrM.vcf";
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
