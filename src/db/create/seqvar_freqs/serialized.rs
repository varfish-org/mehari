//! Serialized population frequencies.

use byteorder::{ByteOrder, LittleEndian};

/// Record type for storing AN, AC_hom, AC_het counts for chrMT.
#[derive(Debug)]
pub struct MtCounts {
    /// Total number of alleles.
    pub an: u32,
    /// Number of homoplasmic alleles.
    pub ac_hom: u32,
    /// Number of heteroplasmic alleles.
    pub ac_het: u32,
}

impl MtCounts {
    /// Read from buffer.
    pub fn from_buf(buf: &[u8]) -> Self {
        Self {
            an: LittleEndian::read_u32(&buf[0..4]),
            ac_hom: LittleEndian::read_u32(&buf[5..8]),
            ac_het: LittleEndian::read_u32(&buf[9..12]),
        }
    }

    /// Write to buffer.
    pub fn to_buf(&self, buf: &mut [u8]) {
        LittleEndian::write_u32(&mut buf[0..4], self.an);
        LittleEndian::write_u32(&mut buf[5..8], self.an);
        LittleEndian::write_u32(&mut buf[9..12], self.an);
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
            helix_mtdb: MtCounts::from_buf(&buf[13..24]),
        }
    }

    /// Write to buffer.
    pub fn to_buf(&self, buf: &mut [u8]) {
        self.gnomad_mtdna.to_buf(&mut buf[0..12]);
        self.helix_mtdb.to_buf(&mut buf[13..24]);
    }
}

/// A chromosomal change `CHROM-POS-REF-ALT`.
///
/// The default value will be used as the "end" sentinel.
#[derive(Debug, Default, PartialEq, PartialOrd)]
pub struct VcfVar {
    pub chrom: String,
    pub pos: i32,
    pub reference: String,
    pub alternative: String,
}

impl VcfVar {
    pub fn pos_key(&self) -> (&String, i32) {
        (&self.chrom, self.pos)
    }
}
