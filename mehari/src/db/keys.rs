//! Mehari-native assembly-agnostic key layouts for RocksDB.

#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord, Clone)]
pub struct Var {
    pub chrom: String,
    pub pos: i32,
    pub reference: String,
    pub alternative: String,
}

impl Var {
    pub fn new(chrom: String, pos: i32, reference: String, alternative: String) -> Self {
        Self {
            chrom,
            pos,
            reference,
            alternative,
        }
    }

    /// Helper to convert from a noodles RecordBuf directly.
    pub fn from_vcf_allele(value: &noodles::vcf::variant::RecordBuf, allele_no: usize) -> Self {
        let chrom = value.reference_sequence_name().to_string();
        let pos = value
            .variant_start()
            .expect("Variant start position is required")
            .get();

        Self {
            chrom,
            pos: pos as i32,
            reference: value.reference_bases().to_string(),
            alternative: value.alternate_bases().as_ref()[allele_no].to_string(),
        }
    }
}

/// Serialize Var into a binary key optimized for RocksDB byte sorting.
impl From<Var> for Vec<u8> {
    fn from(val: Var) -> Self {
        // Pre-allocate precisely to prevent vector reallocation overhead
        let estimated_capacity =
            val.chrom.len() + 1 + 4 + val.reference.len() + 1 + val.alternative.len();
        let mut result = Vec::with_capacity(estimated_capacity);

        // 1. Write raw chromosome string
        result.extend_from_slice(val.chrom.as_bytes());
        // 2. Delimit with a Null Byte. This guarantees that "chr1\0" sorts before "chr10\0"
        result.push(0x00);
        // 3. Write position as Big-Endian so it sorts numerically
        result.extend_from_slice(&val.pos.to_be_bytes());
        // 4. Append reference and alternative bases
        result.extend_from_slice(val.reference.as_bytes());
        result.push(b'>');
        result.extend_from_slice(val.alternative.as_bytes());

        result
    }
}

/// Deserialize from raw RocksDB bytes back into a human-readable Var struct.
impl From<&[u8]> for Var {
    fn from(value: &[u8]) -> Self {
        // Locate the null-byte separator bounding the chromosome name
        let null_idx = value
            .iter()
            .position(|&b| b == 0x00)
            .expect("Corrupted database key: missing chromosome null-terminator");

        let chrom = std::str::from_utf8(&value[0..null_idx])
            .expect("Invalid UTF-8 sequence in chromosome name")
            .to_string();

        // Extract the 4 big-endian position bytes immediately following the null byte
        let pos_start = null_idx + 1;
        let pos_end = pos_start + 4;
        let pos = i32::from_be_bytes(value[pos_start..pos_end].try_into().unwrap());

        // Parse allele structures bounded by the '>' symbol
        let alleles_buf = &value[pos_end..];
        let separator_idx = alleles_buf
            .iter()
            .position(|&b| b == b'>')
            .expect("Corrupted database key: missing allele '>' separator");

        let reference = std::str::from_utf8(&alleles_buf[0..separator_idx])
            .expect("Invalid UTF-8 sequence in reference allele")
            .to_string();
        let alternative = std::str::from_utf8(&alleles_buf[separator_idx + 1..])
            .expect("Invalid UTF-8 sequence in alternative allele")
            .to_string();

        Self {
            chrom,
            pos,
            reference,
            alternative,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip_arbitrary_scaffold() {
        let original = Var::new(
            "NW_021159990.1_unplaced_scaffold".to_string(),
            123456,
            "ATCG".to_string(),
            "G".to_string(),
        );

        let serialized: Vec<u8> = original.clone().into();
        let deserialized = Var::from(serialized.as_slice());

        assert_eq!(original, deserialized);
    }

    #[test]
    fn test_lexicographical_sorting_order() {
        let var_chrom1 = Var::new("chr1".into(), 100, "A".into(), "T".into());
        let var_chrom10 = Var::new("chr10".into(), 5, "A".into(), "T".into());

        let key1: Vec<u8> = var_chrom1.into();
        let key10: Vec<u8> = var_chrom10.into();

        // key1 must sort BEFORE key10 even though '5' < '100',
        // because "chr1\0" comes alphabetically before "chr10\0"
        assert!(key1 < key10);
    }
}
