//! UCSC binning scheme.

// Offsets for UCSC binning.
//
// cf. http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
static BIN_OFFSETS: &[i32] = &[512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0];

const BIN_FIRST_SHIFT: i32 = 17;

const BIN_NEXT_SHIFT: i32 = 3;

// Compute UCSC bin from 0-based half-open interval.
pub fn bin_from_range(begin: i32, end: i32) -> Result<i32, anyhow::Error> {
    let mut begin_bin = begin >> BIN_FIRST_SHIFT;
    let mut end_bin = std::cmp::max(begin, end - 1) >> BIN_FIRST_SHIFT;

    for offset in BIN_OFFSETS {
        if begin_bin == end_bin {
            return Ok(offset + begin_bin);
        }
        begin_bin >>= BIN_NEXT_SHIFT;
        end_bin >>= BIN_NEXT_SHIFT;
    }

    anyhow::bail!(
        "begin {}, end {} out of range in bin_from_range (max is 512M",
        begin,
        end
    );
}
