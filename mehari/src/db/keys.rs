//! VCF-style variant used to build lookup keys for the annotation databases.

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

impl From<Var> for annonars::common::keys::Var {
    fn from(var: Var) -> Self {
        annonars::common::keys::Var {
            chrom: var.chrom,
            pos: var.pos,
            reference: var.reference,
            alternative: var.alternative,
        }
    }
}

impl From<annonars::common::keys::Var> for Var {
    fn from(var: annonars::common::keys::Var) -> Self {
        Var {
            chrom: var.chrom,
            pos: var.pos,
            reference: var.reference,
            alternative: var.alternative,
        }
    }
}

impl From<&Var> for annonars::common::keys::Var {
    fn from(var: &Var) -> Self {
        annonars::common::keys::Var {
            chrom: var.chrom.clone(),
            pos: var.pos,
            reference: var.reference.clone(),
            alternative: var.alternative.clone(),
        }
    }
}

impl From<&annonars::common::keys::Var> for Var {
    fn from(var: &annonars::common::keys::Var) -> Self {
        Var {
            chrom: var.chrom.clone(),
            pos: var.pos,
            reference: var.reference.clone(),
            alternative: var.alternative.clone(),
        }
    }
}
