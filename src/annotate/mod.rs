//! Annotation of VCF files.

use noodles_vcf::header::FileFormat;
use noodles_vcf::variant::record::samples::series::value::genotype::Phasing;
use noodles_vcf::variant::record_buf::samples::sample::value::Genotype;

pub mod seqvars;
pub mod strucvars;

trait VcfRecord {
    fn get_format_value_by_sample_index(
        &self,
        key: &str,
        sample_idx: usize,
    ) -> Option<Option<&noodles_vcf::variant::record_buf::samples::sample::Value>>;
}

impl VcfRecord for noodles_vcf::variant::RecordBuf {
    fn get_format_value_by_sample_index(
        &self,
        key: &str,
        sample_idx: usize,
    ) -> Option<Option<&noodles_vcf::variant::record_buf::samples::sample::Value>> {
        self.samples()
            .get_index(sample_idx)
            .and_then(|v| v.get(key))
    }
}

const VCF_4_4: FileFormat = FileFormat::new(4, 4);

// noodles genotype writing is private: https://github.com/zaeleus/noodles/blob/master/noodles-vcf/src/io/writer/record/samples/sample/value/genotype.rs
// so we recreate that here
pub(crate) fn genotype_string(gt: &Genotype, vcf_version: FileFormat) -> String {
    let gt = gt
        .as_ref()
        .iter()
        .map(|allele| {
            let pos = allele.position();
            let pos = pos.map(|p| format!("{}", p)).unwrap_or(".".to_string());
            let phase = match allele.phasing() {
                Phasing::Phased => "|",
                Phasing::Unphased => "/",
            };
            (phase, pos)
        })
        .fold(String::new(), |mut a, (phase, pos)| {
            a.push_str(phase);
            a.push_str(&pos);
            a
        });

    if vcf_version < VCF_4_4 {
        // Stripping leading phasing information for versions < 4.4
        if gt.starts_with(['/', '|']) {
            gt.strip_prefix(['/', '|']).unwrap().to_string()
        } else {
            gt
        }
    } else {
        gt
    }
}
