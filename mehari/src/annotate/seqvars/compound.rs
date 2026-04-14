use crate::annotate::cli::PhasingStrategy;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record_buf::samples::sample::Value;
use std::collections::{HashMap, HashSet};

/// Phasing profile for a variant allele.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum PhaseGroup {
    /// Phased allele containing an optional phase set and the haplotype index.
    Phased {
        phase_set: Option<i32>,
        haplotype_idx: usize,
    },

    /// Allele present on all haplotypes.
    Homozygous,

    /// Unphased allele.
    Unphased,
}

/// Wrapper holding a VCF record, its mapped Mehari variant, and its genomic context.
#[derive(Debug, Clone)]
pub struct BufferedVariant {
    pub vcf_var: crate::annotate::seqvars::csq::VcfVariant,
    pub record: RecordBuf,
    pub tx_accessions: HashSet<String>,
    pub min_tx_start: i32,
    pub max_tx_end: i32,
    pub phase_groups: Vec<PhaseGroup>,
}

/// Sliding window buffer that accumulates variants and groups them by transcript and phasing.
pub struct VariantBuffer {
    pub current_chrom: String,
    pub min_tx_start: i32,
    pub max_tx_end: i32,
    pub buffer: Vec<BufferedVariant>,
    strategy: PhasingStrategy,
}

impl VariantBuffer {
    pub fn new(strategy: PhasingStrategy) -> Self {
        Self {
            current_chrom: String::new(),
            min_tx_start: i32::MAX,
            max_tx_end: -1,
            buffer: Vec::new(),
            strategy,
        }
    }

    /// Evaluate whether the buffer should be flushed based on the next variant's coordinates.
    pub fn should_flush(&self, next_chrom: &str, next_pos: i32) -> bool {
        if self.buffer.is_empty() {
            return false;
        }
        self.current_chrom != next_chrom
            || next_pos > self.max_tx_end + crate::annotate::seqvars::csq::PADDING
    }

    #[allow(clippy::too_many_arguments)]
    /// Push a variant into the sliding window and update the genomic boundaries.
    pub fn push(
        &mut self,
        vcf_var: crate::annotate::seqvars::csq::VcfVariant,
        record: RecordBuf,
        tx_accessions: HashSet<String>,
        min_tx_start: i32,
        max_tx_end: i32,
        sample_idx: usize,
        alt_allele_idx: usize,
    ) {
        if self.buffer.is_empty() {
            self.current_chrom = vcf_var.chromosome.clone();
            self.min_tx_start = min_tx_start;
            self.max_tx_end = max_tx_end;
        } else {
            self.min_tx_start = std::cmp::min(self.min_tx_start, min_tx_start);
            self.max_tx_end = std::cmp::max(self.max_tx_end, max_tx_end);
        }

        let phase_groups = self.extract_phasing(&record, sample_idx, alt_allele_idx);

        self.buffer.push(BufferedVariant {
            vcf_var,
            record,
            tx_accessions,
            min_tx_start,
            max_tx_end,
            phase_groups,
        });
    }

    /// Flush the buffer, return ordered variants and a list of grouped indices
    pub fn flush(&mut self) -> (Vec<BufferedVariant>, Vec<Vec<usize>>) {
        let all_variants = std::mem::take(&mut self.buffer);
        let mut transcript_buckets: HashMap<String, Vec<usize>> = HashMap::new();

        for (idx, b_var) in all_variants.iter().enumerate() {
            tracing::trace!(
                "Evaluating buffered variant {}:{} against transcript bounds {}-{}",
                b_var.vcf_var.chromosome,
                b_var.vcf_var.position,
                b_var.min_tx_start,
                b_var.max_tx_end
            );
            for tx in &b_var.tx_accessions {
                transcript_buckets.entry(tx.clone()).or_default().push(idx);
            }
        }

        let mut compound_groups = Vec::new();
        let mut unique_signatures = HashSet::new();

        for (_, indices) in transcript_buckets {
            if indices.len() <= 1 {
                continue;
            }

            let mut phase_buckets: HashMap<PhaseGroup, Vec<usize>> = HashMap::new();
            let mut homozygous_indices = Vec::new();

            for &idx in &indices {
                for pg in &all_variants[idx].phase_groups {
                    if *pg == PhaseGroup::Homozygous {
                        homozygous_indices.push(idx);
                    } else if matches!(pg, PhaseGroup::Phased { .. }) {
                        phase_buckets.entry(pg.clone()).or_default().push(idx);
                    }
                }
            }

            for grouped_indices in phase_buckets.values_mut() {
                grouped_indices.extend(homozygous_indices.iter().copied());
            }

            if phase_buckets.is_empty() && homozygous_indices.len() > 1 {
                phase_buckets.insert(
                    PhaseGroup::Phased {
                        phase_set: None,
                        haplotype_idx: 0,
                    },
                    homozygous_indices,
                );
            }

            for (phase_group, mut grouped_indices) in phase_buckets {
                if grouped_indices.len() > 1 && matches!(phase_group, PhaseGroup::Phased { .. }) {
                    grouped_indices.sort_unstable();
                    grouped_indices.dedup();

                    let signature = grouped_indices
                        .iter()
                        .map(|i| i.to_string())
                        .collect::<Vec<_>>()
                        .join("|");

                    if unique_signatures.insert(signature) {
                        compound_groups.push(grouped_indices);
                    }
                }
            }
        }

        self.min_tx_start = i32::MAX;
        self.max_tx_end = -1;

        (all_variants, compound_groups)
    }

    /// Extract Genotype and Phase Set to determine phasing of an allele.
    // Notice we don't even need the `header` argument anymore!
    fn extract_phasing(
        &self,
        record: &RecordBuf,
        sample_idx: usize,
        alt_allele_idx: usize,
    ) -> Vec<PhaseGroup> {
        use noodles::vcf::variant::record::samples::series::value::Genotype as GenotypeTrait;
        use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
        use noodles::vcf::variant::record_buf::samples::sample::value::Genotype as BufGenotype;
        use std::str::FromStr;

        let samples = record.samples();

        let mut phase_set = None;
        if let Some(col) = samples.select("PS")
            && let Some(val) = col.get(sample_idx)
            && let Some(Value::Integer(i)) = val
        {
            phase_set = Some(*i);
        }

        let mut parsed_gt: Vec<(Option<usize>, Phasing)> = Vec::new();

        if let Some(col) =
            samples.select(noodles::vcf::variant::record::samples::keys::key::GENOTYPE)
            && let Some(val) = col.get(sample_idx)
        {
            match val {
                Some(Value::String(s)) => {
                    if let Ok(gt) = BufGenotype::from_str(s) {
                        parsed_gt.extend((&gt).iter().filter_map(|res| res.ok()));
                    }
                }
                Some(Value::Genotype(gt)) => {
                    parsed_gt.extend(gt.iter().filter_map(|res| res.ok()));
                }
                _ => {}
            }
        }

        let mut groups = Vec::new();

        if !parsed_gt.is_empty() {
            let is_homozygous_alt = parsed_gt.len() > 1
                && parsed_gt
                    .iter()
                    .all(|(idx, _)| *idx == Some(alt_allele_idx));

            let is_phased = parsed_gt
                .iter()
                .any(|(_, phasing)| *phasing == Phasing::Phased);

            if !is_phased
                && self.strategy == PhasingStrategy::Strict
                && !is_homozygous_alt
            {
                return vec![PhaseGroup::Unphased];
            }

            for (hap_idx, (allele_idx_opt, _phasing)) in parsed_gt.into_iter().enumerate() {
                if allele_idx_opt == Some(alt_allele_idx) {
                    if is_phased {
                        groups.push(PhaseGroup::Phased {
                            phase_set,
                            haplotype_idx: hap_idx,
                        });
                    } else if self.strategy == PhasingStrategy::Ignore {
                        groups.push(PhaseGroup::Phased {
                            phase_set: Some(1),
                            haplotype_idx: hap_idx,
                        });
                    } else if is_homozygous_alt {
                        groups.push(PhaseGroup::Homozygous);
                    } else {
                        groups.push(PhaseGroup::Unphased);
                    }
                }
            }
        }

        if groups.is_empty() {
            if self.strategy == PhasingStrategy::Ignore {
                groups.push(PhaseGroup::Phased {
                    phase_set: Some(1),
                    haplotype_idx: 0,
                });
            } else {
                groups.push(PhaseGroup::Unphased);
            }
        } else {
            groups.sort_unstable_by_key(|g| format!("{:?}", g));
            groups.dedup();
        }

        groups
    }
}
