use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record_buf::samples::sample::Value;
use std::collections::{HashMap, HashSet};

/// Strategy used to evaluate compound variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum PhasingStrategy {
    /// Variants are only grouped if explicitly phased ('|') and sharing a Phase Set (PS).
    #[default]
    Strict,
    /// Variants in the same transcript on the same haplotype are grouped, ignoring missing or incompatible phasing metadata.
    Ignore,
}

/// Phasing profile for a variant allele.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum PhaseGroup {
    /// Phased allele containing an optional phase set and the haplotype index.
    Phased {
        phase_set: Option<i32>,
        haplotype_idx: usize,
    },
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
        header: &noodles::vcf::Header,
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

        let phase_groups = self.extract_phasing(header, &record, sample_idx, alt_allele_idx);

        self.buffer.push(BufferedVariant {
            vcf_var,
            record,
            tx_accessions,
            min_tx_start,
            max_tx_end,
            phase_groups,
        });
    }

    /// Flush the buffer, group eligible variants into compound variants, and return remaining single variants.
    pub fn flush(&mut self) -> Vec<Vec<BufferedVariant>> {
        let all_variants = std::mem::take(&mut self.buffer);
        let mut transcript_buckets: HashMap<String, Vec<BufferedVariant>> = HashMap::new();

        for b_var in &all_variants {
            for tx in &b_var.tx_accessions {
                transcript_buckets
                    .entry(tx.clone())
                    .or_default()
                    .push(b_var.clone());
            }
        }

        let mut final_groups: Vec<Vec<BufferedVariant>> = Vec::new();
        let mut consumed_singles: HashSet<String> = HashSet::new();
        let mut unique_mnv_signatures: HashSet<String> = HashSet::new();

        for (_, variants) in transcript_buckets {
            if variants.len() == 1 {
                continue;
            }

            let mut phase_buckets: HashMap<PhaseGroup, Vec<BufferedVariant>> = HashMap::new();
            for var in variants {
                for pg in &var.phase_groups {
                    phase_buckets
                        .entry(pg.clone())
                        .or_default()
                        .push(var.clone());
                }
            }

            for (phase_group, cis_variants) in phase_buckets {
                if cis_variants.len() > 1 && matches!(phase_group, PhaseGroup::Phased { .. }) {
                    let mut group_ids: Vec<String> = cis_variants
                        .iter()
                        .map(|v| {
                            format!(
                                "{}:{}:{}:{}",
                                v.vcf_var.chromosome,
                                v.vcf_var.position,
                                v.vcf_var.reference,
                                v.vcf_var.alternative
                            )
                        })
                        .collect();
                    group_ids.sort();
                    let signature = group_ids.join("|");

                    if !unique_mnv_signatures.contains(&signature) {
                        unique_mnv_signatures.insert(signature);

                        for id in group_ids {
                            consumed_singles.insert(id);
                        }

                        let mut unique_vars = Vec::new();
                        let mut seen = HashSet::new();
                        for v in cis_variants {
                            let id = format!(
                                "{}:{}:{}:{}",
                                v.vcf_var.chromosome,
                                v.vcf_var.position,
                                v.vcf_var.reference,
                                v.vcf_var.alternative
                            );
                            if !seen.contains(&id) {
                                seen.insert(id);
                                unique_vars.push(v);
                            }
                        }

                        unique_vars.sort_by_key(|v| v.vcf_var.position);
                        final_groups.push(unique_vars);
                    }
                }
            }
        }

        for b_var in all_variants {
            let id = format!(
                "{}:{}:{}:{}",
                b_var.vcf_var.chromosome,
                b_var.vcf_var.position,
                b_var.vcf_var.reference,
                b_var.vcf_var.alternative
            );
            if !consumed_singles.contains(&id) {
                final_groups.push(vec![b_var]);
                consumed_singles.insert(id);
            }
        }

        self.min_tx_start = i32::MAX;
        self.max_tx_end = -1;
        final_groups.sort_by_key(|g| g[0].vcf_var.position);
        final_groups
    }

    /// Extract Genotype and Phase Set to determine phasing of an allele.
    fn extract_phasing(
        &self,
        header: &noodles::vcf::Header,
        record: &RecordBuf,
        sample_idx: usize,
        alt_allele_idx: usize,
    ) -> Vec<PhaseGroup> {
        let samples = record.samples();

        let mut phase_set = None;
        if let Some(col) = samples.select("PS")
            && let Some(val) = col.get(sample_idx)
            && let Some(Value::Integer(i)) = val
        {
            phase_set = Some(*i);
        }

        let mut gt_str = None;
        if let Some(col) =
            samples.select(noodles::vcf::variant::record::samples::keys::key::GENOTYPE)
            && let Some(val) = col.get(sample_idx)
        {
            match val {
                Some(Value::String(s)) => {
                    gt_str = Some(s.to_owned());
                }
                Some(Value::Genotype(gt)) => {
                    gt_str = Some(crate::annotate::genotype_string(gt, header.file_format()));
                }
                _ => {}
            }
        }

        let mut groups = Vec::new();
        let alt_str = alt_allele_idx.to_string();

        if let Some(gt) = gt_str {
            let is_phased = gt.contains('|');

            if !is_phased && self.strategy == PhasingStrategy::Strict {
                return vec![PhaseGroup::Unphased];
            }

            for (hap_idx, allele_str) in gt.split(&['|', '/'][..]).enumerate() {
                if allele_str == alt_str {
                    if is_phased || self.strategy == PhasingStrategy::Ignore {
                        groups.push(PhaseGroup::Phased {
                            phase_set: if self.strategy == PhasingStrategy::Ignore {
                                Some(1)
                            } else {
                                phase_set
                            },
                            haplotype_idx: hap_idx,
                        });
                    } else {
                        groups.push(PhaseGroup::Unphased);
                    }
                }
            }
        }

        if groups.is_empty() {
            groups.push(PhaseGroup::Unphased);
        }

        groups
    }
}
