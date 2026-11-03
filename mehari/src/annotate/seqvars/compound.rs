use crate::annotate::cli::PhasingStrategy;
use crate::annotate::seqvars::consequence::logic::PADDING;
use crate::annotate::seqvars::consequence::span_end;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record::info::field::key::{END_POSITION, SV_LENGTHS};
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::info::field::value::Array;
use noodles::vcf::variant::record_buf::samples::sample::Value;
use std::collections::{BTreeMap, HashSet};

/// Phasing profile for a variant allele.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
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
    pub vcf_var: crate::annotate::seqvars::consequence::VcfVariant,
    pub record: RecordBuf,
    pub tx_accessions: HashSet<String>,
    pub min_tx_start: i32,
    pub max_tx_end: i32,
    pub phase_groups: Vec<PhaseGroup>,
}

/// A variant that is not evaluated together with others, e.g., a symbolic allele or a breakend.
///
/// Its transcripts can reach far behind it. So it stays known on its own, and does not hold back
/// the variants in between, whose annotations could fill the memory.
#[derive(Debug, Clone)]
struct Event {
    chromosome: String,
    start: i32,
    end: i32,
    tx_accessions: HashSet<String>,
    max_tx_end: i32,
    phase_groups: Vec<PhaseGroup>,
}

/// What grouping needs to know about a variant or an event.
struct Member<'a> {
    start: i64,
    end: i64,
    tx_accessions: &'a HashSet<String>,
    phase_groups: &'a [PhaseGroup],
}

impl<'a> From<&'a BufferedVariant> for Member<'a> {
    fn from(variant: &'a BufferedVariant) -> Self {
        let start = i64::from(variant.vcf_var.position);
        Self {
            start,
            end: start + variant.vcf_var.reference.len() as i64 - 1,
            tx_accessions: &variant.tx_accessions,
            phase_groups: &variant.phase_groups,
        }
    }
}

/// The `SVLEN` of the allele of `record`, if it has one value: `Number=.` of older headers keeps
/// the values of all alleles of a split record.
pub fn sv_length(record: &RecordBuf) -> Option<i32> {
    integer(record, SV_LENGTHS)
}

/// The value of the INFO field `key` of `record`, if it is one string.
pub(crate) fn string(record: &RecordBuf, key: &str) -> Option<String> {
    match record.info().get(key).flatten() {
        Some(InfoValue::String(value)) => Some(value.clone()),
        Some(InfoValue::Array(Array::String(values))) => match &values[..] {
            [Some(value)] => Some(value.clone()),
            _ => None,
        },
        _ => None,
    }
}

/// The value of the INFO field `key` of `record`, if it is one integer. Files without a header
/// line for `key` give a string below VCF 4.3.
pub(crate) fn integer(record: &RecordBuf, key: &str) -> Option<i32> {
    match record.info().get(key).flatten() {
        Some(InfoValue::Integer(value)) => Some(*value),
        Some(InfoValue::String(value)) => value.parse().ok(),
        Some(InfoValue::Array(Array::Integer(values))) => match values[..] {
            [Some(value)] => Some(value),
            _ => None,
        },
        _ => None,
    }
}

/// The last reference position of the event in `record` at `pos`.
///
/// The own `SVLEN` of an allele decides, as in VCF 4.5. `END` holds the longest allele of a record,
/// and stays so when `bcftools norm` splits the record, which also leaves an `SVLEN` with
/// `Number=.` unsplit. An `SVLEN` of 0 is the length difference of VCF 4.3 and earlier, e.g. of an
/// `<INV>`, so `END` decides. An insertion and a breakend lie at their position.
pub fn event_end(record: &RecordBuf, pos: i32) -> i32 {
    let alternate_bases = record.alternate_bases();
    let spans_bases = (alternate_bases.as_ref().first())
        .is_some_and(|alt| alt.starts_with('<') && !alt.starts_with("<INS"));
    if !spans_bases {
        return pos;
    }
    span_end(pos, sv_length(record), integer(record, END_POSITION))
}

/// Sliding window buffer that accumulates variants and groups them by transcript and phasing.
pub struct VariantBuffer {
    pub current_chrom: String,
    pub min_tx_start: i32,
    pub max_tx_end: i32,
    pub buffer: Vec<BufferedVariant>,
    strategy: PhasingStrategy,
    window_size: Option<u32>,
    events: Vec<Event>,
}

impl VariantBuffer {
    pub fn new(strategy: PhasingStrategy, window_size: Option<u32>) -> Self {
        Self {
            current_chrom: String::new(),
            min_tx_start: i32::MAX,
            max_tx_end: -1,
            buffer: Vec::new(),
            strategy,
            window_size,
            events: Vec::new(),
        }
    }

    /// Evaluate whether the buffer should be flushed based on the next variant's coordinates.
    pub fn should_flush(&self, next_chrom: &str, next_pos: i32) -> bool {
        if self.buffer.is_empty() {
            return false;
        }
        self.current_chrom != next_chrom
            || next_pos > self.max_tx_end + crate::annotate::seqvars::consequence::logic::PADDING
    }

    #[allow(clippy::too_many_arguments)]
    /// Push a variant into the sliding window and update the genomic boundaries.
    pub fn push(
        &mut self,
        vcf_var: crate::annotate::seqvars::consequence::VcfVariant,
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

    /// Push a variant that is not evaluated together with others. Its record only keeps its place
    /// in the output, so it takes no part in the sliding window.
    #[allow(clippy::too_many_arguments)]
    pub fn push_event(
        &mut self,
        vcf_var: crate::annotate::seqvars::consequence::VcfVariant,
        record: RecordBuf,
        tx_accessions: HashSet<String>,
        max_tx_end: i32,
        end: i32,
        sample_idx: usize,
        alt_allele_idx: usize,
    ) {
        if !tx_accessions.is_empty() {
            self.events.push(Event {
                chromosome: vcf_var.chromosome.clone(),
                start: vcf_var.position,
                end,
                tx_accessions,
                max_tx_end,
                phase_groups: self.extract_phasing(&record, sample_idx, alt_allele_idx),
            });
        }
        self.push(
            vcf_var,
            record,
            HashSet::new(),
            i32::MAX,
            -1,
            sample_idx,
            alt_allele_idx,
        );
    }

    /// Whether the allele `alt_allele_idx` of sample `sample_idx` in `record` lies on a
    /// haplotype that the phasing strategy evaluates together.
    pub fn is_phased(&self, record: &RecordBuf, sample_idx: usize, alt_allele_idx: usize) -> bool {
        (self
            .extract_phasing(record, sample_idx, alt_allele_idx)
            .iter())
        .any(|group| !matches!(group, PhaseGroup::Unphased))
    }

    /// Forget the events whose transcripts no variant from `next_pos` on can reach.
    pub fn forget_events(&mut self, next_chrom: &str, next_pos: i32) {
        self.events.retain(|event| {
            event.chromosome == next_chrom && next_pos <= event.max_tx_end + PADDING
        });
    }

    /// For each buffered variant, the transcripts on which it shares a haplotype with an event.
    /// Its annotations there do not cover all its phased partners.
    pub fn transcripts_with_events(&self) -> Vec<HashSet<String>> {
        let mut transcripts = vec![HashSet::new(); self.buffer.len()];
        if self.events.is_empty() {
            return transcripts;
        }

        // A large event has thousands of transcripts, and this runs for each flush, so only the
        // transcripts of the buffered variants count.
        let buffered = (self.buffer.iter())
            .flat_map(|variant| &variant.tx_accessions)
            .collect::<HashSet<_>>();
        let event_txs = (self.events.iter())
            .map(|event| {
                (buffered.iter())
                    .filter(|tx| event.tx_accessions.contains(**tx))
                    .map(|tx| (*tx).clone())
                    .collect::<HashSet<_>>()
            })
            .collect::<Vec<_>>();
        let events = self
            .events
            .iter()
            .zip(&event_txs)
            .map(|(event, txs)| Member {
                start: i64::from(event.start),
                end: i64::from(event.end),
                tx_accessions: txs,
                phase_groups: &event.phase_groups,
            });

        // In the order of the positions, with the buffer index of variants and none for events.
        let mut members = (self.buffer.iter().enumerate())
            .map(|(idx, variant)| (Some(idx), Member::from(variant)))
            .chain(events.map(|member| (None, member)))
            .collect::<Vec<_>>();
        members.sort_by_key(|(_, member)| member.start);
        let (variants, members): (Vec<_>, Vec<_>) = members.into_iter().unzip();

        for (tx, run) in self.runs_by_transcript(&members) {
            if run.iter().any(|&i| variants[i].is_none()) {
                for variant in run.iter().filter_map(|&i| variants[i]) {
                    transcripts[variant].insert(tx.to_string());
                }
            }
        }
        transcripts
    }

    /// Flush the buffer, return ordered variants and a list of grouped indices
    pub fn flush(&mut self) -> (Vec<BufferedVariant>, Vec<Vec<usize>>) {
        let all_variants = std::mem::take(&mut self.buffer);
        for b_var in &all_variants {
            tracing::trace!(
                "Evaluating buffered variant {}:{} against transcript bounds {}-{}",
                b_var.vcf_var.chromosome,
                b_var.vcf_var.position,
                b_var.min_tx_start,
                b_var.max_tx_end
            );
        }

        let members = all_variants.iter().map(Member::from).collect::<Vec<_>>();
        let mut compound_groups = self
            .runs_by_transcript(&members)
            .into_iter()
            .map(|(_, run)| run)
            .collect::<Vec<_>>();
        drop(members);

        // Sort the groups by variant index, so that the compound IDs follow the input order.
        compound_groups.sort_unstable();
        compound_groups.dedup();

        self.min_tx_start = i32::MAX;
        self.max_tx_end = -1;

        (all_variants, compound_groups)
    }

    /// The runs of at least two `members`, sorted by start, that lie on one transcript and one
    /// haplotype and within the window, with their transcript.
    fn runs_by_transcript<'a>(&self, members: &[Member<'a>]) -> Vec<(&'a str, Vec<usize>)> {
        let mut transcript_buckets: BTreeMap<&str, Vec<usize>> = BTreeMap::new();

        for (idx, member) in members.iter().enumerate() {
            for tx in member.tx_accessions {
                transcript_buckets.entry(tx).or_default().push(idx);
            }
        }

        let mut runs = Vec::new();

        for (tx, indices) in transcript_buckets {
            if indices.len() <= 1 {
                continue;
            }

            let mut phase_buckets: BTreeMap<PhaseGroup, Vec<usize>> = BTreeMap::new();
            let mut homozygous_indices = Vec::new();

            for &idx in &indices {
                for pg in members[idx].phase_groups {
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

                    for run in self.split_by_window(members, grouped_indices) {
                        if run.len() > 1 {
                            runs.push((tx, run));
                        }
                    }
                }
            }
        }

        runs
    }

    /// Split position-sorted indices into runs where each variant starts at most
    /// `window_size` bp after the end of the run so far. Without a window size, all
    /// indices form one run.
    fn split_by_window(&self, members: &[Member], indices: Vec<usize>) -> Vec<Vec<usize>> {
        let mut runs: Vec<Vec<usize>> = Vec::new();
        let mut run_end = 0;
        for idx in indices {
            let Member { start, end, .. } = members[idx];
            match runs.last_mut() {
                Some(run)
                    if self
                        .window_size
                        .is_none_or(|w| start - run_end <= i64::from(w)) =>
                {
                    run.push(idx);
                    run_end = run_end.max(end);
                }
                _ => {
                    runs.push(vec![idx]);
                    run_end = end;
                }
            }
        }
        runs
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

            if !is_phased && self.strategy == PhasingStrategy::Strict && !is_homozygous_alt {
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
            groups.sort_unstable();
            groups.dedup();
        }

        groups
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::annotate::seqvars::consequence::VcfVariant;
    use noodles::core::Position;
    use noodles::vcf::variant::record_buf::samples::Keys;
    use noodles::vcf::variant::record_buf::{AlternateBases, Samples};
    use rstest::rstest;

    /// Groups come out in input order.
    #[test]
    fn flush_orders_compound_groups() {
        let mut buffer = VariantBuffer::new(PhasingStrategy::Ignore, None);
        // Variant `i` lies on transcripts `tx{i-1}` and `tx{i}`, so neighbours form a group.
        for i in 0..8 {
            let vcf_var = VcfVariant {
                chromosome: "1".into(),
                position: 1000 + i,
                reference: "A".into(),
                alternative: "C".into(),
            };
            let tx_accessions = [format!("tx{}", i - 1), format!("tx{i}")].into();
            buffer.push(vcf_var, Default::default(), tx_accessions, 0, 2000, 0, 1);
        }

        let (_, compound_groups) = buffer.flush();

        let expected = (0..7).map(|i| vec![i, i + 1]).collect::<Vec<_>>();
        assert_eq!(compound_groups, expected);
    }

    /// A gap larger than the window between neighbours starts a new group.
    #[rstest]
    #[case(Some(10), vec![vec![0, 1, 2, 3], vec![5, 6]])]
    #[case(None, vec![(0..7).collect()])]
    fn flush_splits_groups_by_window(
        #[case] window_size: Option<u32>,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        let mut buffer = VariantBuffer::new(PhasingStrategy::Ignore, window_size);
        // The deletion at 1009 ends at 1012, so 1022 is still within the window.
        let variants = [
            (1000, "A", "C"),
            (1005, "A", "C"),
            (1009, "AAAA", "A"),
            (1022, "A", "C"),
            (1040, "A", "C"),
            (1060, "A", "C"),
            (1065, "A", "C"),
        ];
        for (position, reference, alternative) in variants {
            let vcf_var = VcfVariant {
                chromosome: "1".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            };
            buffer.push(
                vcf_var,
                Default::default(),
                ["tx".into()].into(),
                0,
                2000,
                0,
                1,
            );
        }

        let (_, compound_groups) = buffer.flush();

        assert_eq!(compound_groups, expected);
    }

    /// Push an SNV with genotype `0|1` in phase set `ps` on transcript `tx`.
    fn push_phased_snv(buffer: &mut VariantBuffer, tx: &str, pos: i32, ps: i32) {
        let keys: Keys = ["GT", "PS"].into_iter().map(String::from).collect();
        let samples = Samples::new(
            keys,
            vec![vec![
                Some(Value::String("0|1".into())),
                Some(Value::Integer(ps)),
            ]],
        );
        let record = RecordBuf::builder()
            .set_reference_sequence_name("17")
            .set_variant_start(Position::try_from(pos as usize).unwrap())
            .set_reference_bases("A")
            .set_alternate_bases(AlternateBases::from(vec![String::from("C")]))
            .set_samples(samples)
            .build();
        let vcf_var = VcfVariant {
            chromosome: "17".into(),
            position: pos,
            reference: "A".into(),
            alternative: "C".into(),
        };
        let tx_accessions = HashSet::from([tx.to_string()]);
        buffer.push(vcf_var, record, tx_accessions, 0, 1_000_000, 0, 1);
    }

    #[test]
    fn flush_returns_groups_in_input_order() {
        // 4 transcripts with 4 phase sets each, two variants per phase set:
        // 16 groups, so a hash-map order matches input order only by chance.
        let mut buffer = VariantBuffer::new(PhasingStrategy::Strict, None);
        let mut pos = 1000;
        for tx in ["tx1", "tx2", "tx3", "tx4"] {
            for ps in 1..=4 {
                for _ in 0..2 {
                    push_phased_snv(&mut buffer, tx, pos, ps);
                    pos += 10;
                }
            }
        }

        let (_, groups) = buffer.flush();

        let expected: Vec<Vec<usize>> = (0..16).map(|g| vec![2 * g, 2 * g + 1]).collect();
        assert_eq!(groups, expected);
    }

    /// Push an event at `pos..=end` with genotype `genotype` in phase set 1 on transcript `tx`.
    fn push_test_event(
        buffer: &mut VariantBuffer,
        tx: &str,
        genotype: &str,
        (pos, end): (i32, i32),
        max_tx_end: i32,
    ) {
        let keys: Keys = ["GT", "PS"].into_iter().map(String::from).collect();
        let samples = Samples::new(
            keys,
            vec![vec![
                Some(Value::String(genotype.into())),
                Some(Value::Integer(1)),
            ]],
        );
        let record = RecordBuf::builder()
            .set_reference_sequence_name("17")
            .set_variant_start(Position::try_from(pos as usize).unwrap())
            .set_reference_bases("N")
            .set_alternate_bases(AlternateBases::from(vec![String::from("<DEL>")]))
            .set_samples(samples)
            .build();
        let vcf_var = VcfVariant {
            chromosome: "17".into(),
            position: pos,
            reference: "N".into(),
            alternative: "<DEL>".into(),
        };
        let tx_accessions = HashSet::from([tx.to_string()]);
        buffer.push_event(vcf_var, record, tx_accessions, max_tx_end, end, 0, 1);
    }

    /// An event marks the variants on its transcripts and haplotype, and stays out of their group.
    #[test]
    fn event_marks_variants_and_joins_no_group() {
        let mut buffer = VariantBuffer::new(PhasingStrategy::Strict, None);
        push_phased_snv(&mut buffer, "tx1", 1000, 1);
        push_phased_snv(&mut buffer, "tx2", 1010, 1);
        push_test_event(&mut buffer, "tx1", "0|1", (1020, 1030), 5000);
        push_phased_snv(&mut buffer, "tx1", 1040, 1);
        push_phased_snv(&mut buffer, "tx1", 1050, 2);

        let tx1 = HashSet::from(["tx1".to_string()]);
        let none = HashSet::new();
        assert_eq!(
            buffer.transcripts_with_events(),
            [tx1.clone(), none.clone(), none.clone(), tx1, none]
        );
        let (_, groups) = buffer.flush();
        assert_eq!(groups, [vec![0, 3]]);
    }

    /// A run without an event stays unmarked, also while an event is known on its transcript.
    #[test]
    fn event_does_not_mark_runs_without_it() {
        let mut buffer = VariantBuffer::new(PhasingStrategy::Strict, None);
        push_phased_snv(&mut buffer, "tx", 1000, 2);
        push_phased_snv(&mut buffer, "tx", 1010, 2);
        push_test_event(&mut buffer, "tx", "0|1", (1020, 1030), 5000);

        let transcripts = buffer.transcripts_with_events();
        assert!(transcripts.iter().all(HashSet::is_empty), "{transcripts:?}");
    }

    /// Only an event within the window of a variant marks it. An event that spans the variant
    /// lies within.
    #[rstest]
    #[case::near((1005, 1006), true)]
    #[case::far((1020, 1021), false)]
    #[case::spans((900, 1500), true)]
    fn event_marks_variants_within_window(#[case] span: (i32, i32), #[case] marked: bool) {
        let mut buffer = VariantBuffer::new(PhasingStrategy::Strict, Some(10));
        push_phased_snv(&mut buffer, "tx", 1000, 1);
        push_test_event(&mut buffer, "tx", "0|1", span, 5000);

        assert_eq!(buffer.transcripts_with_events()[0].contains("tx"), marked);
    }

    /// An event marks the variants that it spans, but does not connect those that the window
    /// separates.
    #[test]
    fn event_does_not_bridge_window() {
        let mut buffer = VariantBuffer::new(PhasingStrategy::Strict, Some(10));
        push_phased_snv(&mut buffer, "tx", 1000, 1);
        push_test_event(&mut buffer, "tx", "0|1", (1005, 1100), 5000);
        push_phased_snv(&mut buffer, "tx", 1050, 1);

        let transcripts = buffer.transcripts_with_events();
        assert!(transcripts[0].contains("tx"));
        assert!(transcripts[2].contains("tx"));
        let (_, groups) = buffer.flush();
        assert!(groups.is_empty(), "{groups:?}");
    }

    /// An event stays known behind its own flush, until no variant can reach its transcripts.
    #[rstest]
    #[case::in_reach("17", 2000, true)]
    #[case::behind_transcripts("17", 5000 + PADDING + 1, false)]
    #[case::other_chromosome("18", 2000, false)]
    fn event_stays_known_while_in_reach(
        #[case] next_chrom: &str,
        #[case] next_pos: i32,
        #[case] known: bool,
    ) {
        let mut buffer = VariantBuffer::new(PhasingStrategy::Strict, None);
        push_test_event(&mut buffer, "tx", "0|1", (1000, 1100), 5000);
        buffer.flush();
        buffer.forget_events(next_chrom, next_pos);
        push_phased_snv(&mut buffer, "tx", next_pos, 1);

        assert_eq!(buffer.transcripts_with_events()[0].contains("tx"), known);
    }

    /// An event does not hold back the variants behind it, however far its transcripts reach.
    #[test]
    fn event_does_not_extend_window() {
        let mut buffer = VariantBuffer::new(PhasingStrategy::Strict, None);
        let vcf_var = VcfVariant {
            chromosome: "17".into(),
            position: 1000,
            reference: "A".into(),
            alternative: "C".into(),
        };
        let tx_accessions = HashSet::from(["tx".to_string()]);
        buffer.push(vcf_var, Default::default(), tx_accessions, 1000, 2000, 0, 1);
        push_test_event(&mut buffer, "tx", "0|1", (1500, 1600), 100_000_000);

        assert!(buffer.should_flush("17", 2000 + PADDING + 1));
    }

    #[rstest]
    #[case::svlen_negative("<DEL>", Some(InfoValue::Integer(-50)), None, 1050)]
    #[case::split_record_prefers_own_svlen(
        "<DUP>",
        Some(InfoValue::Array(Array::Integer(vec![Some(3000)]))),
        Some(9000),
        4000
    )]
    #[case::unsplit_svlen_falls_back_to_end(
        "<DUP>",
        Some(InfoValue::Array(Array::Integer(vec![Some(7000), Some(3000)]))),
        Some(9000),
        9000
    )]
    #[case::end_only("<DEL>", None, Some(1500), 1500)]
    #[case::end_before_position("<DEL>", None, Some(500), 1000)]
    #[case::neither("<DEL>", None, None, 1000)]
    #[case::insertion_length_is_no_span("<INS>", Some(InfoValue::Integer(300)), None, 1000)]
    #[case::mobile_element("<INS:ME:ALU>", Some(InfoValue::Integer(300)), None, 1000)]
    #[case::breakend("G]13:100]", Some(InfoValue::Integer(300)), Some(9000), 1000)]
    fn event_end_of_record(
        #[case] alternative: &str,
        #[case] svlen: Option<InfoValue>,
        #[case] end: Option<i32>,
        #[case] expected: i32,
    ) {
        let mut info = Vec::new();
        if let Some(svlen) = svlen {
            info.push((SV_LENGTHS.to_string(), Some(svlen)));
        }
        if let Some(end) = end {
            info.push((END_POSITION.to_string(), Some(InfoValue::Integer(end))));
        }
        let record = RecordBuf::builder()
            .set_alternate_bases(AlternateBases::from(vec![alternative.to_string()]))
            .set_info(info.into_iter().collect())
            .build();

        assert_eq!(event_end(&record, 1000), expected);
    }
}
