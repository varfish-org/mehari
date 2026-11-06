//! The breakend records of one event, grouped by `EVENT` or `MATEID` while the input streams by,
//! see `MateGrouper`, each group with its event, see `BreakendEvent`.

use std::collections::VecDeque;
use std::sync::{Arc, OnceLock};

use itertools::Itertools;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record::info::field::key::{
    BREAKEND_EVENT_ID, EVENT_TYPE, MATE_BREAKEND_IDS,
};
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::info::field::value::Array;
use rustc_hash::{FxHashMap, FxHashSet};

use super::consequence::event::{BreakendEvent, DerivedVariant, EventMember};
use super::consequence::logic::ConsequencePredictor;
use super::consequence::{Annotations, Variant, VariantKind};
use super::{Classified, compound};

/// The event of a closed group, with the annotations of its derived variants, predicted once for
/// all of its records, see `ConsequencePredictor::predict_event`.
#[derive(Debug)]
struct GroupEvent {
    event: BreakendEvent,
    annotations: OnceLock<Result<Vec<Annotations>, String>>,
}

/// The event of a grouped record: the derived variant `derived` of the event of its group, which
/// its pairs describe, see `BreakendEvent::derived`.
#[derive(Debug, Clone)]
pub struct RecordEvent {
    group: Arc<GroupEvent>,
    derived: usize,
}

impl RecordEvent {
    /// The variant that the pairs of the record describe.
    fn variant(&self) -> &DerivedVariant {
        &self.group.event.derived[self.derived]
    }

    /// The last changed base of the event.
    pub fn end(&self) -> i32 {
        self.variant().sv.end
    }

    /// The annotations of the event, see `ConsequencePredictor::predict_event`, predicted once
    /// for all records of the group.
    pub fn annotations(&self, predictor: &ConsequencePredictor) -> anyhow::Result<&Annotations> {
        let annotations = (self.group.annotations.get_or_init(|| {
            (predictor.predict_event(&self.group.event)).map_err(|e| e.to_string())
        }))
        .as_ref()
        .map_err(|e| {
            anyhow::anyhow!(
                "failed to annotate {}, the event of breakends: {e}",
                self.variant().var
            )
        })?;
        Ok(&annotations[self.derived])
    }
}

/// A record that leaves `MateGrouper`, with its allele and its event, if any.
pub type Grouped = (RecordBuf, Classified, Option<RecordEvent>);

/// Holds the input records in order while a group of breakend records among them is open, so
/// that each record of a group gets the event of its pair, see `push`.
///
/// A breakend record with a mate and with `MATEID` or `EVENT` opens or joins a group: records
/// with one `EVENT` form a group, and so do records that `MATEID` links. A group is open until
/// the input passes its deadline, see `deadline`, or the contig changes, or the input ends. Then
/// its pairs give their events, see `BreakendEvent::new`, and its records whose mates never arrived
/// count as unresolved. A mate on another contig, e.g. of a translocation, is not waited for. A
/// closed group is gone: a record that names its `EVENT` or one of its IDs opens a new group.
#[derive(Default)]
pub struct MateGrouper {
    /// The held records in input order, each with its allele and its event, `None` while its
    /// group is open.
    #[allow(clippy::type_complexity)]
    held: VecDeque<(RecordBuf, Classified, Option<Option<RecordEvent>>)>,
    /// The input number of the first held record.
    first: usize,
    /// The open groups by number.
    groups: FxHashMap<usize, Group>,
    /// The number of the next group.
    next: usize,
    /// The open group of each `EVENT`.
    events: FxHashMap<String, usize>,
    /// The open group of each `ID` and `MATEID`.
    ids: FxHashMap<String, usize>,
    contig: String,
    /// No open group has an earlier deadline.
    earliest: i32,
    /// The number of records whose mates never arrived.
    pub unresolved: usize,
}

/// An open group of `MateGrouper`.
#[derive(Default)]
struct Group {
    /// The input numbers of the records of `members`.
    records: Vec<usize>,
    members: Vec<EventMember>,
    /// The `EVENT`s of the members, its keys in `MateGrouper::events`.
    events: Vec<String>,
    /// The `ID`s and `MATEID`s of the members, its keys in `MateGrouper::ids`.
    ids: Vec<String>,
    /// See `deadline`.
    deadline: i32,
}

impl MateGrouper {
    /// Takes the next record of the input with its allele `variant`, see `super::classify`, and
    /// gives the records that are ready, in input order.
    pub fn push(&mut self, record: RecordBuf, variant: Classified) -> Vec<Grouped> {
        let mut ready = Vec::new();
        if record.reference_sequence_name() != self.contig {
            ready = self.finish();
            self.contig = record.reference_sequence_name().to_string();
        }
        let position = record.variant_start().map_or(0, |pos| pos.get());
        let position = i32::try_from(position).unwrap_or(i32::MAX);
        self.close_passed(position);
        let number = self.first + self.held.len();
        let member = member(&record, &variant);
        // A record in no group is ready at once.
        self.held
            .push_back((record, variant, member.is_none().then_some(None)));
        if let Some((member, event)) = member {
            self.join(number, member, event);
            // A group whose mates lie in front of the record and are gone closes at once.
            self.close_passed(position);
        }
        ready.extend(self.drain());
        ready
    }

    /// Closes the open groups, and gives all held records, as at the end of the input.
    pub fn finish(&mut self) -> Vec<Grouped> {
        for number in self.groups.keys().copied().collect_vec() {
            self.close(number);
        }
        self.drain()
    }

    /// The held records in front of the first one of an open group.
    fn drain(&mut self) -> Vec<Grouped> {
        let mut ready = Vec::new();
        while let Some((_, _, Some(_))) = self.held.front() {
            if let Some((record, variant, Some(event))) = self.held.pop_front() {
                self.first += 1;
                ready.push((record, variant, event));
            }
        }
        ready
    }

    /// Adds `member`, the record of the input with the number `record`, with its `EVENT` to the
    /// open group that one of its keys names, or to a new group. The groups that it names merge.
    fn join(&mut self, record: usize, member: EventMember, event: Option<String>) {
        let named = (event.iter().filter_map(|event| self.events.get(event)))
            .chain((member.ids.iter().chain(&member.mate_ids)).filter_map(|id| self.ids.get(id)))
            .copied()
            .unique()
            .collect_vec();
        let number = named.first().copied().unwrap_or_else(|| {
            self.next += 1;
            self.next - 1
        });
        let mut group = self.groups.remove(&number).unwrap_or_default();
        for other in named.iter().skip(1) {
            if let Some(other) = self.groups.remove(other) {
                for key in &other.events {
                    self.events.insert(key.clone(), number);
                }
                for key in &other.ids {
                    self.ids.insert(key.clone(), number);
                }
                group.records.extend(other.records);
                group.members.extend(other.members);
                group.events.extend(other.events);
                group.ids.extend(other.ids);
            }
        }
        if let Some(key) = event {
            self.events.insert(key.clone(), number);
            group.events.push(key);
        }
        for key in member.ids.iter().chain(&member.mate_ids) {
            self.ids.insert(key.clone(), number);
            group.ids.push(key.clone());
        }
        group.records.push(record);
        group.members.push(member);
        group.deadline = deadline(&self.contig, &group.members);
        self.earliest = self.earliest.min(group.deadline);
        self.groups.insert(number, group);
    }

    /// Closes the groups whose deadline lies in front of `position`.
    fn close_passed(&mut self, position: i32) {
        if position <= self.earliest {
            return;
        }
        let passed = (self.groups.iter())
            .filter(|(_, group)| group.deadline < position)
            .map(|(&number, _)| number)
            .collect_vec();
        for number in passed {
            self.close(number);
        }
        self.earliest = (self.groups.values())
            .map(|group| group.deadline)
            .min()
            .unwrap_or(i32::MAX);
    }

    /// Closes the group `number`: the records of its pairs get their events, see
    /// `BreakendEvent::new`, and those whose mates never arrived count as unresolved.
    fn close(&mut self, number: usize) {
        let Some(group) = self.groups.remove(&number) else {
            return;
        };
        for key in &group.events {
            self.events.remove(key);
        }
        for key in &group.ids {
            self.ids.remove(key);
        }
        self.unresolved += missing_mates(&self.contig, &group.members);
        let group_event = Arc::new(GroupEvent {
            event: BreakendEvent::new(&self.contig, group.members),
            annotations: OnceLock::new(),
        });
        let mut events = vec![None; group.records.len()];
        for (derived, variant) in group_event.event.derived.iter().enumerate() {
            for &i in &variant.members {
                events[i] = Some(RecordEvent {
                    group: group_event.clone(),
                    derived,
                });
            }
        }
        // Records leave only once their group is closed, so those of the members are held.
        for (record, event) in group.records.iter().zip(events) {
            if let Some((_, _, slot)) = self.held.get_mut(record - self.first) {
                *slot = Some(event);
            }
        }
    }
}

/// The member of a group that `record` with the allele `variant` is, with its `EVENT`: a
/// breakend with a mate and with `MATEID` or `EVENT`.
fn member(record: &RecordBuf, variant: &Classified) -> Option<(EventMember, Option<String>)> {
    let event = compound::string(record, BREAKEND_EVENT_ID);
    let mate_ids = mate_ids(record);
    if event.is_none() && mate_ids.is_empty() {
        return None;
    }
    let (var, breakend) = match variant {
        Ok(Some(Variant {
            var,
            kind: VariantKind::Breakend(breakend),
        })) if breakend.mate.is_some() => (var.clone(), breakend.clone()),
        _ => return None,
    };
    let member = EventMember {
        var,
        breakend,
        ids: record.ids().as_ref().iter().cloned().collect(),
        mate_ids,
        event_type: compound::string(record, EVENT_TYPE)
            .and_then(|event_type| event_type.parse().ok()),
    };
    Some((member, event))
}

/// The `MATEID` values of `record`.
fn mate_ids(record: &RecordBuf) -> Vec<String> {
    match record.info().get(MATE_BREAKEND_IDS).flatten() {
        Some(InfoValue::String(ids)) => ids.split(',').map(String::from).collect(),
        Some(InfoValue::Array(Array::String(ids))) => ids.iter().flatten().cloned().collect(),
        _ => Vec::new(),
    }
}

/// The position after which a group of `members` on `chromosome` stops waiting for mates: the
/// largest mate position on `chromosome`, plus the largest `CIPOS` maximum among them, as an
/// imprecise mate can lie behind the position that its partner names, see
/// `BreakendEvent::new`. `i32::MIN` without a mate on `chromosome`.
fn deadline(chromosome: &str, members: &[EventMember]) -> i32 {
    let slack = (members.iter())
        .filter_map(|member| member.breakend.cipos)
        .map(|(_, max)| max.max(0))
        .max()
        .unwrap_or(0);
    (members.iter())
        .filter_map(|member| member.breakend.mate.as_ref())
        .filter(|mate| mate.chromosome == chromosome)
        .map(|mate| mate.position)
        .max()
        .map_or(i32::MIN, |position| position.saturating_add(slack))
}

/// The number of `members` whose mates on `chromosome` are not among them: those that `MATEID`
/// names, else those at the mate positions of their breakends.
fn missing_mates(chromosome: &str, members: &[EventMember]) -> usize {
    let ids = members
        .iter()
        .flat_map(|member| &member.ids)
        .collect::<FxHashSet<_>>();
    (members.iter())
        .filter(|member| match &member.breakend.mate {
            Some(mate) if mate.chromosome != chromosome => false,
            Some(mate) if member.mate_ids.is_empty() => !members
                .iter()
                .any(|other| other.var.position == mate.position),
            _ => !member.mate_ids.iter().all(|id| ids.contains(id)),
        })
        .count()
}

#[cfg(test)]
mod test {
    use super::{Grouped, MateGrouper};
    use crate::annotate::seqvars::classify;
    use itertools::Itertools;
    use noodles::vcf::variant::RecordBuf;

    /// A record of `contig` at `position` with `REF` `N`, `alt`, the ID `id` and the INFO fields
    /// `info`, e.g. `MATEID=b;IMPRECISE;CIPOS=-10,10`.
    fn record(contig: &str, position: usize, alt: &str, id: &str, info: &str) -> RecordBuf {
        use noodles::core::Position;
        use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
        use noodles::vcf::variant::record_buf::info::field::value::Array;
        use noodles::vcf::variant::record_buf::{AlternateBases, Ids};

        let info = (info.split(';').filter(|field| !field.is_empty()))
            .map(|field| match field.split_once('=') {
                None => (field.to_string(), Some(InfoValue::Flag)),
                Some(("CIPOS", values)) => {
                    let values = values.split(',').map(|value| value.parse().ok());
                    let values = InfoValue::Array(Array::Integer(values.collect()));
                    ("CIPOS".to_string(), Some(values))
                }
                Some((key, value)) => {
                    let value = InfoValue::Array(Array::String(vec![Some(value.to_string())]));
                    (key.to_string(), Some(value))
                }
            })
            .collect();
        RecordBuf::builder()
            .set_reference_sequence_name(contig)
            .set_variant_start(Position::try_from(position).unwrap())
            .set_ids(Ids::from_iter([id.to_string()]))
            .set_reference_bases("N")
            .set_alternate_bases(AlternateBases::from(vec![alt.to_string()]))
            .set_info(info)
            .build()
    }

    /// `MateGrouper::push` of `record` with its allele, see `classify`.
    fn push(grouper: &mut MateGrouper, record: RecordBuf) -> Vec<Grouped> {
        let variant = classify(&record, Default::default());
        grouper.push(record, variant)
    }

    /// A breakend record of `contig` at `position` with `alt`, the ID `id` and the `MATEID`
    /// `mate_id`.
    fn bnd_record(contig: &str, position: usize, alt: &str, id: &str, mate_id: &str) -> RecordBuf {
        record(contig, position, alt, id, &format!("MATEID={mate_id}"))
    }

    /// Only the records of a group whose mates never arrived count as unresolved: not those of
    /// a translocation, whose mate lies on another contig, nor a pair that gives no event.
    #[rstest::rstest]
    #[case::missing_mate(vec![("17", 100, "G[17:200[", "a", "gone")], 1)]
    #[case::missing_mate_at_contig_change(vec![("13", 100, "G[13:200[", "a", "gone"), ("17", 100, "G[17:200[", "b", "c"), ("17", 200, "]17:100]T", "c", "b")], 1)]
    #[case::translocation(vec![("13", 100, "]17:200]T", "a", "b"), ("17", 200, "G[13:100[", "b", "a")], 0)]
    #[case::head_to_head_pair(vec![("17", 100, "G]17:200]", "a", "b"), ("17", 200, "T]17:100]", "b", "a")], 0)]
    #[case::deletion(vec![("17", 100, "G[17:200[", "a", "b"), ("17", 200, "]17:100]T", "b", "a")], 0)]
    #[case::inserted_contig(vec![("17", 100, "G[<ctg1>:1[", "a", "c"), ("17", 101, "]<ctg1>:329]T", "b", "d")], 0)]
    fn unresolved_records(
        #[case] records: Vec<(&str, usize, &str, &str, &str)>,
        #[case] unresolved: usize,
    ) {
        let mut grouper = MateGrouper::default();
        let mut ready = 0;
        for (contig, position, alt, id, mate_id) in &records {
            ready += push(
                &mut grouper,
                bnd_record(contig, *position, alt, id, mate_id),
            )
            .len();
        }
        ready += grouper.finish().len();
        assert_eq!(ready, records.len());
        assert_eq!(grouper.unresolved, unresolved);
    }

    /// A group holds the records from its first one until the input passes its deadline, the
    /// largest mate position plus the largest `CIPOS` maximum, so a missing mate costs nothing
    /// beyond that: the records behind it flow. A record whose mate lies in front of it in a
    /// closed group does not wait.
    #[rstest::rstest]
    #[case::missing_mate(
        &[(100, "G[17:200[", "MATEID=gone"), (150, "C", ""), (200, "C", ""), (201, "C", ""), (300, "C", "")],
        &[0, 0, 0, 4, 1],
    )]
    #[case::imprecise_missing_mate(
        &[(100, "G[17:200[", "MATEID=gone;IMPRECISE;CIPOS=-10,10"), (201, "C", ""), (210, "C", ""), (211, "C", "")],
        &[0, 0, 0, 4],
    )]
    #[case::mate_in_closed_group(
        &[(100, "G[17:200[", "MATEID=b"), (200, "]17:100]T", "MATEID=a"), (300, "]17:100]T", "MATEID=a")],
        &[0, 0, 3],
    )]
    fn records_wait_until_deadline(
        #[case] records: &[(usize, &str, &str)],
        #[case] ready: &[usize],
    ) {
        let mut grouper = MateGrouper::default();
        let ids = ["a", "b", "c", "d", "e"];
        let pushed = (records.iter().zip(ids))
            .map(|(&(position, alt, info), id)| {
                push(&mut grouper, record("17", position, alt, id, info)).len()
            })
            .collect_vec();
        assert_eq!(pushed, ready);
        assert!(grouper.finish().is_empty());
        assert_eq!(grouper.unresolved, 1);
    }

    /// The event of each record as its POS, end and whether it is imprecise: each pair of a
    /// group gives its own, also a pair whose IDs an earlier pair had, and two pairs of one
    /// `EVENT`, one after the other or interleaved. Two records that `MATEID` links pair if the
    /// mate position of each lies within the `CIPOS` of the other, VCF 4.5 5.4.8. A breakend with
    /// two mates, split into two records with its ID, VCF 4.5 5.4.3, pairs with each mate on its
    /// own, and an `EVENTTYPE` that the pair does not have gives no event.
    #[rstest::rstest]
    #[case::reused_ids(
        &[(100, "G[17:200[", "a", "MATEID=b"), (200, "]17:100]T", "b", "MATEID=a"), (300, "G[17:400[", "a", "MATEID=b"), (400, "]17:300]T", "b", "MATEID=a")],
        &[Some((100, 199, false)), Some((100, 199, false)), Some((300, 399, false)), Some((300, 399, false))],
    )]
    #[case::event_of_pairs_in_turn(
        &[(100, "G[17:200[", "a", "EVENT=e"), (200, "]17:100]T", "b", "EVENT=e"), (300, "G[17:400[", "c", "EVENT=e"), (400, "]17:300]T", "d", "EVENT=e")],
        &[Some((100, 199, false)), Some((100, 199, false)), Some((300, 399, false)), Some((300, 399, false))],
    )]
    #[case::event_of_interleaved_pairs(
        &[(100, "G[17:300[", "a", "EVENT=e"), (200, "G[17:400[", "b", "EVENT=e"), (300, "]17:100]T", "c", "EVENT=e"), (400, "]17:200]T", "d", "EVENT=e")],
        &[Some((100, 299, false)), Some((200, 399, false)), Some((100, 299, false)), Some((200, 399, false))],
    )]
    #[case::imprecise_pair(
        &[(100, "G[17:205[", "a", "MATEID=b;IMPRECISE;CIPOS=-10,10"), (200, "]17:95]T", "b", "MATEID=a;IMPRECISE;CIPOS=-10,10")],
        &[Some((100, 199, true)), Some((100, 199, true))],
    )]
    #[case::imprecise_pair_beyond_cipos(
        &[(100, "G[17:215[", "a", "MATEID=b;IMPRECISE;CIPOS=-10,10"), (200, "]17:95]T", "b", "MATEID=a;IMPRECISE;CIPOS=-10,10")],
        &[None, None],
    )]
    #[case::imprecise_pair_without_mate_ids(
        &[(100, "G[17:205[", "a", "EVENT=e;IMPRECISE;CIPOS=-10,10"), (200, "]17:95]T", "b", "EVENT=e;IMPRECISE;CIPOS=-10,10")],
        &[None, None],
    )]
    #[case::multiple_mates(
        &[(100, "G[17:200[", "u", "MATEID=v"), (100, "G[13:300[", "u", "MATEID=z"), (200, "]17:100]T", "v", "MATEID=u")],
        &[Some((100, 199, false)), None, Some((100, 199, false))],
    )]
    #[case::multiple_mates_in_other_order(
        &[(100, "G[13:300[", "u", "MATEID=z"), (100, "G[17:200[", "u", "MATEID=v"), (200, "]17:100]T", "v", "MATEID=u")],
        &[None, Some((100, 199, false)), Some((100, 199, false))],
    )]
    #[case::double_minute(
        &[(100, "]17:200]G", "a", "MATEID=b;EVENTTYPE=DOUBLEMINUTE"), (200, "T[17:100[", "b", "MATEID=a;EVENTTYPE=DOUBLEMINUTE")],
        &[None, None],
    )]
    fn events_of_pairs(
        #[case] records: &[(usize, &str, &str, &str)],
        #[case] expected: &[Option<(i32, i32, bool)>],
    ) {
        let mut grouper = MateGrouper::default();
        let mut grouped = Vec::new();
        for &(position, alt, id, info) in records {
            grouped.extend(push(&mut grouper, record("17", position, alt, id, info)));
        }
        grouped.extend(grouper.finish());
        let events = (grouped.iter())
            .map(|(_, _, event)| {
                (event.as_ref()).map(|event| {
                    let derived = event.variant();
                    (derived.var.position, derived.sv.end, derived.sv.imprecise)
                })
            })
            .collect_vec();
        assert_eq!(events, expected);
    }
}
