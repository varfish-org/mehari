//! The breakend records of one event, grouped by `EVENT` or `MATEID` while the input streams by,
//! see `MateGrouper`, each group with its event, see `BreakendEvent`.

use std::collections::VecDeque;

use itertools::Itertools;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record::info::field::key::{
    BREAKEND_EVENT_ID, EVENT_TYPE, MATE_BREAKEND_IDS,
};
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::info::field::value::Array;
use rustc_hash::{FxHashMap, FxHashSet};

use super::consequence::event::{BreakendEvent, EventMember};
use super::consequence::{Variant, VariantKind};
use super::{Classified, compound};

/// A record that leaves `MateGrouper`, with its allele, and with the event of its group if it is
/// the last record of a group whose event has a derived variant, see `BreakendEvent::derived`.
pub type Grouped = (RecordBuf, Classified, Option<BreakendEvent>);

/// A place in the input: the rank of a contig, see `MateGrouper::ranks`, and a position on it.
type Locus = (usize, i32);

/// Holds the input records in order while the last record of a group of breakend records waits
/// for the event of its group, see `push`.
///
/// A breakend record with a mate and with `MATEID` or `EVENT` opens or joins a group: records
/// with one `EVENT` form a group, and so do records that `MATEID` links. A group is open until
/// the input passes its deadline, see `deadline`, or ends. Then its records whose mates never
/// arrived count as unresolved, and its latest record gets its event, see `BreakendEvent::new`,
/// if a pair of the group derives a variant. Only that record waits, with the records behind
/// it; the others leave at once. A closed group is gone: a record that names its `EVENT` or one
/// of its IDs opens a new group.
#[derive(Default)]
pub struct MateGrouper {
    /// The held records in input order, each with its allele and its event, `None` while it
    /// waits for its group.
    #[allow(clippy::type_complexity)]
    held: VecDeque<(RecordBuf, Classified, Option<Option<BreakendEvent>>)>,
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
    /// The rank of each contig: its index among the `##contig` lines of the header, then, for
    /// the contigs that the header lacks, the order in which the input reaches them.
    ranks: FxHashMap<String, usize>,
    /// No open group has an earlier deadline.
    earliest: Locus,
    /// The number of records whose mates never arrived.
    pub unresolved: usize,
}

/// An open group of `MateGrouper`.
#[derive(Default)]
struct Group {
    /// The members, with the variants that they derive.
    event: BreakendEvent,
    /// The `EVENT`s of the members, its keys in `MateGrouper::events`.
    events: Vec<String>,
    /// The `ID`s and `MATEID`s of the members, its keys in `MateGrouper::ids`.
    ids: Vec<String>,
    /// See `deadline`.
    deadline: Locus,
    /// The input number of the latest member while the group derives a variant: it waits for
    /// the event.
    carrier: Option<usize>,
}

impl MateGrouper {
    /// A grouper for an input whose header has the `##contig` lines `contigs`, in their order.
    pub fn new<S: Into<String>>(contigs: impl IntoIterator<Item = S>) -> Self {
        let ranks = (contigs.into_iter().enumerate())
            .map(|(rank, contig)| (contig.into(), rank))
            .collect();
        Self {
            ranks,
            ..Default::default()
        }
    }

    /// Takes the next record of the input with its allele `variant`, see `super::classify`, and
    /// gives the records that are ready, in input order.
    pub fn push(&mut self, record: RecordBuf, variant: Classified) -> Vec<Grouped> {
        let contig = record.reference_sequence_name();
        let rank = match self.ranks.get(contig) {
            Some(&rank) => rank,
            None => {
                let rank = self.ranks.len();
                self.ranks.insert(contig.to_string(), rank);
                rank
            }
        };
        let position = record.variant_start().map_or(0, |pos| pos.get());
        let locus = (rank, i32::try_from(position).unwrap_or(i32::MAX));
        self.close_passed(locus);
        let number = self.first + self.held.len();
        let member = member(&record, &variant);
        // A record in no group is ready at once, a member once it joins, see `join`.
        self.held
            .push_back((record, variant, member.is_none().then_some(None)));
        if let Some((member, event)) = member {
            self.join(number, member, event);
            // A group whose mates lie in front of the record and are gone closes at once.
            self.close_passed(locus);
        }
        self.drain()
    }

    /// Closes the open groups, and gives all held records, as at the end of the input.
    pub fn finish(&mut self) -> Vec<Grouped> {
        for number in self.groups.keys().copied().collect_vec() {
            self.close(number);
        }
        self.drain()
    }

    /// The held records in front of the first one that waits for its group.
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
    /// The record waits for the event if the group with it derives a variant, and the member
    /// that waited before leaves.
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
        let mut members = std::mem::take(&mut group.event.members);
        for other in named.iter().skip(1) {
            if let Some(other) = self.groups.remove(other) {
                for key in &other.events {
                    self.events.insert(key.clone(), number);
                }
                for key in &other.ids {
                    self.ids.insert(key.clone(), number);
                }
                members.extend(other.event.members);
                group.events.extend(other.events);
                group.ids.extend(other.ids);
                self.release(other.carrier);
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
        members.push(member);
        group.event = BreakendEvent::new(members);
        group.deadline = self.deadline(&group.event.members);
        self.earliest = self.earliest.min(group.deadline);
        self.release(group.carrier.take());
        if group.event.derived.is_empty() {
            self.release(Some(record));
        } else {
            group.carrier = Some(record);
        }
        self.groups.insert(number, group);
    }

    /// Lets the held record with the input number `record`, if any, leave without an event.
    fn release(&mut self, record: Option<usize>) {
        if let Some((_, _, slot)) = record.and_then(|record| self.held.get_mut(record - self.first))
        {
            *slot = Some(None);
        }
    }

    /// Closes the groups whose deadline lies in front of `locus`.
    fn close_passed(&mut self, locus: Locus) {
        if locus <= self.earliest {
            return;
        }
        let passed = (self.groups.iter())
            .filter(|(_, group)| group.deadline < locus)
            .map(|(&number, _)| number)
            .collect_vec();
        for number in passed {
            self.close(number);
        }
        self.earliest = (self.groups.values())
            .map(|group| group.deadline)
            .min()
            .unwrap_or((usize::MAX, i32::MAX));
    }

    /// Closes the group `number`: its records whose mates never arrived count as unresolved, and
    /// the member that waits gets the event.
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
        self.unresolved += missing_mates(&group.event.members);
        if let Some((_, _, slot)) =
            (group.carrier).and_then(|record| self.held.get_mut(record - self.first))
        {
            *slot = Some(Some(group.event));
        }
    }

    /// The locus after which a group of `members` stops waiting for mates: that of the farthest
    /// mate, with the largest `CIPOS` maximum among them added to its position, as an imprecise
    /// mate can lie behind the position that its partner names, see `BreakendEvent::new`. A mate
    /// on a contig that the header lacks and that the input has not reached ranks last, and one
    /// on an assembly contig is not waited for, see `is_assembly`. Without a mate to wait for,
    /// the first locus.
    fn deadline(&self, members: &[EventMember]) -> Locus {
        let slack = (members.iter())
            .filter_map(|member| member.breakend.cipos)
            .map(|(_, max)| max.max(0))
            .max()
            .unwrap_or(0);
        (members.iter())
            .filter_map(|member| member.breakend.mate.as_ref())
            .filter(|mate| !is_assembly(&mate.chromosome))
            .map(|mate| {
                let rank = self.ranks.get(&mate.chromosome).copied();
                (
                    rank.unwrap_or(usize::MAX),
                    mate.position.saturating_add(slack),
                )
            })
            .max()
            .unwrap_or((0, i32::MIN))
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

/// Whether `contig` is that of an assembly, e.g. `<ctg1>` of an inserted contig in VCF 4.5,
/// which has no records in the input.
fn is_assembly(contig: &str) -> bool {
    contig.starts_with('<')
}

/// The number of `members` whose mates are not among them: those that `MATEID` names, else those
/// at the mate positions of their breakends. A mate on an assembly contig is no record.
fn missing_mates(members: &[EventMember]) -> usize {
    let ids = members
        .iter()
        .flat_map(|member| &member.ids)
        .collect::<FxHashSet<_>>();
    (members.iter())
        .filter(|member| match &member.breakend.mate {
            Some(mate) if is_assembly(&mate.chromosome) => false,
            Some(mate) if member.mate_ids.is_empty() => !members.iter().any(|other| {
                other.var.chromosome == mate.chromosome && other.var.position == mate.position
            }),
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

    /// A record leaves at once, unless it is the latest member of a group that derives a
    /// variant: then it waits with the records behind it until the input passes the deadline of
    /// the group, the farthest mate plus the largest `CIPOS` maximum, and gets the event. The
    /// contigs rank in the order of the header, 2, 13 and 17, so a translocation completes on its
    /// last record. A group whose mate never arrives closes at its deadline, at once if the mate
    /// lies in front of the record, and at the end of the input if the mate lies on a contig
    /// that the header lacks, here 21; its records whose mates never arrived count, also one of
    /// a group whose event the latest member then gets. `unresolved` is that count after the
    /// input and after `finish`.
    #[rstest::rstest]
    #[case::deletion_pair(
        &[("17", 100, "G[17:200[", "a", "MATEID=b"), ("17", 150, "C", "s", ""), ("17", 200, "]17:100]T", "b", "MATEID=a"), ("17", 201, "C", "t", "")],
        &[1, 1, 0, 2],
        &["b"],
        (0, 0),
    )]
    #[case::imprecise_pair(
        &[("17", 100, "G[17:205[", "a", "MATEID=b;IMPRECISE;CIPOS=-10,10"), ("17", 200, "]17:95]T", "b", "MATEID=a;IMPRECISE;CIPOS=-10,10"), ("17", 215, "C", "s", ""), ("17", 216, "C", "t", "")],
        &[1, 0, 0, 3],
        &["b"],
        (0, 0),
    )]
    #[case::missing_mate(
        &[("17", 100, "G[17:200[", "a", "MATEID=gone"), ("17", 150, "C", "s", ""), ("17", 200, "C", "t", ""), ("17", 201, "C", "u", "")],
        &[1, 1, 1, 1],
        &[],
        (1, 1),
    )]
    #[case::imprecise_missing_mate(
        &[("17", 100, "G[17:200[", "a", "MATEID=gone;IMPRECISE;CIPOS=-10,10"), ("17", 201, "C", "s", ""), ("17", 210, "C", "t", ""), ("17", 211, "C", "u", "")],
        &[1, 1, 1, 1],
        &[],
        (1, 1),
    )]
    #[case::mate_in_closed_group(
        &[("17", 100, "G[17:200[", "a", "MATEID=b"), ("17", 200, "]17:100]T", "b", "MATEID=a"), ("17", 300, "]17:100]T", "c", "MATEID=a")],
        &[1, 0, 2],
        &["b"],
        (1, 1),
    )]
    #[case::dangling_member(
        &[("17", 100, "G[17:200[", "a", "EVENT=e"), ("17", 150, "G[17:400[", "c", "EVENT=e"), ("17", 200, "]17:100]T", "b", "EVENT=e"), ("17", 300, "C", "s", ""), ("17", 401, "C", "t", "")],
        &[1, 1, 0, 0, 3],
        &["b"],
        (1, 1),
    )]
    #[case::translocation(
        &[("2", 100, "G[13:200[", "a", "MATEID=b"), ("2", 300, "C", "s", ""), ("13", 200, "]2:100]T", "b", "MATEID=a"), ("13", 300, "C", "t", "")],
        &[1, 1, 1, 1],
        &[],
        (0, 0),
    )]
    #[case::event_across_contigs(
        &[("13", 100, "G[13:200[", "a", "EVENT=e"), ("13", 150, "G[17:500[", "c", "EVENT=e"), ("13", 200, "]13:100]T", "b", "EVENT=e"), ("17", 500, "]13:150]T", "d", "EVENT=e")],
        &[1, 1, 0, 1],
        &["d"],
        (0, 0),
    )]
    #[case::mate_on_contig_without_header_line(
        &[("17", 100, "G[21:200[", "a", "MATEID=gone"), ("17", 300, "C", "s", "")],
        &[1, 1],
        &[],
        (0, 1),
    )]
    #[case::mate_on_passed_contig(
        &[("17", 100, "G[13:200[", "a", "MATEID=gone")],
        &[1],
        &[],
        (1, 1),
    )]
    fn records_leave_at_once_but_last(
        #[case] records: &[(&str, usize, &str, &str, &str)],
        #[case] ready: &[usize],
        #[case] carriers: &[&str],
        #[case] unresolved: (usize, usize),
    ) {
        let mut grouper = MateGrouper::new(["2", "13", "17"]);
        let mut grouped = Vec::new();
        let mut pushed = Vec::new();
        for &(contig, position, alt, id, info) in records {
            let out = push(&mut grouper, record(contig, position, alt, id, info));
            pushed.push(out.len());
            grouped.extend(out);
        }
        assert_eq!(pushed, ready);
        let after_input = grouper.unresolved;
        grouped.extend(grouper.finish());
        assert_eq!((after_input, grouper.unresolved), unresolved);
        let with_event = (grouped.iter())
            .filter(|(_, _, event)| event.is_some())
            .map(|(record, _, _)| record.ids().as_ref().iter().join(";"))
            .collect_vec();
        assert_eq!(with_event, carriers);
        assert_eq!(grouped.len(), records.len());
    }

    /// The variants of the event that the records get, as the index of the record, POS, end and
    /// whether it is imprecise: the last record of a group gets them all. Each pair of a group gives its own variant, also
    /// a pair whose IDs an earlier pair had, and two pairs of one `EVENT`, one after the other or
    /// interleaved. Two records that `MATEID` links pair if the mate position of each lies
    /// within the `CIPOS` of the other, VCF 4.5 5.4.8. A breakend with two mates, split into two
    /// records with its ID, VCF 4.5 5.4.3, pairs with each mate on its own, and an `EVENTTYPE`
    /// that the pair does not have gives no variant.
    #[rstest::rstest]
    #[case::reused_ids(
        &[(100, "G[17:200[", "a", "MATEID=b"), (200, "]17:100]T", "b", "MATEID=a"), (300, "G[17:400[", "a", "MATEID=b"), (400, "]17:300]T", "b", "MATEID=a")],
        &[(1, 100, 199, false), (3, 300, 399, false)],
    )]
    #[case::event_of_pairs_in_turn(
        &[(100, "G[17:200[", "a", "EVENT=e"), (200, "]17:100]T", "b", "EVENT=e"), (300, "G[17:400[", "c", "EVENT=e"), (400, "]17:300]T", "d", "EVENT=e")],
        &[(1, 100, 199, false), (3, 300, 399, false)],
    )]
    #[case::event_of_interleaved_pairs(
        &[(100, "G[17:300[", "a", "EVENT=e"), (200, "G[17:400[", "b", "EVENT=e"), (300, "]17:100]T", "c", "EVENT=e"), (400, "]17:200]T", "d", "EVENT=e")],
        &[(3, 100, 299, false), (3, 200, 399, false)],
    )]
    #[case::imprecise_pair(
        &[(100, "G[17:205[", "a", "MATEID=b;IMPRECISE;CIPOS=-10,10"), (200, "]17:95]T", "b", "MATEID=a;IMPRECISE;CIPOS=-10,10")],
        &[(1, 100, 199, true)],
    )]
    #[case::imprecise_pair_beyond_cipos(
        &[(100, "G[17:215[", "a", "MATEID=b;IMPRECISE;CIPOS=-10,10"), (200, "]17:95]T", "b", "MATEID=a;IMPRECISE;CIPOS=-10,10")],
        &[],
    )]
    #[case::imprecise_pair_without_mate_ids(
        &[(100, "G[17:205[", "a", "EVENT=e;IMPRECISE;CIPOS=-10,10"), (200, "]17:95]T", "b", "EVENT=e;IMPRECISE;CIPOS=-10,10")],
        &[],
    )]
    #[case::multiple_mates(
        &[(100, "G[17:200[", "u", "MATEID=v"), (100, "G[13:300[", "u", "MATEID=z"), (200, "]17:100]T", "v", "MATEID=u")],
        &[(2, 100, 199, false)],
    )]
    #[case::multiple_mates_in_other_order(
        &[(100, "G[13:300[", "u", "MATEID=z"), (100, "G[17:200[", "u", "MATEID=v"), (200, "]17:100]T", "v", "MATEID=u")],
        &[(2, 100, 199, false)],
    )]
    #[case::double_minute(
        &[(100, "]17:200]G", "a", "MATEID=b;EVENTTYPE=DOUBLEMINUTE"), (200, "T[17:100[", "b", "MATEID=a;EVENTTYPE=DOUBLEMINUTE")],
        &[],
    )]
    fn events_of_pairs(
        #[case] records: &[(usize, &str, &str, &str)],
        #[case] expected: &[(usize, i32, i32, bool)],
    ) {
        let mut grouper = MateGrouper::default();
        let mut grouped = Vec::new();
        for &(position, alt, id, info) in records {
            grouped.extend(push(&mut grouper, record("17", position, alt, id, info)));
        }
        grouped.extend(grouper.finish());
        let events = (grouped.iter().enumerate())
            .flat_map(|(i, (_, _, event))| {
                (event.iter().flat_map(|event| &event.derived)).map(move |derived| {
                    let sv = &derived.sv;
                    (i, derived.var.position, sv.end, sv.imprecise)
                })
            })
            .collect_vec();
        assert_eq!(events, expected);
    }
}
