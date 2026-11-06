//! The breakend records of one event, grouped by `EVENT` or `MATEID`, and the deletions, tandem
//! duplications and inversions that their pairs describe on one contig, see `derive_events`.

use std::collections::VecDeque;
use std::sync::{Arc, OnceLock};

use itertools::Itertools;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record::info::field::key::{BREAKEND_EVENT_ID, MATE_BREAKEND_IDS};
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::info::field::value::Array;
use rustc_hash::{FxHashMap, FxHashSet};

use super::consequence::logic::ConsequencePredictor;
use super::consequence::{Annotations, Breakend, StructuralVariant, SvKind, VcfVariant};
use super::{breakend, compound};

/// The event of a pair of breakend records, or of two pairs, see `derive_events`.
#[derive(Debug)]
pub struct GroupEvent {
    var: VcfVariant,
    sv: StructuralVariant,
    /// The bases inserted at the junctions of the event, in the order of
    /// `ConsequencePredictor::predict_sv_with_inserted`.
    inserted: [usize; 2],
    annotations: OnceLock<Result<Annotations, String>>,
}

impl GroupEvent {
    /// The last changed base of the event.
    pub fn end(&self) -> i32 {
        self.sv.end
    }

    /// The annotations of the event, see `ConsequencePredictor::predict_sv_with_inserted`,
    /// predicted once for all records of the event.
    pub fn annotations(&self, predictor: &ConsequencePredictor) -> anyhow::Result<&Annotations> {
        (self.annotations.get_or_init(|| {
            (predictor.predict_sv_with_inserted(&self.var, &self.sv, self.inserted))
                .map(Option::unwrap_or_default)
                .map_err(|e| e.to_string())
        }))
        .as_ref()
        .map_err(|e| {
            anyhow::anyhow!(
                "failed to annotate {}, the event of breakends: {e}",
                self.var
            )
        })
    }
}

/// A record that leaves `MateGrouper`, with its event, if any.
pub type Grouped = (RecordBuf, Option<Arc<GroupEvent>>);

/// Holds the input records in order while a group of breakend records among them is open, so
/// that each record of a group gets the event of its pair, see `push`.
///
/// A breakend record with a mate and with `MATEID` or `EVENT` opens or joins a group: records
/// with one `EVENT` form a group, and so do records that `MATEID` links. A group is open until
/// the input passes its deadline, see `deadline`, or the contig changes, or the input ends. Then
/// its pairs give their events, see `derive_events`, and its records whose mates never arrived
/// count as unresolved. A mate on another contig, e.g. of a translocation, is not waited for. A
/// closed group is gone: a record that names its `EVENT` or one of its IDs opens a new group.
#[derive(Default)]
pub struct MateGrouper {
    /// The held records in input order, each with its event, `None` while its group is open.
    held: VecDeque<(RecordBuf, Option<Option<Arc<GroupEvent>>>)>,
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
    members: Vec<Member>,
    /// The `EVENT`s of the members, its keys in `MateGrouper::events`.
    events: Vec<String>,
    /// The `ID`s and `MATEID`s of the members, its keys in `MateGrouper::ids`.
    ids: Vec<String>,
    /// See `deadline`.
    deadline: i32,
}

impl MateGrouper {
    /// Takes the next record of the input, and gives the records that are ready, in input order.
    pub fn push(&mut self, record: RecordBuf) -> Vec<Grouped> {
        let mut ready = Vec::new();
        if record.reference_sequence_name() != self.contig {
            ready = self.finish();
            self.contig = record.reference_sequence_name().to_string();
        }
        let position = record.variant_start().map_or(0, |pos| pos.get());
        let position = i32::try_from(position).unwrap_or(i32::MAX);
        self.close_passed(position);
        let member = member(&record, self.first + self.held.len());
        // A record in no group is ready at once.
        self.held
            .push_back((record, member.is_none().then_some(None)));
        if let Some((member, event)) = member {
            self.join(member, event);
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
        while let Some((_, Some(_))) = self.held.front() {
            if let Some((record, Some(event))) = self.held.pop_front() {
                self.first += 1;
                ready.push((record, event));
            }
        }
        ready
    }

    /// Adds `member` with its `EVENT` to the open group that one of its keys names, or to a new
    /// group. The groups that it names merge.
    fn join(&mut self, member: Member, event: Option<String>) {
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
    /// `derive_events`, and those whose mates never arrived count as unresolved.
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
        let mut events = vec![None; group.members.len()];
        for (members, event) in derive_events(&self.contig, &group.members) {
            let event = Arc::new(event);
            for i in members {
                events[i] = Some(event.clone());
            }
        }
        self.unresolved += missing_mates(&self.contig, &group.members);
        // Records leave only once their group is closed, so those of the members are held.
        for (member, event) in group.members.iter().zip(events) {
            if let Some((_, slot)) = self.held.get_mut(member.record - self.first) {
                *slot = Some(event);
            }
        }
    }
}

/// The member of a group that `record`, the record of the input with the number `number`, is,
/// with its `EVENT`: a breakend with a mate and with `MATEID` or `EVENT`.
fn member(record: &RecordBuf, number: usize) -> Option<(Member, Option<String>)> {
    let event = compound::string(record, BREAKEND_EVENT_ID);
    let mate_ids = mate_ids(record);
    if event.is_none() && mate_ids.is_empty() {
        return None;
    }
    let breakend = breakend(record).filter(|breakend| breakend.mate.is_some())?;
    let position = i32::try_from(record.variant_start().map_or(0, |pos| pos.get())).ok()?;
    let member = Member {
        record: number,
        position,
        reference: record.reference_bases().to_string(),
        breakend,
        ids: record.ids().as_ref().iter().cloned().collect(),
        mate_ids,
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
/// imprecise mate can lie behind the position that its partner names, see `pairs`. `i32::MIN`
/// without a mate on `chromosome`.
fn deadline(chromosome: &str, members: &[Member]) -> i32 {
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
fn missing_mates(chromosome: &str, members: &[Member]) -> usize {
    let ids = members
        .iter()
        .flat_map(|member| &member.ids)
        .collect::<FxHashSet<_>>();
    (members.iter())
        .filter(|member| match &member.breakend.mate {
            Some(mate) if mate.chromosome != chromosome => false,
            Some(mate) if member.mate_ids.is_empty() => {
                !members.iter().any(|other| other.position == mate.position)
            }
            _ => !member.mate_ids.iter().all(|id| ids.contains(id)),
        })
        .count()
}

/// A breakend record of a group.
#[derive(Debug, Clone)]
struct Member {
    /// The number of the record in the input.
    record: usize,
    position: i32,
    /// `REF`, whose first base is the padding base of an event that starts behind it.
    reference: String,
    breakend: Breakend,
    ids: Vec<String>,
    mate_ids: Vec<String>,
}

/// The pairs of `members` whose breakends join each other on `chromosome`, the one at the lower
/// position first. Two members that `MATEID` links both ways pair if the mate position of each
/// lies within the `CIPOS` of the other, as VCF 4.5 writes the mate position of an imprecise
/// breakend where the mate would lie if the own position were exact. Other members pair if the
/// mate position of each is the position of the other.
fn pairs(chromosome: &str, members: &[Member]) -> Vec<(usize, usize)> {
    let mate = |i: usize| {
        (members[i].breakend.mate.as_ref())
            .filter(|mate| mate.chromosome == chromosome)
            .map(|mate| mate.position)
    };
    let links =
        |i: usize, j: usize| (members[i].mate_ids.iter()).any(|id| members[j].ids.contains(id));
    let within = |position: i32, i: usize| {
        let (min, max) = members[i].breakend.cipos.unwrap_or((0, 0));
        let at = members[i].position;
        (at.saturating_add(min)..=at.saturating_add(max)).contains(&position)
    };
    let join = |i: usize, j: usize| match (mate(i), mate(j)) {
        (Some(mate_i), Some(mate_j)) if links(i, j) && links(j, i) => {
            within(mate_i, j) && within(mate_j, i)
        }
        (Some(mate_i), Some(mate_j)) => {
            mate_i == members[j].position && mate_j == members[i].position
        }
        _ => false,
    };
    let mut paired = vec![false; members.len()];
    let mut pairs = Vec::new();
    for i in 0..members.len() {
        let partner = (i + 1..members.len()).find(|&j| !paired[i] && !paired[j] && join(i, j));
        if let Some(j) = partner {
            (paired[i], paired[j]) = (true, true);
            let lower_first = members[i].position <= members[j].position;
            pairs.push(if lower_first { (i, j) } else { (j, i) });
        }
    }
    pairs
}

/// The events that the pairs of breakends among `members` on `chromosome` describe, see
/// `pairs`, each with the members of its pairs, as symbolic alleles with their padding base.
/// With POS a < b, the pair `t[b[` at a and `]a]t` at b is the deletion of a+1 to b-1, the pair
/// `]b]t` at a and `t[a[` at b the tandem duplication of a to b, and the pairs `t]b]` at a with
/// `t]a]` at b and `[b+1[t` at a+1 with `[a+1[t` at b+1 the inversion of a+1 to b, as in VCF
/// 4.5. Bases inserted at a junction do not change the event, only its reading frame, see
/// `GroupEvent::inserted`. Other pairs give none: a head-to-head or tail-to-tail pair alone,
/// other joins, or a deletion of no base.
fn derive_events(chromosome: &str, members: &[Member]) -> Vec<(Vec<usize>, GroupEvent)> {
    let positions = |(lo, hi): (usize, usize)| (members[lo].position, members[hi].position);
    // The bracket forms of the breakends of a pair, as whether the breakend keeps the left side
    // and whether its joined piece extends right of the mate.
    let form = |(lo, hi): (usize, usize)| {
        let form = |i: usize| {
            let breakend = &members[i].breakend;
            let extends_right = (breakend.mate.as_ref()).is_some_and(|mate| mate.extends_right);
            (breakend.keeps_left, extends_right)
        };
        (form(lo), form(hi))
    };
    const T_OPEN: (bool, bool) = (true, true); // t[p[
    const CLOSE_T: (bool, bool) = (false, false); // ]p]t
    const T_CLOSE: (bool, bool) = (true, false); // t]p]
    const OPEN_T: (bool, bool) = (false, true); // [p[t

    // The pairs of each event with its kind, POS and end, and the member at POS, whose base pads
    // the event.
    let mut found = Vec::new();
    let (mut heads, mut tails) = (Vec::new(), Vec::new());
    for pair in pairs(chromosome, members) {
        let (a, b) = positions(pair);
        match form(pair) {
            (T_OPEN, CLOSE_T) if b - a >= 2 => {
                found.push((vec![pair], SvKind::Deletion, a, b - 1, Some(pair.0)));
            }
            // The base in front of a is not in the records.
            (CLOSE_T, T_OPEN) if a >= 1 => {
                let kind = SvKind::Duplication { tandem: true };
                found.push((vec![pair], kind, a - 1, b, None));
            }
            (T_CLOSE, T_CLOSE) => heads.push(pair),
            (OPEN_T, OPEN_T) => tails.push(pair),
            _ => {}
        }
    }
    for head in heads {
        let (a, b) = positions(head);
        let next = a.checked_add(1).zip(b.checked_add(1));
        if let Some(i) = tails.iter().position(|&tail| Some(positions(tail)) == next) {
            let tail = tails.swap_remove(i);
            found.push((vec![head, tail], SvKind::Inversion, a, b, Some(head.0)));
        }
    }

    (found.into_iter())
        .map(|(pairs, kind, position, end, padding)| {
            let indices = pairs.iter().flat_map(|&(lo, hi)| [lo, hi]).collect_vec();
            let breakends = || indices.iter().map(|&i| &members[i].breakend);
            // The intervals around the start and the end are those of the breakends there. In
            // an imprecise event, a precise breakend is exact, and an imprecise one without
            // `CIPOS` leaves the side unknown.
            let imprecise = breakends().any(|breakend| breakend.imprecise);
            let side = |member_of_pair: fn((usize, usize)) -> usize| {
                let mut interval = None;
                for &pair in &pairs {
                    let breakend = &members[member_of_pair(pair)].breakend;
                    let (min, max) = match (breakend.imprecise, breakend.cipos) {
                        (true, cipos) => cipos?,
                        (false, Some(cipos)) => cipos,
                        (false, None) if imprecise => (0, 0),
                        (false, None) => continue,
                    };
                    interval = Some(match interval {
                        None => (min, max),
                        Some((min_a, max_a)) => (min.min(min_a), max.max(max_a)),
                    });
                }
                interval
            };
            let intervals = (side(|(lo, _)| lo), side(|(_, hi)| hi));
            let reference = padding
                .and_then(|i| members[i].reference.get(..1))
                .unwrap_or("N");
            let alternative = match kind {
                SvKind::Deletion => "<DEL>",
                SvKind::Duplication { .. } => "<DUP:TANDEM>",
                _ => "<INV>",
            };
            // The bases inserted at the junction of each pair, as its lower breakend has them.
            let inserted = [0, 1]
                .map(|n| (pairs.get(n)).map_or(0, |&(lo, _)| members[lo].breakend.inserted.len()));
            let event = GroupEvent {
                var: VcfVariant {
                    chromosome: chromosome.into(),
                    position,
                    reference: reference.into(),
                    alternative: alternative.into(),
                },
                sv: StructuralVariant {
                    kind,
                    end,
                    imprecise,
                    intervals,
                    junction: true,
                },
                inserted,
                annotations: OnceLock::new(),
            };
            (indices, event)
        })
        .collect()
}

#[cfg(test)]
mod test {
    use super::{MateGrouper, Member, derive_events};
    use crate::annotate::seqvars::consequence::{StructuralVariant, SvKind, VcfVariant};
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
    fn unresolved_records(
        #[case] records: Vec<(&str, usize, &str, &str, &str)>,
        #[case] unresolved: usize,
    ) {
        let mut grouper = MateGrouper::default();
        let mut ready = 0;
        for (contig, position, alt, id, mate_id) in &records {
            ready += (grouper.push(bnd_record(contig, *position, alt, id, mate_id))).len();
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
                (grouper.push(record("17", position, alt, id, info))).len()
            })
            .collect_vec();
        assert_eq!(pushed, ready);
        assert!(grouper.finish().is_empty());
        assert_eq!(grouper.unresolved, 1);
    }

    /// The event of each record as its POS, end and whether it is imprecise: each pair of a
    /// group gives its own, also a pair whose IDs an earlier pair had, and two pairs of one
    /// `EVENT`, one after the other or interleaved. Two records that `MATEID` links pair if the
    /// mate position of each lies within the `CIPOS` of the other, VCF 4.5 5.4.8.
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
    fn events_of_pairs(
        #[case] records: &[(usize, &str, &str, &str)],
        #[case] expected: &[Option<(i32, i32, bool)>],
    ) {
        let mut grouper = MateGrouper::default();
        let mut grouped = Vec::new();
        for &(position, alt, id, info) in records {
            grouped.extend(grouper.push(record("17", position, alt, id, info)));
        }
        grouped.extend(grouper.finish());
        let events = (grouped.iter())
            .map(|(_, event)| {
                (event.as_ref()).map(|event| (event.var.position, event.sv.end, event.sv.imprecise))
            })
            .collect_vec();
        assert_eq!(events, expected);
    }

    /// A member at `position` with `REF` `G` and the breakend `alt`, and `CIPOS` if given.
    fn member(position: i32, alt: &str, cipos: Option<(i32, i32)>) -> Member {
        let mut breakend: super::Breakend = alt.parse().unwrap();
        breakend.cipos = cipos;
        breakend.imprecise = cipos.is_some_and(|cipos| cipos != (0, 0));
        Member {
            record: 0,
            position,
            reference: "G".into(),
            breakend,
            ids: Vec::new(),
            mate_ids: Vec::new(),
        }
    }

    fn event(
        position: i32,
        reference: &str,
        alternative: &str,
        kind: SvKind,
        end: i32,
    ) -> Option<(VcfVariant, StructuralVariant)> {
        Some((
            VcfVariant {
                chromosome: "17".into(),
                position,
                reference: reference.into(),
                alternative: alternative.into(),
            },
            StructuralVariant {
                kind,
                end,
                imprecise: false,
                intervals: (None, None),
                junction: true,
            },
        ))
    }

    const DUP: SvKind = SvKind::Duplication { tandem: true };

    /// The members in input order and reversed give the same event. A member without partner
    /// gives none, while its group can.
    #[rstest::rstest]
    #[case::deletion(
        &[(100, "G[17:200["), (200, "]17:100]T")],
        event(100, "G", "<DEL>", SvKind::Deletion, 199),
    )]
    #[case::deletion_with_inserted_bases(
        &[(100, "GAC[17:200["), (200, "]17:100]CAT")],
        event(100, "G", "<DEL>", SvKind::Deletion, 199),
    )]
    #[case::deletion_of_one_base(
        &[(100, "G[17:102["), (102, "]17:100]T")],
        event(100, "G", "<DEL>", SvKind::Deletion, 101),
    )]
    #[case::adjacent_bases(&[(100, "G[17:101["), (101, "]17:100]T")], None)]
    #[case::tandem_duplication(
        &[(100, "]17:200]G"), (200, "T[17:100[")],
        event(99, "N", "<DUP:TANDEM>", DUP, 200),
    )]
    #[case::inversion(
        &[(100, "G]17:200]"), (101, "[17:201[T"), (200, "A]17:100]"), (201, "[17:101[C")],
        event(100, "G", "<INV>", SvKind::Inversion, 200),
    )]
    #[case::inversion_of_one_base(
        &[(100, "G]17:101]"), (101, "[17:102[T"), (101, "A]17:100]"), (102, "[17:101[C")],
        event(100, "G", "<INV>", SvKind::Inversion, 101),
    )]
    #[case::inversion_at_largest_position(
        &[(100, "G]17:2147483647]"), (101, "[17:2147483647[T"), (2147483647, "A]17:100]"), (2147483647, "[17:101[C")],
        None,
    )]
    #[case::mismatched_pair(&[(100, "G[17:200["), (200, "]17:101]T")], None)]
    #[case::mismatched_orientation(&[(100, "G[17:200["), (200, "T[17:100[")], None)]
    #[case::head_to_head_pair(&[(100, "G]17:200]"), (200, "A]17:100]")], None)]
    #[case::tail_to_tail_pair(&[(101, "[17:201[T"), (201, "[17:101[C")], None)]
    #[case::inversion_junctions_apart(
        &[(100, "G]17:200]"), (102, "[17:201[T"), (200, "A]17:100]"), (201, "[17:102[C")],
        None,
    )]
    #[case::mates_on_another_contig(&[(100, "G[13:200["), (200, "]13:100]T")], None)]
    #[case::member_without_partner(
        &[(100, "G[17:200["), (200, "]17:100]T"), (300, "]17:100]T")],
        event(100, "G", "<DEL>", SvKind::Deletion, 199),
    )]
    fn event_of_breakends(
        #[case] records: &[(i32, &str)],
        #[case] expected: Option<(VcfVariant, StructuralVariant)>,
    ) {
        let mut members = (records.iter())
            .map(|&(position, alt)| member(position, alt, None))
            .collect::<Vec<_>>();
        let events = |members: &[Member]| {
            (derive_events("17", members).into_iter())
                .map(|(_, event)| (event.var, event.sv))
                .collect_vec()
        };
        let expected = expected.into_iter().collect_vec();
        assert_eq!(events(&members), expected);
        members.reverse();
        assert_eq!(events(&members), expected);
    }

    /// The event is imprecise if a breakend is, with the intervals of the breakends at its
    /// start and at its end.
    #[test]
    fn imprecise_event_of_breakends() {
        let members = [
            member(100, "G[17:200[", Some((-5, 5))),
            member(200, "]17:100]T", None),
        ];
        let events = derive_events("17", &members);
        assert_eq!(
            (events.iter())
                .map(|(_, event)| (event.sv.imprecise, event.sv.intervals))
                .collect_vec(),
            [(true, (Some((-5, 5)), Some((0, 0))))]
        );
    }

    /// The bases inserted at each junction of an event, as the breakend at the lower position of
    /// its pair has them, with the records of the event.
    #[rstest::rstest]
    #[case::deletion(&[(100, "GAC[17:200["), (200, "]17:100]ACT")], [2, 0], &[0, 1])]
    #[case::inversion(
        &[(100, "GA]17:200]"), (101, "[17:201[CCT"), (200, "AT]17:100]"), (201, "[17:101[GGC")],
        [1, 2],
        &[0, 2, 1, 3],
    )]
    fn inserted_bases_of_events(
        #[case] records: &[(i32, &str)],
        #[case] inserted: [usize; 2],
        #[case] indices: &[usize],
    ) {
        let members = (records.iter())
            .map(|&(position, alt)| member(position, alt, None))
            .collect_vec();
        let events = (derive_events("17", &members).into_iter())
            .map(|(members, event)| (members, event.inserted))
            .collect_vec();
        assert_eq!(events, [(indices.to_vec(), inserted)]);
    }
}
