//! The event of a set of breakends: the deletions, tandem duplications and inversions that their
//! pairs describe on one contig, see `BreakendEvent::new`.

use itertools::Itertools;

use super::{Breakend, StructuralVariant, SvKind, VcfVariant};

/// The breakends of one event, e.g. the records of one `EVENT`, with the variants that their
/// pairs describe.
#[derive(Debug, Clone, Default)]
pub struct BreakendEvent {
    pub members: Vec<EventMember>,
    /// The variants that pairs of `members` describe, see `BreakendEvent::new`.
    pub derived: Vec<DerivedVariant>,
}

/// A breakend of a `BreakendEvent`.
#[derive(Debug, Clone)]
pub struct EventMember {
    /// Its record, whose `REF` holds the padding base of an event that starts behind it.
    pub var: VcfVariant,
    pub breakend: Breakend,
    /// The `ID`s of its record.
    pub ids: Vec<String>,
    /// The `MATEID`s of its record.
    pub mate_ids: Vec<String>,
    /// The `EVENTTYPE` of its record, e.g. `DEL`.
    pub event_type: Option<EventType>,
}

/// `EVENTTYPE` of VCF 4.5: the type of the event of a breakend. `Other` keeps a type that the
/// spec does not reserve.
#[derive(Debug, PartialEq, Eq, Hash, Clone, strum::EnumString, strum::Display)]
#[strum(serialize_all = "UPPERCASE")]
pub enum EventType {
    Del,
    #[strum(serialize = "DEL:ME")]
    DelMe,
    Ins,
    #[strum(serialize = "INS:ME")]
    InsMe,
    Dup,
    #[strum(serialize = "DUP:TANDEM")]
    DupTandem,
    #[strum(serialize = "DUP:DISPERSED")]
    DupDispersed,
    Inv,
    Tra,
    #[strum(serialize = "TRA:BALANCED")]
    TraBalanced,
    #[strum(serialize = "TRA:UNBALANCED")]
    TraUnbalanced,
    Chromothripsis,
    Chromoplexy,
    Bfb,
    DoubleMinute,
    #[strum(default)]
    Other(String),
}

/// A variant that pairs of breakends of a `BreakendEvent` describe, as a symbolic allele with
/// its padding base.
#[derive(Debug, Clone)]
pub struct DerivedVariant {
    /// The members of its pairs, as indices into `BreakendEvent::members`.
    pub members: Vec<usize>,
    pub var: VcfVariant,
    pub sv: StructuralVariant,
    /// The bases inserted at the junctions of the variant, in the order of
    /// `ConsequencePredictor::sv_fusions`.
    pub inserted: [usize; 2],
}

impl BreakendEvent {
    /// The event of `members`, with the variants that their pairs on each contig describe, see
    /// `derive_events`. If a member of a pair has an `EVENTTYPE`, only a variant of that type is
    /// derived from it: `DEL`, `DUP:TANDEM` or `DUP`, or `INV`. Any other type,
    /// e.g. `TRA` or `DOUBLEMINUTE`, and members of different types give none. `EVENTTYPE` is
    /// per allele, so a group that a shared ID forms from two adjacencies checks each on its own.
    pub fn new(members: Vec<EventMember>) -> Self {
        let mut derived = (members.iter().map(|member| &member.var.chromosome))
            .unique()
            .flat_map(|chromosome| derive_events(chromosome, &members))
            .collect_vec();
        derived.retain(|derived| {
            let types = (derived.members.iter())
                .filter_map(|&i| members[i].event_type.as_ref())
                .unique()
                .collect_vec();
            let fits = match types[..] {
                [] => true,
                [event_type] => matches!(
                    (event_type, derived.sv.kind),
                    ("DEL", SvKind::Deletion)
                        | ("DUP:TANDEM" | "DUP", SvKind::Duplication { .. })
                        | ("INV", SvKind::Inversion)
                ),
                _ => false,
            };
            if !fits {
                tracing::debug!(
                    "{}: not derived from breakends with EVENTTYPE {}",
                    derived.var,
                    types.iter().join(",")
                );
            }
            fits
        });
        Self { members, derived }
    }
}

/// The pairs of `members` on `chromosome` whose breakends join each other there, the one at the
/// lower position first. Two members that `MATEID` links both ways pair if the mate position of
/// each lies within the `CIPOS` of the other, as VCF 4.5 writes the mate position of an imprecise
/// breakend where the mate would lie if the own position were exact. Other members pair if the
/// mate position of each is the position of the other.
fn pairs(chromosome: &str, members: &[EventMember]) -> Vec<(usize, usize)> {
    let mate = |i: usize| {
        (members[i].breakend.mate.as_ref())
            .filter(|mate| mate.chromosome == chromosome)
            .filter(|_| members[i].var.chromosome == chromosome)
            .map(|mate| mate.position)
    };
    let links =
        |i: usize, j: usize| (members[i].mate_ids.iter()).any(|id| members[j].ids.contains(id));
    let within = |position: i32, i: usize| {
        let (min, max) = members[i].breakend.cipos.unwrap_or((0, 0));
        let at = members[i].var.position;
        (at.saturating_add(min)..=at.saturating_add(max)).contains(&position)
    };
    let join = |i: usize, j: usize| match (mate(i), mate(j)) {
        (Some(mate_i), Some(mate_j)) if links(i, j) && links(j, i) => {
            within(mate_i, j) && within(mate_j, i)
        }
        (Some(mate_i), Some(mate_j)) => {
            mate_i == members[j].var.position && mate_j == members[i].var.position
        }
        _ => false,
    };
    let mut paired = vec![false; members.len()];
    let mut pairs = Vec::new();
    for i in 0..members.len() {
        let partner = (i + 1..members.len()).find(|&j| !paired[i] && !paired[j] && join(i, j));
        if let Some(j) = partner {
            (paired[i], paired[j]) = (true, true);
            let lower_first = members[i].var.position <= members[j].var.position;
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
/// `DerivedVariant::inserted`. Other pairs give none: a head-to-head or tail-to-tail pair alone,
/// other joins, or a deletion of no base.
fn derive_events(chromosome: &str, members: &[EventMember]) -> Vec<DerivedVariant> {
    let positions = |(lo, hi): (usize, usize)| (members[lo].var.position, members[hi].var.position);
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
                .and_then(|i| members[i].var.reference.get(..1))
                .unwrap_or("N");
            let alternative = match kind {
                SvKind::Deletion => "<DEL>",
                SvKind::Duplication { .. } => "<DUP:TANDEM>",
                _ => "<INV>",
            };
            // The bases inserted at the junction of each pair, as its lower breakend has them.
            let inserted = [0, 1]
                .map(|n| (pairs.get(n)).map_or(0, |&(lo, _)| members[lo].breakend.inserted.len()));
            DerivedVariant {
                members: indices,
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
            }
        })
        .collect()
}

#[cfg(test)]
mod test {
    use super::{BreakendEvent, EventMember, EventType, derive_events};
    use crate::annotate::seqvars::consequence::{StructuralVariant, SvKind, VcfVariant};
    use itertools::Itertools;

    /// A member at `position` with `REF` `G` and the breakend `alt`, and `CIPOS` if given.
    fn member(position: i32, alt: &str, cipos: Option<(i32, i32)>) -> EventMember {
        let mut breakend: super::Breakend = alt.parse().unwrap();
        breakend.cipos = cipos;
        breakend.imprecise = cipos.is_some_and(|cipos| cipos != (0, 0));
        EventMember {
            var: VcfVariant {
                chromosome: "17".into(),
                position,
                reference: "G".into(),
                alternative: alt.into(),
            },
            breakend,
            ids: Vec::new(),
            mate_ids: Vec::new(),
            event_type: None,
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
        let events = |members: &[EventMember]| {
            (derive_events("17", members).into_iter())
                .map(|derived| (derived.var, derived.sv))
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
                .map(|derived| (derived.sv.imprecise, derived.sv.intervals))
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
            .map(|derived| (derived.members, derived.inserted))
            .collect_vec();
        assert_eq!(events, [(indices.to_vec(), inserted)]);
    }

    /// `EVENTTYPE` keeps a derived variant of its type only: a pair that looks like a deletion
    /// or a tandem duplication can be part of a translocation or a double minute.
    #[rstest::rstest]
    #[case::deletion(&[(100, "G[17:200[", Some("DEL")), (200, "]17:100]T", Some("DEL"))], true)]
    #[case::deletion_as_translocation(&[(100, "G[17:200[", Some("TRA")), (200, "]17:100]T", Some("TRA"))], false)]
    #[case::different_types(&[(100, "G[17:200[", Some("DEL")), (200, "]17:100]T", Some("DUP"))], false)]
    #[case::tandem_duplication(&[(100, "]17:200]G", Some("DUP:TANDEM")), (200, "T[17:100[", None)], true)]
    #[case::duplication(&[(100, "]17:200]G", Some("DUP")), (200, "T[17:100[", Some("DUP"))], true)]
    #[case::double_minute(&[(100, "]17:200]G", Some("DOUBLEMINUTE")), (200, "T[17:100[", Some("DOUBLEMINUTE"))], false)]
    #[case::inversion(
        &[(100, "G]17:200]", Some("INV")), (101, "[17:201[T", Some("INV")), (200, "A]17:100]", Some("INV")), (201, "[17:101[C", Some("INV"))],
        true,
    )]
    #[case::per_pair(
        &[(100, "T[17:200[", Some("DEL")), (200, "]17:100]A", Some("DEL")), (300, "T[17:400[", Some("TRA")), (400, "]17:300]A", Some("TRA"))],
        true,
    )]
    fn event_type_restricts_event(
        #[case] records: &[(i32, &str, Option<&str>)],
        #[case] derived: bool,
    ) {
        let members = (records.iter())
            .map(|&(position, alt, event_type)| EventMember {
                event_type: event_type.map(|event_type| event_type.parse().unwrap()),
                ..member(position, alt, None)
            })
            .collect_vec();
        let event = BreakendEvent::new(members);
        assert_eq!(event.derived.len(), usize::from(derived));
    }

    /// `EVENTTYPE` reads and writes the values that VCF 4.5 reserves, and keeps any other.
    #[rstest::rstest]
    #[case("DEL", EventType::Del)]
    #[case("DEL:ME", EventType::DelMe)]
    #[case("INS", EventType::Ins)]
    #[case("INS:ME", EventType::InsMe)]
    #[case("DUP", EventType::Dup)]
    #[case("DUP:TANDEM", EventType::DupTandem)]
    #[case("DUP:DISPERSED", EventType::DupDispersed)]
    #[case("INV", EventType::Inv)]
    #[case("TRA", EventType::Tra)]
    #[case("TRA:BALANCED", EventType::TraBalanced)]
    #[case("TRA:UNBALANCED", EventType::TraUnbalanced)]
    #[case("CHROMOTHRIPSIS", EventType::Chromothripsis)]
    #[case("CHROMOPLEXY", EventType::Chromoplexy)]
    #[case("BFB", EventType::Bfb)]
    #[case("DOUBLEMINUTE", EventType::DoubleMinute)]
    #[case("del", EventType::Other("del".into()))]
    fn event_type_values(#[case] value: &str, #[case] event_type: EventType) {
        assert_eq!(value.parse(), Ok(event_type.clone()));
        assert_eq!(event_type.to_string(), value);
    }
}
