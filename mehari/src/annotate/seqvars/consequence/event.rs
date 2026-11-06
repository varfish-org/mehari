//! The event of a set of breakends: the deletions, inversions and duplications that their
//! adjacencies describe, from a walk over the segments between their breakpoints, see
//! `BreakendEvent::new`.

use itertools::Itertools;

use super::{Breakend, ConfidenceInterval, StructuralVariant, SvKind, VcfVariant};

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
    /// The members of the pairs that bound it, as indices into `BreakendEvent::members`.
    pub members: Vec<usize>,
    pub var: VcfVariant,
    pub sv: StructuralVariant,
    /// The bases inserted at the junctions of the variant, in the order of
    /// `ConsequencePredictor::sv_fusions`.
    pub inserted: [usize; 2],
}

impl BreakendEvent {
    /// The event of `members`, with the variants of the segments between the breakpoints of
    /// their pairs, see `Walk::derived`. If an end of a segment has more than one junction, e.g.
    /// of a breakend with two mates, the walk cannot tell which one a chain takes, and each pair
    /// on a contig gives its own variant, see `derive_events`. If a member of a variant has an
    /// `EVENTTYPE`, the variant is kept only if it fits: `DEL`, `DUP` and `INV` keep a variant of
    /// their kind, `DUP:TANDEM` a tandem dup, `DUP:DISPERSED` and `DOUBLEMINUTE` a dup elsewhere,
    /// and the complex types and translocations, e.g. `CHROMOTHRIPSIS` or `TRA`, every variant.
    /// Any other type, e.g. `INS`, and members of different types keep none. `EVENTTYPE` is per
    /// allele, so a group that a shared ID forms from two adjacencies checks each on its own.
    pub fn new(members: Vec<EventMember>) -> Self {
        let pairs = pairs(&members);
        let mut derived = match Walk::new(&members, &pairs) {
            Ok(walk) => walk.derived(&members, &pairs),
            Err(member) => {
                tracing::debug!(
                    "{}: a segment end with several junctions, so each pair gives its own variant",
                    members[member].var
                );
                (members.iter().map(|member| &member.var.chromosome))
                    .unique()
                    .flat_map(|chromosome| derive_events(chromosome, &members, &pairs))
                    .collect_vec()
            }
        };
        derived.retain(|derived| {
            let types = (derived.members.iter())
                .filter_map(|&i| members[i].event_type.as_ref())
                .unique()
                .collect_vec();
            let fits = match types[..] {
                [] => true,
                [event_type] => matches!(
                    (event_type, derived.sv.kind),
                    (EventType::Del, SvKind::Deletion)
                        | (EventType::Dup, SvKind::Duplication { .. })
                        | (EventType::DupTandem, SvKind::Duplication { tandem: true })
                        | (
                            EventType::DupDispersed | EventType::DoubleMinute,
                            SvKind::Duplication { tandem: false }
                        )
                        | (EventType::Inv, SvKind::Inversion)
                        | (
                            EventType::Chromothripsis
                                | EventType::Chromoplexy
                                | EventType::Bfb
                                | EventType::Tra
                                | EventType::TraBalanced
                                | EventType::TraUnbalanced,
                            _
                        )
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

/// The pairs of `members` whose breakends join each other, the one at the lower position first.
/// Two members that `MATEID` links both ways pair if the mate position of each lies on the
/// contig of the other, within its `CIPOS`, as VCF 4.5 writes the mate position of an imprecise
/// breakend where the mate would lie if the own position were exact. Other members pair if the
/// mate of each is the contig and position of the other.
fn pairs(members: &[EventMember]) -> Vec<(usize, usize)> {
    // The mate position of `i` if its mate lies on the contig of `j`.
    let mate = |i: usize, j: usize| {
        (members[i].breakend.mate.as_ref())
            .filter(|mate| mate.chromosome == members[j].var.chromosome)
            .map(|mate| mate.position)
    };
    let links =
        |i: usize, j: usize| (members[i].mate_ids.iter()).any(|id| members[j].ids.contains(id));
    let within = |position: i32, i: usize| {
        let (min, max) = members[i].breakend.cipos.unwrap_or((0, 0));
        let at = members[i].var.position;
        (at.saturating_add(min)..=at.saturating_add(max)).contains(&position)
    };
    let join = |i: usize, j: usize| match (mate(i, j), mate(j, i)) {
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

/// The breakpoint graph of the pairs of an event with the chains through it, after VCF 4.5: "an
/// arbitrary rearrangement can be summarized as a set of novel adjacencies". The breakpoints cut
/// each contig into segments. A segment has a left and a right end, `2 * segment` and
/// `2 * segment + 1`, which the reference joins to the ends of its neighbours, and a pair joins
/// the ends at which its breakends sit by a junction.
struct Walk<'a> {
    /// The segments of the contigs with breakpoints, in the order of the members, each contig
    /// from its start to its end.
    segments: Vec<Segment<'a>>,
    /// The junction at each end of a segment, if any.
    junctions: Vec<Option<Junction>>,
    /// The chains, see `Walk::new`.
    chains: Vec<Chain>,
}

/// A segment of a `Walk`: the bases between two gaps of a contig, each gap g between the bases
/// g and g + 1, or between a telomere and a gap.
#[derive(Clone, Copy)]
struct Segment<'a> {
    contig: &'a str,
    /// The gap in front of the segment, `None` at the start of the contig.
    front: Option<i32>,
    /// The gap behind the segment, `None` at the end of the contig.
    back: Option<i32>,
}

/// A junction at an end of a segment of a `Walk`.
#[derive(Clone, Copy)]
struct Junction {
    /// Its pair, as an index into the pairs.
    pair: usize,
    /// The breakend at the end, as an index into the members.
    member: usize,
    /// The end at the other side.
    other: usize,
    /// Whether its pair is that of a circle, with the `EVENTTYPE` `DOUBLEMINUTE`.
    circular: bool,
}

/// A chain of a `Walk`: the segments in its order, see `Visit`, and how it ends.
struct Chain {
    visits: Vec<Visit>,
    end: ChainEnd,
}

/// A segment in a `Chain`.
#[derive(Clone, Copy)]
struct Visit {
    segment: usize,
    /// Whether the chain runs through it from its left end to its right end.
    forward: bool,
    /// Whether the chain enters it by a junction, else by the reference or at a telomere.
    by_junction: bool,
}

/// How a `Chain` ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChainEnd {
    /// At a telomere, as it starts at one.
    Telomere,
    /// Back at the junction that starts it, as a circle.
    Cycle,
    /// Somewhere unknown: at a breakpoint whose other side has a junction of its own, or, for a
    /// chain that starts at a junction, at a telomere.
    Open,
}

impl<'a> Walk<'a> {
    /// The graph of the pairs of `members` with its chains, or the member at an end with a
    /// second junction, for which a chain has two ways to go.
    ///
    /// A breakend that keeps the left side, `t[p[` or `t]p]`, cuts the gap behind its position
    /// and sits on the right end of the segment in front of it; one that keeps the right side
    /// cuts the gap in front of its position and sits on the left end of the segment behind
    /// it. A pair is a junction if each breakend joins the side that the other keeps.
    ///
    /// A chain starts at each telomere that no chain reached, and runs through a segment from
    /// the end where it enters to the other end. There it takes the junction if no chain took
    /// it yet, else it goes on along the reference into the next segment of the contig, until a
    /// telomere. The reference leads into no end with a junction of its own, as the junction
    /// replaces it there, see VCF 4.5 "Explicit partners", except into a segment that a junction
    /// not taken yet joins to itself, the first copy of a tandem dup: there the chain stops,
    /// open. Chains from telomeres leave out the junctions of a double minute, as its circle is
    /// no part of the chromosome. Then each junction that no chain took starts a chain at its
    /// other end: a circle if the chain comes back to the junction, else open.
    fn new(members: &'a [EventMember], pairs: &[(usize, usize)]) -> Result<Self, usize> {
        let joins = |i: usize, j: usize| {
            (members[i].breakend.mate.as_ref())
                .is_some_and(|mate| mate.extends_right != members[j].breakend.keeps_left)
        };
        let junctions = (pairs.iter().enumerate())
            .filter(|&(_, &(i, j))| joins(i, j) && joins(j, i))
            .collect_vec();
        let gap = |i: usize| {
            let EventMember { var, breakend, .. } = &members[i];
            if breakend.keeps_left {
                var.position
            } else {
                var.position.saturating_sub(1)
            }
        };
        // The segments, and the end at which each breakend of a junction sits.
        let mut segments = Vec::new();
        let mut ends = vec![0; members.len()];
        for contig in members
            .iter()
            .map(|member| member.var.chromosome.as_str())
            .unique()
        {
            let breakends = (junctions.iter())
                .flat_map(|&(_, &(i, j))| [i, j])
                .filter(|&i| members[i].var.chromosome == contig)
                .collect_vec();
            let gaps = breakends
                .iter()
                .map(|&i| gap(i))
                .sorted()
                .dedup()
                .collect_vec();
            if gaps.is_empty() {
                continue;
            }
            let first = segments.len();
            let fronts = std::iter::once(None).chain(gaps.iter().copied().map(Some));
            let backs = gaps.iter().copied().map(Some).chain([None]);
            segments.extend((fronts.zip(backs)).map(|(front, back)| Segment {
                contig,
                front,
                back,
            }));
            for i in breakends {
                let segment = first + gaps.binary_search(&gap(i)).unwrap_or_else(|k| k);
                ends[i] = if members[i].breakend.keeps_left {
                    2 * segment + 1
                } else {
                    2 * segment + 2
                };
            }
        }
        let mut walk = Self {
            junctions: vec![None; 2 * segments.len()],
            segments,
            chains: Vec::new(),
        };
        for (pair, &(i, j)) in junctions {
            let circular = [i, j]
                .iter()
                .any(|&k| members[k].event_type == Some(EventType::DoubleMinute));
            for (member, other) in [(i, j), (j, i)] {
                let slot = &mut walk.junctions[ends[member]];
                if slot.is_some() {
                    return Err(member);
                }
                *slot = Some(Junction {
                    pair,
                    member,
                    other: ends[other],
                    circular,
                });
            }
        }

        let mut taken = vec![false; pairs.len()];
        let mut reached = vec![false; walk.junctions.len()];
        let mut chains = Vec::new();
        for start in 0..walk.junctions.len() {
            if walk.neighbour(start).is_none() && !reached[start] {
                let (chain, last) = walk.chain(start, None, &mut taken);
                if chain.end == ChainEnd::Telomere {
                    reached[last] = true;
                }
                chains.push(chain);
            }
        }
        for end in 0..walk.junctions.len() {
            if let Some(junction) = walk.junctions[end]
                && !taken[junction.pair]
            {
                taken[junction.pair] = true;
                let (chain, _) = walk.chain(junction.other, Some(junction.pair), &mut taken);
                chains.push(chain);
            }
        }
        walk.chains = chains;
        Ok(walk)
    }

    /// The end across the gap from `end`, `None` at a telomere.
    fn neighbour(&self, end: usize) -> Option<usize> {
        let segment = &self.segments[end / 2];
        if end % 2 == 1 {
            segment.back.map(|_| end + 1)
        } else {
            segment.front.map(|_| end - 1)
        }
    }

    /// The chain that enters a segment at `entry`, from a telomere, or by the junction `first`
    /// that starts it, with the end where it stops, see `Walk::new`. `taken` holds the pairs
    /// whose junctions a chain took.
    fn chain(&self, mut entry: usize, first: Option<usize>, taken: &mut [bool]) -> (Chain, usize) {
        let junction = |end: usize| {
            self.junctions[end].filter(|junction| first.is_some() || !junction.circular)
        };
        let mut visits = Vec::new();
        let mut by_junction = first.is_some();
        loop {
            visits.push(Visit {
                segment: entry / 2,
                forward: entry.is_multiple_of(2),
                by_junction,
            });
            let exit = entry ^ 1;
            let end = match junction(exit) {
                Some(junction) if !taken[junction.pair] => {
                    taken[junction.pair] = true;
                    entry = junction.other;
                    by_junction = true;
                    continue;
                }
                Some(junction) if Some(junction.pair) == first => ChainEnd::Cycle,
                _ => match self.neighbour(exit) {
                    None if first.is_some() => ChainEnd::Open,
                    None => ChainEnd::Telomere,
                    Some(next)
                        if junction(next).is_some_and(|junction| {
                            taken[junction.pair] || junction.other != next ^ 1
                        }) =>
                    {
                        ChainEnd::Open
                    }
                    Some(next) => {
                        entry = next;
                        by_junction = false;
                        continue;
                    }
                },
            };
            return (Chain { visits, end }, exit);
        }
    }

    /// The variants of the segments between two breakpoints by the chains through them: one
    /// that no chain runs through is a deletion, one that a chain runs through backward an
    /// inversion, one that a chain runs through twice in a row, through a junction from one end
    /// to the other, a tandem dup, and one that chains run through twice otherwise, or a circle
    /// once, a dup elsewhere. One that a chain runs through forward stays. A segment at a
    /// telomere gives none, as a breakend there may have a mate that is not observed, see VCF
    /// 4.5 "Single breakends", nor does one of an open chain. A variant has the members of the
    /// pairs at its two gaps, and joins its sides, see `StructuralVariant::junction`, only where
    /// they are the junctions of its kind, e.g. those of a deletion its flanks.
    fn derived(&self, members: &[EventMember], pairs: &[(usize, usize)]) -> Vec<DerivedVariant> {
        let mut visits = vec![Vec::new(); self.segments.len()];
        for (n, chain) in self.chains.iter().enumerate() {
            for (step, visit) in chain.visits.iter().enumerate() {
                visits[visit.segment].push((n, step, *visit));
            }
        }
        let mut derived = Vec::new();
        for (segment, visits) in visits.iter().enumerate() {
            let Segment {
                contig,
                front: Some(front),
                back: Some(back),
            } = self.segments[segment]
            else {
                continue;
            };
            let ends = |end: ChainEnd| visits.iter().any(|&(n, ..)| self.chains[n].end == end);
            if ends(ChainEnd::Open) {
                continue;
            }
            let kind = match visits[..] {
                [] => SvKind::Deletion,
                _ if ends(ChainEnd::Cycle) => SvKind::Duplication { tandem: false },
                [(_, _, visit)] if visit.forward => continue,
                [_] => SvKind::Inversion,
                [(n, step, visit), (m, next, again)]
                    if n == m
                        && next == step + 1
                        && visit.forward == again.forward
                        && again.by_junction =>
                {
                    SvKind::Duplication { tandem: true }
                }
                _ => SvKind::Duplication { tandem: false },
            };
            // The ends of the segment, and those of its neighbours on the contig.
            let (left, right) = (2 * segment, 2 * segment + 1);
            let (before, behind) = (left - 1, right + 1);
            let joined = |a: usize, b: usize| self.junctions[a].is_some_and(|j| j.other == b);
            let inserted = |end: usize| {
                self.junctions[end].map_or(0, |j| members[pairs[j.pair].0].breakend.inserted.len())
            };
            let (junction, inserted) = match kind {
                SvKind::Deletion if joined(before, behind) => (true, [inserted(before), 0]),
                SvKind::Duplication { tandem: true } => (true, [inserted(right), 0]),
                SvKind::Inversion if joined(before, right) && joined(left, behind) => {
                    (true, [inserted(before), inserted(left)])
                }
                SvKind::Duplication { .. } => (true, [0, 0]),
                _ => (false, [0, 0]),
            };
            let bounding = [before, left, right, behind]
                .iter()
                .filter_map(|&end| self.junctions[end].map(|j| j.pair))
                .unique()
                .sorted_by_key(|&pair| {
                    let (lo, hi) = pairs[pair];
                    (members[lo].var.position, members[hi].var.position)
                })
                .flat_map(|pair| [pairs[pair].0, pairs[pair].1])
                .collect_vec();
            let imprecise = (bounding.iter()).any(|&i| members[i].breakend.imprecise);
            let side = |gap: [usize; 2]| {
                let breakends = gap.iter().filter_map(|&end| self.junctions[end]);
                interval(breakends.map(|j| &members[j.member].breakend), imprecise)
            };
            let reference = (members.iter())
                .find(|member| member.var.chromosome == contig && member.var.position == front)
                .and_then(|member| member.var.reference.get(..1))
                .unwrap_or("N");
            derived.push(DerivedVariant {
                members: bounding,
                var: VcfVariant {
                    chromosome: contig.into(),
                    position: front,
                    reference: reference.into(),
                    alternative: alternative(kind).into(),
                },
                sv: StructuralVariant {
                    kind,
                    end: back,
                    imprecise,
                    intervals: (side([before, left]), side([right, behind])),
                    junction,
                },
                inserted,
            });
        }
        derived
    }
}

/// The events that the `pairs` of breakends among `members` on `chromosome` describe, each with
/// the members of its pairs, as symbolic alleles with their padding base.
/// With POS a < b, the pair `t[b[` at a and `]a]t` at b is the deletion of a+1 to b-1, the pair
/// `]b]t` at a and `t[a[` at b the tandem duplication of a to b, and the pairs `t]b]` at a with
/// `t]a]` at b and `[b+1[t` at a+1 with `[a+1[t` at b+1 the inversion of a+1 to b, as in VCF
/// 4.5. Bases inserted at a junction do not change the event, only its reading frame, see
/// `DerivedVariant::inserted`. Other pairs give none: a head-to-head or tail-to-tail pair alone,
/// other joins, or a deletion of no base.
fn derive_events(
    chromosome: &str,
    members: &[EventMember],
    pairs: &[(usize, usize)],
) -> Vec<DerivedVariant> {
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
    let on_contig = |(lo, hi): &(usize, usize)| {
        members[*lo].var.chromosome == chromosome && members[*hi].var.chromosome == chromosome
    };
    for &pair in pairs.iter().filter(|pair| on_contig(pair)) {
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
            // The intervals around the start and the end are those of the breakends there.
            let imprecise = breakends().any(|breakend| breakend.imprecise);
            let side = |member_of_pair: fn((usize, usize)) -> usize| {
                let breakends = pairs
                    .iter()
                    .map(|&pair| &members[member_of_pair(pair)].breakend);
                interval(breakends, imprecise)
            };
            let intervals = (side(|(lo, _)| lo), side(|(_, hi)| hi));
            let reference = padding
                .and_then(|i| members[i].var.reference.get(..1))
                .unwrap_or("N");
            // The bases inserted at the junction of each pair, as its lower breakend has them.
            let inserted = [0, 1]
                .map(|n| (pairs.get(n)).map_or(0, |&(lo, _)| members[lo].breakend.inserted.len()));
            DerivedVariant {
                members: indices,
                var: VcfVariant {
                    chromosome: chromosome.into(),
                    position,
                    reference: reference.into(),
                    alternative: alternative(kind).into(),
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

/// The symbolic allele of a derived variant of `kind`.
fn alternative(kind: SvKind) -> &'static str {
    match kind {
        SvKind::Deletion => "<DEL>",
        SvKind::Duplication { tandem: true } => "<DUP:TANDEM>",
        SvKind::Duplication { tandem: false } => "<DUP>",
        _ => "<INV>",
    }
}

/// The confidence interval of a side of a derived variant from the `breakends` there, all of
/// them, as offsets. In an `imprecise` variant, a precise breakend is exact, and an imprecise
/// one without `CIPOS` leaves the side unknown.
fn interval<'a>(
    breakends: impl Iterator<Item = &'a Breakend>,
    imprecise: bool,
) -> Option<ConfidenceInterval> {
    let mut interval = None;
    for breakend in breakends {
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
}

#[cfg(test)]
mod test {
    use super::{BreakendEvent, ChainEnd, EventMember, EventType, Walk, pairs};
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

    /// The members of `records`, as contig, position and breakend.
    fn members_on(records: &[(&str, i32, &str)]) -> Vec<EventMember> {
        (records.iter())
            .map(|&(contig, position, alt)| {
                let mut member = member(position, alt, None);
                member.var.chromosome = contig.into();
                member
            })
            .collect()
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
    /// gives none, while its group can, and so does a junction of an inversion alone, as the
    /// other end of its inverted bases is unknown.
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
            (BreakendEvent::new(members.to_vec()).derived.into_iter())
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
        let members = vec![
            member(100, "G[17:200[", Some((-5, 5))),
            member(200, "]17:100]T", None),
        ];
        let event = BreakendEvent::new(members);
        assert_eq!(
            (event.derived.iter())
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
        let events = (BreakendEvent::new(members).derived.into_iter())
            .map(|derived| (derived.members, derived.inserted))
            .collect_vec();
        assert_eq!(events, [(indices.to_vec(), inserted)]);
    }

    /// The chains through the segments of VCF 4.5 Figure 1, a 3-break operation, and Figure 7,
    /// a reciprocal translocation, each from a telomere to a telomere, with the segments as
    /// contig, first and last base, and whether the chain runs forward. All segments touch a
    /// telomere, so nothing is derived.
    #[rstest::rstest]
    #[case::figure_1(
        &[("2", 321681, "G]17:198982]"), ("2", 321682, "]13:123456]T"), ("13", 123456, "C[2:321682["), ("13", 123457, "[17:198983[A"), ("17", 198982, "A]2:321681]"), ("17", 198983, "[13:123457[C")],
        &["2:..321681+ 17:..198982-", "2:321682..- 13:..123456-", "13:123457..- 17:198983..+"],
    )]
    #[case::figure_7(
        &[("2", 321681, "G[13:123457["), ("2", 321682, "]13:123456]T"), ("13", 123456, "C[2:321682["), ("13", 123457, "]2:321681]A")],
        &["2:..321681+ 13:123457..+", "2:321682..- 13:..123456-"],
    )]
    fn chains_of_figures(#[case] records: &[(&str, i32, &str)], #[case] expected: &[&str]) {
        let members = members_on(records);
        let pairs = pairs(&members);
        let Ok(walk) = Walk::new(&members, &pairs) else {
            panic!("no walk")
        };
        let chains = (walk.chains.iter())
            .map(|chain| {
                assert_eq!(chain.end, ChainEnd::Telomere);
                (chain.visits.iter())
                    .map(|visit| {
                        let segment = &walk.segments[visit.segment];
                        let first = segment
                            .front
                            .map_or(String::new(), |gap| (gap + 1).to_string());
                        let last = segment.back.map_or(String::new(), |gap| gap.to_string());
                        let strand = if visit.forward { "+" } else { "-" };
                        format!("{}:{first}..{last}{strand}", segment.contig)
                    })
                    .join(" ")
            })
            .collect_vec();
        assert_eq!(chains, expected);
        assert!(BreakendEvent::new(members).derived.is_empty());
    }

    /// The variants of the segments between breakpoints, from the chains through them, as
    /// contig, POS, ALT, end and whether a junction joins their sides: a lost segment is a
    /// deletion, one inverted an inversion, one twice in a row through its own junction a tandem
    /// dup, and one on a circle, e.g. of a double minute, a dup elsewhere. A segment at a
    /// telomere, one moved, and one of a chain whose end is unknown give none, e.g. the
    /// fragment that a deletion takes in. A segment end with two junctions, here those of a
    /// breakend with two mates, leaves each pair to its own variant.
    #[rstest::rstest]
    #[case::unbalanced_translocation(
        &[("2", 321681, "G[13:123457["), ("13", 123457, "]2:321681]A")],
        None,
        &[],
    )]
    #[case::inversion_with_losses(
        &[("17", 100, "G]17:200]"), ("17", 110, "[17:205[T"), ("17", 200, "A]17:100]"), ("17", 205, "[17:110[C")],
        None,
        &[("17", 100, "<DEL>", 109, false), ("17", 109, "<INV>", 200, false), ("17", 200, "<DEL>", 204, false)],
    )]
    #[case::deletion_through_fragment(
        &[("17", 100, "G[13:500["), ("13", 500, "]17:100]T"), ("13", 600, "T[17:200["), ("17", 200, "]13:600]T")],
        None,
        &[("17", 100, "<DEL>", 199, false)],
    )]
    #[case::chromothripsis(
        &[("17", 100, "G[17:301["), ("17", 200, "[17:401[T"), ("17", 300, "G]17:400]"), ("17", 301, "]17:100]T"), ("17", 400, "G]17:300]"), ("17", 401, "[17:200[T")],
        None,
        &[("17", 100, "<DEL>", 199, false), ("17", 199, "<INV>", 300, false)],
    )]
    #[case::tandem_duplication(
        &[("17", 100, "]17:200]G"), ("17", 200, "T[17:100[")],
        None,
        &[("17", 99, "<DUP:TANDEM>", 200, true)],
    )]
    #[case::double_minute(
        &[("17", 100, "]17:200]G"), ("17", 200, "T[17:100[")],
        Some("DOUBLEMINUTE"),
        &[("17", 99, "<DUP>", 200, true)],
    )]
    #[case::multiple_mates(
        &[("17", 100, "G[17:200["), ("17", 100, "G[13:300["), ("17", 200, "]17:100]T"), ("13", 300, "]17:100]T")],
        None,
        &[("17", 100, "<DEL>", 199, true)],
    )]
    fn variants_of_walk(
        #[case] records: &[(&str, i32, &str)],
        #[case] event_type: Option<&str>,
        #[case] expected: &[(&str, i32, &str, i32, bool)],
    ) {
        let mut members = members_on(records);
        for member in &mut members {
            member.event_type = event_type.map(|event_type| event_type.parse().unwrap());
        }
        let derived = (BreakendEvent::new(members).derived.iter())
            .map(|derived| {
                let (var, sv) = (&derived.var, &derived.sv);
                (
                    var.chromosome.clone(),
                    var.position,
                    var.alternative.clone(),
                    sv.end,
                    sv.junction,
                )
            })
            .collect_vec();
        let expected = (expected.iter())
            .map(|&(contig, position, alt, end, junction)| {
                (contig.to_string(), position, alt.to_string(), end, junction)
            })
            .collect_vec();
        assert_eq!(derived, expected);
    }

    /// `EVENTTYPE` keeps the derived variants that fit it, as their `ALT`s: those of its kind
    /// for a simple type, as a pair that looks like a deletion can be part of a dup, a copy
    /// elsewhere for a double minute, whose junction closes a circle, and all of them for a
    /// complex rearrangement or a translocation. Other types, e.g. `INS`, and members of
    /// different types keep none.
    #[rstest::rstest]
    #[case::deletion(&[(100, "G[17:200[", Some("DEL")), (200, "]17:100]T", Some("DEL"))], &["<DEL>"])]
    #[case::deletion_as_duplication(&[(100, "G[17:200[", Some("DUP")), (200, "]17:100]T", Some("DUP"))], &[])]
    #[case::deletion_as_translocation(&[(100, "G[17:200[", Some("TRA")), (200, "]17:100]T", Some("TRA"))], &["<DEL>"])]
    #[case::different_types(&[(100, "G[17:200[", Some("DEL")), (200, "]17:100]T", Some("DUP"))], &[])]
    #[case::tandem_duplication(&[(100, "]17:200]G", Some("DUP:TANDEM")), (200, "T[17:100[", None)], &["<DUP:TANDEM>"])]
    #[case::duplication(&[(100, "]17:200]G", Some("DUP")), (200, "T[17:100[", Some("DUP"))], &["<DUP:TANDEM>"])]
    #[case::dispersed_duplication(&[(100, "]17:200]G", Some("DUP:DISPERSED")), (200, "T[17:100[", Some("DUP:DISPERSED"))], &[])]
    #[case::double_minute(&[(100, "]17:200]G", Some("DOUBLEMINUTE")), (200, "T[17:100[", Some("DOUBLEMINUTE"))], &["<DUP>"])]
    #[case::inversion(
        &[(100, "G]17:200]", Some("INV")), (101, "[17:201[T", Some("INV")), (200, "A]17:100]", Some("INV")), (201, "[17:101[C", Some("INV"))],
        &["<INV>"],
    )]
    #[case::chromothripsis(
        &[(100, "G[17:301[", Some("CHROMOTHRIPSIS")), (200, "[17:401[T", Some("CHROMOTHRIPSIS")), (300, "G]17:400]", Some("CHROMOTHRIPSIS")), (301, "]17:100]T", Some("CHROMOTHRIPSIS")), (400, "G]17:300]", Some("CHROMOTHRIPSIS")), (401, "[17:200[T", Some("CHROMOTHRIPSIS"))],
        &["<DEL>", "<INV>"],
    )]
    #[case::chromothripsis_as_insertion(
        &[(100, "G[17:301[", Some("INS")), (200, "[17:401[T", Some("INS")), (300, "G]17:400]", Some("INS")), (301, "]17:100]T", Some("INS")), (400, "G]17:300]", Some("INS")), (401, "[17:200[T", Some("INS"))],
        &[],
    )]
    #[case::per_pair(
        &[(100, "T[17:200[", Some("DEL")), (200, "]17:100]A", Some("DEL")), (300, "T[17:400[", Some("INS")), (400, "]17:300]A", Some("INS"))],
        &["<DEL>"],
    )]
    fn event_type_restricts_event(
        #[case] records: &[(i32, &str, Option<&str>)],
        #[case] alternatives: &[&str],
    ) {
        let members = (records.iter())
            .map(|&(position, alt, event_type)| EventMember {
                event_type: event_type.map(|event_type| event_type.parse().unwrap()),
                ..member(position, alt, None)
            })
            .collect_vec();
        let event = BreakendEvent::new(members);
        let derived = (event.derived.iter()).map(|derived| derived.var.alternative.as_str());
        assert_eq!(derived.collect_vec(), alternatives);
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
