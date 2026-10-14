//! Transcript selection for `--restrict-transcripts`, `--pick`, `--pick-order` and
//! `--flag-pick`.

use std::cmp::{Ordering, Reverse};

use hgvs::data::interface::TxForRegionRecord;
use indexmap::IndexMap;

use super::plugin::TxAnnotation;
use super::terms::{ANN_PICK, Consequence, Message};
use crate::annotate::cli::{PickCriterion, PickUnit, TagCriterion, TranscriptSettings};
use crate::annotate::seqvars::provider::PbsTranscriptExt;
use crate::errors::SeqvarsError;
use crate::pbs::txs::{Transcript, TranscriptBiotype, TranscriptTag};

/// A transcript near the variant, annotated on demand.
pub(crate) struct Candidate<'a> {
    record: TxForRegionRecord,
    tx: &'a Transcript,
    /// `None` until annotated, then the annotation if the transcript gives one.
    ann: Option<Option<TxAnnotation>>,
}

impl<'a> Candidate<'a> {
    pub(crate) fn new(record: TxForRegionRecord, tx: &'a Transcript) -> Self {
        Self {
            record,
            tx,
            ann: None,
        }
    }

    fn annotate<F>(&mut self, annotate: &mut F) -> Result<(), SeqvarsError>
    where
        F: FnMut(&TxForRegionRecord) -> Result<Option<TxAnnotation>, SeqvarsError>,
    {
        if self.ann.is_none() {
            self.ann = Some(annotate(&self.record)?);
        }
        Ok(())
    }

    /// Whether the candidate was annotated but gives no annotation.
    fn is_void(&self) -> bool {
        matches!(self.ann, Some(None))
    }

    /// Most severe consequence, larger is more severe; `None` if not annotated.
    fn severity(&self) -> Option<Reverse<Consequence>> {
        match &self.ann {
            Some(Some(tx)) => Some(Reverse(
                tx.ann
                    .consequences
                    .first()
                    .copied()
                    .unwrap_or(Consequence::GeneVariant),
            )),
            _ => None,
        }
    }
}

/// The annotations to report, in candidate order.
pub(crate) struct Selection {
    pub anns: Vec<TxAnnotation>,
    /// Number of genes dropped by `--restrict-transcripts`.
    pub dropped_genes: usize,
}

/// Annotate `cands` and select the annotations to report according to `settings`.
///
/// `annotate` is only called where an annotation is needed: for the reported candidates and
/// for those still tied when `--pick-order` reaches `severity`.
pub(crate) fn select<F>(
    settings: &TranscriptSettings,
    mut cands: Vec<Candidate<'_>>,
    mut annotate: F,
) -> Result<Selection, SeqvarsError>
where
    F: FnMut(&TxForRegionRecord) -> Result<Option<TxAnnotation>, SeqvarsError>,
{
    let order = &settings.pick_order;
    let mut report = vec![true; cands.len()];
    let mut fallback = vec![false; cands.len()];
    let mut flagged = vec![false; cands.len()];
    let mut dropped_genes = 0;

    if !settings.restrict_transcripts.is_empty() {
        let all: Vec<_> = (0..cands.len()).collect();
        for (gene, idxs) in group_by(&cands, &all, |cand| gene_key(cand.tx).to_string()) {
            let allowed: Vec<_> = idxs
                .iter()
                .map(|&i| {
                    settings
                        .restrict_transcripts
                        .iter()
                        .any(|&tag| has_tag(cands[i].tx, tag))
                })
                .collect();
            if allowed.contains(&true) {
                for (&i, &allowed) in idxs.iter().zip(&allowed) {
                    report[i] = allowed;
                }
            } else if settings.restrict_fallback {
                let top = best(&mut cands, &idxs, order, &mut annotate)?;
                for &i in &idxs {
                    report[i] = Some(i) == top;
                    fallback[i] = Some(i) == top;
                }
            } else {
                tracing::debug!("no transcript of gene {gene} has the --restrict-transcripts tags");
                dropped_genes += 1;
                for &i in &idxs {
                    report[i] = false;
                }
            }
        }
    }

    if let Some(unit) = settings.pick {
        let pool: Vec<_> = (0..cands.len()).filter(|&i| report[i]).collect();
        let mut picked = vec![false; cands.len()];
        for (_, idxs) in group_by(&cands, &pool, |cand| unit_key(unit, cand)) {
            if let Some(top) = best(&mut cands, &idxs, order, &mut annotate)? {
                picked[top] = true;
            }
        }
        if settings.flag_pick {
            flagged = picked;
        } else {
            report = picked;
        }
    }

    let mut anns = Vec::new();
    for (i, cand) in cands.iter_mut().enumerate() {
        if !report[i] {
            continue;
        }
        cand.annotate(&mut annotate)?;
        let Some(Some(mut tx)) = cand.ann.take() else {
            continue;
        };
        if fallback[i] {
            tx.ann
                .messages
                .get_or_insert_with(Vec::new)
                .push(Message::InfoRestrictFallback);
        }
        if flagged[i] {
            tx.ann
                .custom_fields
                .insert(ANN_PICK.to_string(), Some("1".to_string()));
        }
        anns.push(tx);
    }

    Ok(Selection {
        anns,
        dropped_genes,
    })
}

/// Index of the top-ranked candidate among `idxs` that gives an annotation.
fn best<F>(
    cands: &mut [Candidate<'_>],
    idxs: &[usize],
    order: &[PickCriterion],
    annotate: &mut F,
) -> Result<Option<usize>, SeqvarsError>
where
    F: FnMut(&TxForRegionRecord) -> Result<Option<TxAnnotation>, SeqvarsError>,
{
    let mut pool = idxs.to_vec();
    while !pool.is_empty() {
        let mut tied = pool.clone();
        for &criterion in order {
            if tied.len() <= 1 {
                break;
            }
            if criterion == PickCriterion::Severity {
                for &i in &tied {
                    cands[i].annotate(annotate)?;
                }
                tied.retain(|&i| !cands[i].is_void());
            }
            if let Some(top) = tied
                .iter()
                .copied()
                .max_by(|&a, &b| compare(criterion, &cands[a], &cands[b]))
            {
                tied.retain(|&i| compare(criterion, &cands[i], &cands[top]).is_eq());
            }
        }
        if let Some(top) = tied
            .into_iter()
            .min_by(|&a, &b| cands[a].record.tx_ac.cmp(&cands[b].record.tx_ac))
        {
            cands[top].annotate(annotate)?;
            if !cands[top].is_void() {
                return Ok(Some(top));
            }
        }
        // Rank again without the candidates that turned out to give no annotation.
        pool.retain(|&i| !cands[i].is_void());
    }
    Ok(None)
}

/// Compare `a` with `b` by `criterion`; `Greater` if `a` ranks higher.
fn compare(criterion: PickCriterion, a: &Candidate<'_>, b: &Candidate<'_>) -> Ordering {
    match criterion {
        PickCriterion::Coding => is_coding(a.tx).cmp(&is_coding(b.tx)),
        PickCriterion::Unflagged => a.tx.is_clean().cmp(&b.tx.is_clean()),
        PickCriterion::Severity => a.severity().cmp(&b.severity()),
        PickCriterion::Length => a.tx.tx_len().cmp(&b.tx.tx_len()),
        _ => match as_tag(criterion) {
            Some(tag) => has_tag(a.tx, tag).cmp(&has_tag(b.tx, tag)),
            None => Ordering::Equal,
        },
    }
}

fn as_tag(criterion: PickCriterion) -> Option<TagCriterion> {
    match criterion {
        PickCriterion::ManeSelect => Some(TagCriterion::ManeSelect),
        PickCriterion::ManeSelectBackport => Some(TagCriterion::ManeSelectBackport),
        PickCriterion::ManePlusClinical => Some(TagCriterion::ManePlusClinical),
        PickCriterion::ManePlusClinicalBackport => Some(TagCriterion::ManePlusClinicalBackport),
        PickCriterion::RefseqSelect => Some(TagCriterion::RefseqSelect),
        PickCriterion::RefseqSelectBackport => Some(TagCriterion::RefseqSelectBackport),
        PickCriterion::EnsemblCanonical => Some(TagCriterion::EnsemblCanonical),
        PickCriterion::EnsemblCanonicalBackport => Some(TagCriterion::EnsemblCanonicalBackport),
        PickCriterion::GencodePrimary => Some(TagCriterion::GencodePrimary),
        PickCriterion::GencodePrimaryBackport => Some(TagCriterion::GencodePrimaryBackport),
        PickCriterion::Basic => Some(TagCriterion::Basic),
        PickCriterion::BasicBackport => Some(TagCriterion::BasicBackport),
        PickCriterion::Coding
        | PickCriterion::Unflagged
        | PickCriterion::Severity
        | PickCriterion::Length => None,
    }
}

/// Whether `tx` has `tag`; a tag also matches its backport.
fn has_tag(tx: &Transcript, tag: TagCriterion) -> bool {
    let (native, backport) = match tag {
        TagCriterion::ManeSelect => (
            Some(TranscriptTag::ManeSelect),
            TranscriptTag::ManeSelectBackport,
        ),
        TagCriterion::ManeSelectBackport => (None, TranscriptTag::ManeSelectBackport),
        TagCriterion::ManePlusClinical => (
            Some(TranscriptTag::ManePlusClinical),
            TranscriptTag::ManePlusClinicalBackport,
        ),
        TagCriterion::ManePlusClinicalBackport => (None, TranscriptTag::ManePlusClinicalBackport),
        TagCriterion::RefseqSelect => (
            Some(TranscriptTag::RefSeqSelect),
            TranscriptTag::RefSeqSelectBackport,
        ),
        TagCriterion::RefseqSelectBackport => (None, TranscriptTag::RefSeqSelectBackport),
        TagCriterion::EnsemblCanonical => (
            Some(TranscriptTag::EnsemblCanonical),
            TranscriptTag::EnsemblCanonicalBackport,
        ),
        TagCriterion::EnsemblCanonicalBackport => (None, TranscriptTag::EnsemblCanonicalBackport),
        TagCriterion::GencodePrimary => (
            Some(TranscriptTag::GencodePrimary),
            TranscriptTag::GencodePrimaryBackport,
        ),
        TagCriterion::GencodePrimaryBackport => (None, TranscriptTag::GencodePrimaryBackport),
        TagCriterion::Basic => (Some(TranscriptTag::Basic), TranscriptTag::BasicBackport),
        TagCriterion::BasicBackport => (None, TranscriptTag::BasicBackport),
    };
    tx.tags
        .iter()
        .filter_map(|&tag| TranscriptTag::try_from(tag).ok())
        .any(|tag| Some(tag) == native || tag == backport)
}

fn is_coding(tx: &Transcript) -> bool {
    tx.biotype == TranscriptBiotype::Coding as i32
}

/// The key that groups the transcripts of one gene.
fn gene_key(tx: &Transcript) -> &str {
    &tx.gene_id
}

fn unit_key(unit: PickUnit, cand: &Candidate<'_>) -> String {
    match unit {
        PickUnit::Gene => gene_key(cand.tx).to_string(),
        PickUnit::Transcript => cand.record.tx_ac.clone(),
        PickUnit::Allele => String::new(),
    }
}

/// Group `idxs` by `key`, in order of first appearance.
fn group_by(
    cands: &[Candidate<'_>],
    idxs: &[usize],
    key: impl Fn(&Candidate<'_>) -> String,
) -> IndexMap<String, Vec<usize>> {
    let mut groups: IndexMap<String, Vec<usize>> = IndexMap::new();
    for &i in idxs {
        groups.entry(key(&cands[i])).or_default().push(i);
    }
    groups
}

#[cfg(test)]
mod test {
    use std::cell::RefCell;

    use super::*;
    use crate::annotate::cli::DEFAULT_PICK_ORDER;
    use crate::annotate::seqvars::consequence::plugin::TxContext;
    use crate::annotate::seqvars::consequence::terms::AnnField;
    use crate::pbs::txs::{ExonAlignment, GenomeAlignment};

    use PickCriterion as P;
    use TranscriptTag as T;

    fn tx(id: &str, gene: &str, tags: &[TranscriptTag], len: i32) -> Transcript {
        Transcript {
            id: id.into(),
            gene_id: gene.into(),
            biotype: TranscriptBiotype::Coding as i32,
            tags: tags.iter().map(|&tag| tag as i32).collect(),
            genome_alignments: vec![GenomeAlignment {
                exons: vec![ExonAlignment {
                    alt_cds_start_i: Some(0),
                    alt_cds_end_i: Some(len),
                    ..Default::default()
                }],
                ..Default::default()
            }],
            ..Default::default()
        }
    }

    fn settings(pick: Option<PickUnit>, order: &[PickCriterion]) -> TranscriptSettings {
        TranscriptSettings {
            pick,
            pick_order: order.to_vec(),
            ..Default::default()
        }
    }

    fn annotation(ann: AnnField) -> TxAnnotation {
        TxAnnotation {
            context: TxContext::empty(),
            ann,
        }
    }

    /// Outcome of [`run`]: reported feature IDs and the IDs that were annotated.
    struct Outcome {
        selection: Selection,
        annotated: Vec<String>,
    }

    impl Outcome {
        fn ids(&self) -> Vec<&str> {
            self.selection
                .anns
                .iter()
                .map(|tx| tx.ann.feature_id.as_str())
                .collect()
        }
    }

    /// Select among `txs`; each annotation has the consequence from `csqs`, or none if missing.
    fn run(
        settings: &TranscriptSettings,
        txs: &[Transcript],
        csqs: &[(&str, Consequence)],
    ) -> Outcome {
        let annotated = RefCell::new(Vec::new());
        let cands = txs
            .iter()
            .map(|tx| {
                let record = TxForRegionRecord {
                    tx_ac: tx.id.clone(),
                    alt_ac: "NC_000001.11".into(),
                    alt_strand: 1,
                    alt_aln_method: "splign".into(),
                    start_i: 0,
                    end_i: 100,
                };
                Candidate::new(record, tx)
            })
            .collect();
        let selection = select(settings, cands, |record| {
            annotated.borrow_mut().push(record.tx_ac.clone());
            Ok(csqs
                .iter()
                .find(|(id, _)| *id == record.tx_ac)
                .map(|&(id, csq)| {
                    annotation(AnnField {
                        feature_id: id.into(),
                        consequences: vec![csq],
                        ..Default::default()
                    })
                }))
        })
        .unwrap();
        Outcome {
            selection,
            annotated: annotated.into_inner(),
        }
    }

    const MISSENSE: Consequence = Consequence::MissenseVariant;
    const STOP: Consequence = Consequence::StopGained;

    #[test]
    fn no_options_report_all() {
        let txs = [tx("NM_2.1", "G", &[], 10), tx("NM_1.1", "G", &[], 10)];
        let outcome = run(
            &TranscriptSettings::default(),
            &txs,
            &[("NM_1.1", MISSENSE), ("NM_2.1", MISSENSE)],
        );
        assert_eq!(outcome.ids(), ["NM_2.1", "NM_1.1"]);
    }

    #[test]
    fn later_criteria_break_ties() {
        let txs = [
            tx("NM_1.1", "G", &[T::Basic], 100),
            tx("NM_2.1", "G", &[T::Basic, T::EnsemblCanonical], 50),
        ];
        let csqs = [("NM_1.1", MISSENSE), ("NM_2.1", MISSENSE)];
        let order = [P::Basic, P::EnsemblCanonical, P::Length];
        let outcome = run(&settings(Some(PickUnit::Gene), &order), &txs, &csqs);
        assert_eq!(outcome.ids(), ["NM_2.1"]);
    }

    #[test]
    fn lowest_accession_breaks_remaining_ties() {
        let txs = [tx("NM_2.1", "G", &[], 10), tx("NM_1.1", "G", &[], 10)];
        let csqs = [("NM_1.1", MISSENSE), ("NM_2.1", MISSENSE)];
        let outcome = run(&settings(Some(PickUnit::Gene), &[P::Length]), &txs, &csqs);
        assert_eq!(outcome.ids(), ["NM_1.1"]);
    }

    #[test]
    fn order_decides_between_severity_and_mane() {
        let txs = [
            tx("NM_1.1", "G", &[T::ManeSelect], 10),
            tx("NM_2.1", "G", &[], 10),
        ];
        let csqs = [("NM_1.1", MISSENSE), ("NM_2.1", STOP)];
        let mane_first = [P::ManeSelect, P::Severity];
        let severity_first = [P::Severity, P::ManeSelect];
        let gene = Some(PickUnit::Gene);
        assert_eq!(
            run(&settings(gene, &mane_first), &txs, &csqs).ids(),
            ["NM_1.1"]
        );
        assert_eq!(
            run(&settings(gene, &severity_first), &txs, &csqs).ids(),
            ["NM_2.1"]
        );
    }

    #[test]
    fn tag_only_order_annotates_only_the_pick() {
        let txs = [
            tx("NM_1.1", "G", &[], 10),
            tx("NM_2.1", "G", &[T::ManeSelect], 10),
            tx("NM_3.1", "G", &[], 10),
        ];
        let csqs = [("NM_1.1", STOP), ("NM_2.1", MISSENSE), ("NM_3.1", STOP)];
        let outcome = run(
            &settings(Some(PickUnit::Gene), &[P::ManeSelect]),
            &txs,
            &csqs,
        );
        assert_eq!(outcome.ids(), ["NM_2.1"]);
        assert_eq!(outcome.annotated, ["NM_2.1"]);
    }

    #[test]
    fn severity_annotates_only_the_tied_candidates() {
        let txs = [
            tx("NM_1.1", "G", &[T::ManeSelect], 10),
            tx("NM_2.1", "G", &[T::ManeSelect], 10),
            tx("NM_3.1", "G", &[], 10),
        ];
        let csqs = [("NM_1.1", MISSENSE), ("NM_2.1", STOP), ("NM_3.1", STOP)];
        let order = [P::ManeSelect, P::Severity];
        let outcome = run(&settings(Some(PickUnit::Gene), &order), &txs, &csqs);
        assert_eq!(outcome.ids(), ["NM_2.1"]);
        assert_eq!(outcome.annotated, ["NM_1.1", "NM_2.1"]);
    }

    #[test]
    fn one_pick_per_gene_and_allele() {
        let txs = [
            tx("NM_1.1", "G1", &[], 20),
            tx("NM_2.1", "G1", &[], 10),
            tx("NM_3.1", "G2", &[], 10),
        ];
        let csqs = [("NM_1.1", MISSENSE), ("NM_2.1", MISSENSE), ("NM_3.1", STOP)];
        let by_gene = run(
            &settings(Some(PickUnit::Gene), DEFAULT_PICK_ORDER),
            &txs,
            &csqs,
        );
        assert_eq!(by_gene.ids(), ["NM_1.1", "NM_3.1"]);
        let by_allele = run(
            &settings(Some(PickUnit::Allele), DEFAULT_PICK_ORDER),
            &txs,
            &csqs,
        );
        assert_eq!(by_allele.ids(), ["NM_3.1"]);
    }

    #[test]
    fn candidate_without_annotation_falls_through() {
        let txs = [
            tx("NM_1.1", "G", &[T::ManeSelect], 10),
            tx("NM_2.1", "G", &[], 10),
        ];
        let outcome = run(
            &settings(Some(PickUnit::Gene), &[P::ManeSelect]),
            &txs,
            &[("NM_2.1", MISSENSE)],
        );
        assert_eq!(outcome.ids(), ["NM_2.1"]);
    }

    #[test]
    fn transcript_unit_keeps_one_annotation_per_transcript() {
        let tx = tx("NM_1.1", "G", &[], 10);
        let record = |start_i| TxForRegionRecord {
            tx_ac: "NM_1.1".into(),
            alt_ac: "NC_000001.11".into(),
            alt_strand: 1,
            alt_aln_method: "splign".into(),
            start_i,
            end_i: 100,
        };
        let cands = vec![
            Candidate::new(record(0), &tx),
            Candidate::new(record(1), &tx),
        ];
        let selection = select(
            &settings(Some(PickUnit::Transcript), &[P::Severity]),
            cands,
            |record| {
                Ok(Some(annotation(AnnField {
                    feature_id: "NM_1.1".into(),
                    consequences: vec![if record.start_i == 0 { MISSENSE } else { STOP }],
                    ..Default::default()
                })))
            },
        )
        .unwrap();
        assert_eq!(selection.anns.len(), 1);
        assert_eq!(selection.anns[0].ann.consequences, [STOP]);
    }

    #[test]
    fn restrict_drops_genes_without_the_tag() {
        let txs = [
            tx("NM_1.1", "G1", &[T::ManeSelect], 10),
            tx("NM_2.1", "G1", &[], 10),
            tx("NM_3.1", "G2", &[], 10),
        ];
        let csqs = [("NM_1.1", MISSENSE), ("NM_2.1", MISSENSE), ("NM_3.1", STOP)];
        let settings = TranscriptSettings {
            restrict_transcripts: vec![TagCriterion::ManeSelect],
            ..Default::default()
        };
        let outcome = run(&settings, &txs, &csqs);
        assert_eq!(outcome.ids(), ["NM_1.1"]);
        assert_eq!(outcome.selection.dropped_genes, 1);
    }

    #[test]
    fn restrict_fallback_keeps_and_marks_the_top_ranked() {
        let txs = [
            tx("NM_1.1", "G1", &[T::ManeSelect], 10),
            tx("NM_2.1", "G2", &[], 10),
            tx("NM_3.1", "G2", &[], 20),
        ];
        let csqs = [
            ("NM_1.1", MISSENSE),
            ("NM_2.1", MISSENSE),
            ("NM_3.1", MISSENSE),
        ];
        let settings = TranscriptSettings {
            restrict_transcripts: vec![TagCriterion::ManeSelect],
            restrict_fallback: true,
            ..Default::default()
        };
        let outcome = run(&settings, &txs, &csqs);
        assert_eq!(outcome.ids(), ["NM_1.1", "NM_3.1"]);
        assert_eq!(outcome.selection.dropped_genes, 0);
        assert_eq!(outcome.selection.anns[0].ann.messages, None);
        assert_eq!(
            outcome.selection.anns[1].ann.messages,
            Some(vec![Message::InfoRestrictFallback])
        );
    }

    #[test]
    fn flag_pick_keeps_all_and_marks_one_per_unit() {
        let txs = [
            tx("NM_1.1", "G1", &[], 10),
            tx("NM_2.1", "G1", &[T::ManeSelect], 10),
            tx("NM_3.1", "G2", &[], 10),
        ];
        let csqs = [
            ("NM_1.1", MISSENSE),
            ("NM_2.1", MISSENSE),
            ("NM_3.1", MISSENSE),
        ];
        let settings = TranscriptSettings {
            flag_pick: true,
            ..settings(Some(PickUnit::Gene), DEFAULT_PICK_ORDER)
        };
        let outcome = run(&settings, &txs, &csqs);
        let picks: Vec<_> = outcome
            .selection
            .anns
            .iter()
            .map(|tx| tx.ann.custom_fields.get(ANN_PICK).cloned().flatten())
            .collect();
        assert_eq!(outcome.ids(), ["NM_1.1", "NM_2.1", "NM_3.1"]);
        assert_eq!(picks, [None, Some("1".to_string()), Some("1".to_string())]);
    }

    #[test]
    fn tag_matches_its_backport_but_not_vice_versa() {
        let native = tx("NM_1.1", "G", &[T::ManeSelect], 10);
        let backport = tx("NM_2.1", "G", &[T::ManeSelectBackport], 10);
        assert!(has_tag(&backport, TagCriterion::ManeSelect));
        assert!(has_tag(&backport, TagCriterion::ManeSelectBackport));
        assert!(!has_tag(&native, TagCriterion::ManeSelectBackport));
    }
}
