use criterion::{black_box, criterion_group, criterion_main, Criterion};

use anyhow::Result;
use extism::convert::Json;
use extism::{FromBytesOwned, Manifest, Plugin, ToBytes, Wasm};
use mehari_plugins::*;
use std::path::PathBuf;

fn call_nmd_plugin(tva: &TranscriptVariationAllele, plugin: &mut Plugin) -> Result<bool> {
    if let Some(Json(result)) = plugin
        .call::<Json<TranscriptVariationAllele>, Option<Json<Annotation>>>(
            "process",
            Json(tva.clone()),
        )?
    {
        Ok(result.annotation == "NMD")
    } else {
        Ok(false)
    }
}

fn call_nmd_native_with_conversion(tva: &TranscriptVariationAllele) -> Result<bool> {
    let tva_json = Json(tva).to_bytes().unwrap();
    let result = call_nmd_native(tva);
    let Json(_tva) =
        black_box(Json::<TranscriptVariationAllele>::from_bytes_owned(&tva_json).unwrap());
    result
}

fn call_nmd_native(tva: &TranscriptVariationAllele) -> Result<bool> {
    fn transcript_escapes_nmd(tva: &TranscriptVariationAllele) -> bool {
        // To qualify for NMD, at least one of the following consequences is required:
        //   "stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant",
        if !tva
            .overlap_consequences()
            .iter()
            .any(|oc| INCLUDE_SO.contains(&oc.so_term().as_ref()))
        {
            return false;
        }

        let transcript = tva.transcript();
        let exons = transcript.exons();
        let strand = transcript.strand();
        let variant_feature_end_position = tva.variation_feature().seq_region().end();

        // Rules for NMD prediction
        variant_within_last_exon(variant_feature_end_position, exons)
            || variant_upstream_of_penultimate_exon(variant_feature_end_position, strand, exons)
            || variant_within_first_100_coding_bases(tva.transcript_variation())
            || variant_is_intronless(transcript)
    }

    // Included SO terms
    const INCLUDE_SO: [&str; 4] = [
        "stop_gained",
        "frameshift_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
    ];

    /// Checks whether the variant location falls within the last exon of the transcript
    fn variant_within_last_exon(variant_feature_end_position: u64, exons: &[Exon]) -> bool {
        let vf_end = variant_feature_end_position;
        exons
            .last()
            .map(|exon| vf_end >= exon.start() && vf_end <= exon.end())
            .unwrap_or(false)
    }

    fn variant_upstream_of_penultimate_exon(
        variant_feature_end_position: u64,
        strand: isize,
        exons: &[Exon],
    ) -> bool {
        let vf_end = variant_feature_end_position;
        exons
            .get(exons.len() - 2)
            .map(|second_last_exon| {
                let coding_region_end = if strand == -1 {
                    second_last_exon.start().saturating_add(51)
                } else {
                    second_last_exon.end().saturating_sub(51)
                };
                vf_end >= coding_region_end && vf_end <= second_last_exon.end()
            })
            .unwrap_or(false)
    }

    fn variant_within_first_100_coding_bases(tv: &TranscriptVariation) -> bool {
        tv.cds().end() <= 101
    }

    fn variant_is_intronless(transcript: &Transcript) -> bool {
        !transcript.has_introns()
    }
    Ok(transcript_escapes_nmd(tva))
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("wasm-vs-native");
    group.significance_level(0.1).sample_size(1000);

    // setup input data
    let tva = TranscriptVariationAlleleBuilder::default()
        .transcript(
            TranscriptBuilder::default()
                .introns(vec![Intron::new(10, 20), Intron::new(30, 40)])
                .exons(vec![
                    Exon::new(0, 10),
                    Exon::new(20, 30),
                    Exon::new(40, 100),
                ])
                .strand(-1)
                .build()
                .unwrap(),
        )
        .transcript_variation(TranscriptVariation::new(CodingSequence::new(0, 100)))
        .variation_feature(VariationFeature::new(SeqRegion::new(15, 25)))
        .overlap_consequences(vec![OverlapConsequence::new("stop_gained".into())])
        .build()
        .unwrap();

    // initialize plugin only once
    let mut plugin_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    plugin_path.push("../target/wasm32-unknown-unknown/release/mehari_plugin_vep_nmd.wasm");
    let url = Wasm::file(plugin_path);
    let manifest = Manifest::new([url]);
    let mut plugin = Plugin::new(manifest, [], true).unwrap();

    group.bench_function("call nmd plugin wasm optimized", |b| {
        b.iter(|| call_nmd_plugin(black_box(&tva), black_box(&mut plugin)))
    });

    group.bench_function("call nmd plugin native", |b| {
        // cloning here because `tva` is cloned inside `call_nmd_plugin`
        // (and we do not want to measure the difference clone makes)
        b.iter(|| call_nmd_native(black_box(&(tva.clone()))))
    });

    group.bench_function("call nmd plugin with conversion", |b| {
        // cloning here because `tva` is cloned inside `call_nmd_plugin`
        // (and we do not want to measure the difference clone makes)
        b.iter(|| call_nmd_native_with_conversion(black_box(&(tva.clone()))))
    });

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
