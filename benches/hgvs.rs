use biocommons_bioutils::assemblies::Assembly;
use criterion::{criterion_group, criterion_main, Criterion};
use hgvs::mapper::assembly::Mapper;
use mehari::annotate::seqvars::load_tx_db;
use mehari::annotate::seqvars::provider::ConfigBuilder as MehariProviderConfigBuilder;
use mehari::annotate::seqvars::provider::Provider as MehariProvider;
use std::sync::Arc;

fn assembly_mapper(c: &mut Criterion) {
    let tx_path = "/mnt/data/mehari/0.21.0/db/grch37/txs.bin.zst";
    let tx_db = load_tx_db(tx_path).unwrap();
    let provider = Arc::new(MehariProvider::new(
        tx_db,
        Assembly::Grch37p10,
        MehariProviderConfigBuilder::default()
            .transcript_picking(false)
            .build()
            .unwrap(),
    ));

    c.bench_function("assembly::Mapper::new", |b| {
        b.iter(|| Mapper::new(Default::default(), provider.clone()))
    });
}

criterion_group!(benches, assembly_mapper);
criterion_main!(benches);
