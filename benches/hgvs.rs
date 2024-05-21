use biocommons_bioutils::assemblies::Assembly;
use criterion::{criterion_group, criterion_main, Criterion};
use hgvs::mapper::alignment;
use hgvs::mapper::{assembly, variant};
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

    c.bench_function("instantiate-assembly-mapper", |b| {
        b.iter(|| assembly::Mapper::new(Default::default(), provider.clone()))
    });
}

fn alignment_mapper(c: &mut Criterion) {
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
    let tx_ac = "NM_178434.2";
    let alt_ac = "NC_000001.10";
    c.bench_function("instantiate-alignment-mapper", |b| {
        b.iter(|| {
            alignment::Mapper::new(
                &Default::default(),
                provider.clone(),
                tx_ac,
                alt_ac,
                "not-transcript",
            )
        })
    });
    c.bench_function("instantiate-alignment-mapper-transcript", |b| {
        b.iter(|| {
            alignment::Mapper::new(
                &Default::default(),
                provider.clone(),
                tx_ac,
                alt_ac,
                "transcript",
            )
        })
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default().with_profiler(perf::FlamegraphProfiler::new(100));
    targets = assembly_mapper, alignment_mapper
);
criterion_main!(benches);

pub mod perf {
    use std::{fs::File, os::raw::c_int, path::Path};

    use criterion::profiler::Profiler;
    use pprof::ProfilerGuard;

    /// Small custom profiler that can be used with Criterion to create a flamegraph for benchmarks.
    /// Also see [the Criterion documentation on this][custom-profiler].
    ///
    /// ## Example on how to enable the custom profiler:
    ///
    /// ```
    /// mod perf;
    /// use perf::FlamegraphProfiler;
    ///
    /// fn fibonacci_profiled(criterion: &mut Criterion) {
    ///     // Use the criterion struct as normal here.
    /// }
    ///
    /// fn custom() -> Criterion {
    ///     Criterion::default().with_profiler(FlamegraphProfiler::new())
    /// }
    ///
    /// criterion_group! {
    ///     name = benches;
    ///     config = custom();
    ///     targets = fibonacci_profiled
    /// }
    /// ```
    ///
    /// The neat thing about this is that it will sample _only_ the benchmark, and not other stuff like
    /// the setup process.
    ///
    /// Further, it will only kick in if `--profile-time <time>` is passed to the benchmark binary.
    /// A flamegraph will be created for each individual benchmark in its report directory under
    /// `profile/flamegraph.svg`.
    ///
    /// [custom-profiler]: https://bheisler.github.io/criterion.rs/book/user_guide/profiling.html#implementing-in-process-profiling-hooks
    pub struct FlamegraphProfiler<'a> {
        frequency: c_int,
        active_profiler: Option<ProfilerGuard<'a>>,
    }

    impl<'a> FlamegraphProfiler<'a> {
        #[allow(dead_code)]
        pub fn new(frequency: c_int) -> Self {
            FlamegraphProfiler {
                frequency,
                active_profiler: None,
            }
        }
    }

    impl<'a> Profiler for FlamegraphProfiler<'a> {
        fn start_profiling(&mut self, _benchmark_id: &str, _benchmark_dir: &Path) {
            self.active_profiler = Some(ProfilerGuard::new(self.frequency).unwrap());
        }

        fn stop_profiling(&mut self, _benchmark_id: &str, benchmark_dir: &Path) {
            std::fs::create_dir_all(benchmark_dir).unwrap();
            let flamegraph_path = benchmark_dir.join("flamegraph.svg");
            let flamegraph_file = File::create(&flamegraph_path)
                .expect("File system error while creating flamegraph.svg");
            if let Some(profiler) = self.active_profiler.take() {
                profiler
                    .report()
                    .build()
                    .unwrap()
                    .flamegraph(flamegraph_file)
                    .expect("Error writing flamegraph");
            }
        }
    }
}
