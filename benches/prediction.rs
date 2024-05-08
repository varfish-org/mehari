use std::sync::Arc;

use biocommons_bioutils::assemblies::Assembly;
use criterion::{criterion_group, criterion_main, BatchSize, Criterion};

use mehari::annotate::seqvars::csq::{ConfigBuilder, ConsequencePredictor, VcfVariant};
use mehari::annotate::seqvars::load_tx_db;
use mehari::annotate::seqvars::provider::ConfigBuilder as MehariProviderConfigBuilder;
use mehari::annotate::seqvars::provider::Provider as MehariProvider;

fn csq_prediction(c: &mut Criterion) {
    let mut group = c.benchmark_group("csq");
    group.sample_size(500);
    group.confidence_level(0.99);

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

    let vars = vec![
        "17:41197701:G:C",
        "17:41196309:G:C",
        "17:41196310:G:C",
        "17:41196311:G:C",
        "17:41196312:G:C",
        "17:41196313:G:C",
        "17:41197818:G:C",
        "17:41197819:G:C",
        "17:41197820:G:C",
        "17:41197821:G:C",
        "17:41197822:G:C",
        "17:41197823:G:C",
        "17:41277379:A:C",
        "17:41277380:G:C",
        "17:41277381:G:T",
        "17:41277382:G:C",
        "17:41277383:A:C",
        "17:41277384:G:C",
        "2:179631246:G:A",
    ];
    let vars: Vec<VcfVariant> = vars
        .iter()
        .map(|spdi| {
            let spdi = spdi.split(':').map(|s| s.to_string()).collect::<Vec<_>>();
            VcfVariant {
                chromosome: spdi[0].clone(),
                position: spdi[1].parse().unwrap(),
                reference: spdi[2].clone(),
                alternative: spdi[3].clone(),
            }
        })
        .collect();

    let predictor = ConsequencePredictor::new(
        provider,
        Assembly::Grch37p10,
        ConfigBuilder::default()
            .report_all_transcripts(true)
            .build()
            .unwrap(),
    );

    group.bench_function("predict", |b| {
        b.iter_batched(
            || vars.clone(),
            |vars| {
                vars.iter().for_each(|var| {
                    predictor.predict(var).unwrap();
                })
            },
            BatchSize::SmallInput,
        )
    });
    group.finish()
}

criterion_group!(
    name = benches;
    config = Criterion::default().with_profiler(perf::FlamegraphProfiler::new(100));
    targets = csq_prediction
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
            let flamegraph_file = File::create(flamegraph_path)
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
