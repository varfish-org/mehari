use std::str::FromStr;
use std::sync::Arc;

use annonars::common::cli::is_canonical;
use annonars::common::keys;
use anyhow::Result;
use biocommons_bioutils::assemblies::Assembly;
use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use noodles::vcf::record::info::field;
use noodles::vcf::Record;
use rocksdb::{BoundColumnFamily, DBWithThreadMode, MultiThreaded};

use mehari::annotate::seqvars::csq::{
    ConfigBuilder as ConsequencePredictorConfigBuilder, ConsequencePredictor, TranscriptSource,
    VcfVariant,
};
use mehari::annotate::seqvars::provider::{
    ConfigBuilder as MehariProviderConfigBuilder, Provider as MehariProvider,
};
use mehari::annotate::seqvars::{
    annotate_record_auto, annotate_record_clinvar, annotate_record_mt, annotate_record_xy,
    load_tx_db, CHROM_AUTO, CHROM_MT, CHROM_XY,
};
use mehari::pbs::txs::TxSeqDatabase;

const ROCKSDB_SEQVARS_FREQS_PATH: &str = "/mnt/data/mehari/0.21.0/db/grch37/seqvars/freqs";
const ROCKSDB_SEQVARS_CLINVAR_PATH: &str = "/mnt/data/mehari/0.21.0/db/grch37/seqvars/clinvar";
const TXSEQ_DB_PATH: &str = "/mnt/data/mehari/0.21.0/db/grch37/txs.bin.zst";

type DbCol<'a> = Arc<BoundColumnFamily<'a>>;

struct AnnotationSources {
    pub db_freq: DBWithThreadMode<MultiThreaded>,
    pub db_clinvar: DBWithThreadMode<MultiThreaded>,
    pub assembly: Assembly,
    pub predictor: ConsequencePredictor,
}

fn records(n: usize) -> Vec<Record> {
    let mut reader = noodles::vcf::reader::Builder::default()
        .build_from_path("tests/data/annotate/seqvars/NA-12878WGS_dragen.first10k.vcf.gz")
        .unwrap();
    let header = reader.read_header().unwrap();
    let records: Vec<_> = reader
        .records(&header)
        .flatten()
        .filter(|r| r.alternate_bases().len() == 1)
        .filter(|r| {
            let vcf_var = keys::Var::from_vcf_allele(r, 0);
            vcf_var.alternative != "*"
        })
        .take(n)
        .collect();
    records
}

fn seqvar_annotation(c: &mut Criterion) {
    let records = records(100);

    fn freq_db() -> Result<DBWithThreadMode<MultiThreaded>> {
        let options = rocksdb::Options::default();
        let db_freq = rocksdb::DB::open_cf_for_read_only(
            &options,
            ROCKSDB_SEQVARS_FREQS_PATH,
            ["meta", "autosomal", "gonosomal", "mitochondrial"],
            false,
        )?;

        Ok(db_freq)
    }

    fn clinvar_db() -> Result<DBWithThreadMode<MultiThreaded>> {
        let options = rocksdb::Options::default();
        let db_clinvar = rocksdb::DB::open_cf_for_read_only(
            &options,
            ROCKSDB_SEQVARS_CLINVAR_PATH,
            ["meta", "clinvar"],
            false,
        )?;
        Ok(db_clinvar)
    }

    fn transcript_db() -> Result<TxSeqDatabase> {
        load_tx_db(TXSEQ_DB_PATH)
    }

    fn setup() -> Result<AnnotationSources> {
        let db_freq = freq_db()?;
        let db_clinvar = clinvar_db()?;
        let tx_db = transcript_db()?;
        let assembly = Assembly::Grch37;

        let provider = Arc::new(MehariProvider::new(
            tx_db,
            assembly,
            MehariProviderConfigBuilder::default()
                .transcript_picking(true)
                .build()
                .unwrap(),
        ));
        let predictor = ConsequencePredictor::new(
            provider,
            assembly,
            ConsequencePredictorConfigBuilder::default()
                .report_all_transcripts(true)
                .transcript_source(TranscriptSource::Both)
                .build()
                .unwrap(),
        );
        Ok(AnnotationSources {
            db_freq,
            db_clinvar,
            assembly,
            predictor,
        })
    }

    fn handles(dbs: &AnnotationSources) -> (DbCol, DbCol, DbCol, DbCol) {
        let cf_autosomal = dbs.db_freq.cf_handle("autosomal").unwrap();
        let cf_gonosomal = dbs.db_freq.cf_handle("gonosomal").unwrap();
        let cf_mtdna = dbs.db_freq.cf_handle("mitochondrial").unwrap();
        let cf_clinvar = dbs.db_clinvar.cf_handle("clinvar").unwrap();
        (cf_autosomal, cf_gonosomal, cf_mtdna, cf_clinvar)
    }

    fn annotate_record(
        mut vcf_record: &mut Record,
        dbs: &AnnotationSources,
        cf_autosomal: &DbCol,
        cf_gonosomal: &DbCol,
        cf_mtdna: &DbCol,
        cf_clinvar: &DbCol,
    ) -> Result<()> {
        let db_freq = &dbs.db_freq;

        // We currently can only process records with one alternate allele.
        if vcf_record.alternate_bases().len() != 1 {
            // skip for now
            return Ok(());
        }

        // Get first alternate allele record.
        let vcf_var = keys::Var::from_vcf_allele(vcf_record, 0);

        // Skip records with a deletion as alternative allele.
        if vcf_var.alternative == "*" {
            // skip for now
            return Ok(());
        }

        // Only attempt lookups into RocksDB for canonical contigs.
        if is_canonical(vcf_var.chrom.as_str()) {
            // Build key for RocksDB database from `vcf_var`.
            let key: Vec<u8> = vcf_var.clone().into();

            // Annotate with frequency.
            if CHROM_AUTO.contains(vcf_var.chrom.as_str()) {
                annotate_record_auto(db_freq, cf_autosomal, &key, vcf_record)?;
            } else if CHROM_XY.contains(vcf_var.chrom.as_str()) {
                annotate_record_xy(db_freq, cf_gonosomal, &key, vcf_record)?;
            } else if CHROM_MT.contains(vcf_var.chrom.as_str()) {
                annotate_record_mt(db_freq, cf_mtdna, &key, vcf_record)?;
            } else {
                tracing::trace!(
                    "Record @{:?} on non-canonical chromosome, skipping.",
                    &vcf_var
                );
            }

            // Annotate with ClinVar information.
            annotate_record_clinvar(&dbs.db_clinvar, cf_clinvar, &key, vcf_record)?;
        }

        let keys::Var {
            chrom,
            pos,
            reference,
            alternative,
        } = vcf_var;

        // Annotate with variant effect.
        if let Some(ann_fields) = dbs.predictor.predict(&VcfVariant {
            chromosome: chrom,
            position: pos,
            reference,
            alternative,
        })? {
            if !ann_fields.is_empty() {
                vcf_record.info_mut().insert(
                    field::Key::from_str("ANN").unwrap(),
                    Some(field::Value::Array(field::value::Array::String(
                        ann_fields.iter().map(|ann| Some(ann.to_string())).collect(),
                    ))),
                );
            }
        }
        Ok(())
    }

    fn annotate_records(
        records: &mut Vec<Record>,
        dbs: &AnnotationSources,
        cf_autosomal: &DbCol,
        cf_gonosomal: &DbCol,
        cf_mtdna: &DbCol,
        cf_clinvar: &DbCol,
    ) -> Result<()> {
        for record in records {
            annotate_record(
                record,
                dbs,
                cf_autosomal,
                cf_gonosomal,
                cf_mtdna,
                cf_clinvar,
            )?
        }
        Ok(())
    }

    let dbs = setup().unwrap();
    let (cf_autosomal, cf_gonosomal, cf_mtdna, cf_clinvar) = handles(&dbs);
    let mut group = c.benchmark_group("annotate-records");
    group.bench_function("annotate-record", |b| {
        b.iter_batched(
            || records.clone(),
            |mut records| {
                annotate_records(
                    &mut records,
                    &dbs,
                    &cf_autosomal,
                    &cf_gonosomal,
                    &cf_mtdna,
                    &cf_clinvar,
                )
            },
            BatchSize::SmallInput,
        );
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default().with_profiler(perf::FlamegraphProfiler::new(100));
    targets = seqvar_annotation
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
