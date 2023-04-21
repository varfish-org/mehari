//! Copying of RocksDB database.

use std::{path::PathBuf, time::Instant};

use clap::Parser;

use crate::db::create::seqvar_clinvar::rocksdb_tuning;

/// Command line arguments for `db copy` sub command.
#[derive(Parser, Debug)]
#[command(about = "Copy rocksdb databases (parameter normalization)", long_about = None)]
pub struct Args {
    /// Path to input directory.
    #[arg(long)]
    pub path_in: PathBuf,
    /// Path to output directory.
    #[arg(long)]
    pub path_out: PathBuf,
}

/// Main entry point for `db copy` sub command.
pub fn run(common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!(
        "Copying RocksDB database\ncommon args: {:#?}\nargs: {:#?}",
        common,
        args
    );

    tracing::info!("Opening input database");
    // List all column families in database.
    let cf_names = rocksdb::DB::list_cf(&rocksdb::Options::default(), &args.path_in)?;
    let db_read = rocksdb::DB::open_cf_for_read_only(
        &rocksdb::Options::default(),
        &args.path_in,
        &cf_names,
        false,
    )?;

    // Get all column families from `db_read`.
    let cfs_read = cf_names
        .iter()
        .map(|cf| db_read.cf_handle(cf).unwrap())
        .collect::<Vec<_>>();

    tracing::info!("Opening output database");
    let mut options = rocksdb_tuning(rocksdb::Options::default());
    options.create_if_missing(true);
    options.create_missing_column_families(true);
    options.prepare_for_bulk_load();
    options.set_disable_auto_compactions(true);
    let db_write = rocksdb::DB::open_cf(&options, &args.path_out, &cf_names)?;

    // Get all column families from `db_write`.
    let cfs_write = cf_names
        .iter()
        .map(|cf| db_write.cf_handle(cf).unwrap())
        .collect::<Vec<_>>();

    // Perform the main work of copying over data.
    tracing::info!("Copying data");
    for (cf_name, (cf_read, cf_write)) in cf_names.iter().zip(cfs_read.iter().zip(cfs_write.iter()))
    {
        tracing::info!("Copying data from column family {}", cf_name);
        let mut iter = db_read.iterator_cf(cf_read, rocksdb::IteratorMode::Start);
        while let Some((key, value)) = iter.next().transpose()? {
            db_write.put_cf(cf_write, key, value)?;
        }
    }

    // Finally, compact manually.
    tracing::info!("Enforcing manual compaction");
    cfs_write
        .iter()
        .for_each(|cf| db_write.compact_range_cf(cf, None::<&[u8]>, None::<&[u8]>));

    let compaction_start = Instant::now();
    let mut last_printed = compaction_start;
    while db_write
        .property_int_value(rocksdb::properties::COMPACTION_PENDING)?
        .unwrap()
        > 0
        || db_write
            .property_int_value(rocksdb::properties::NUM_RUNNING_COMPACTIONS)?
            .unwrap()
            > 0
    {
        std::thread::sleep(std::time::Duration::from_millis(100));
        if last_printed.elapsed() > std::time::Duration::from_millis(1000) {
            log::info!(
                "... waiting for compaction for {:?}",
                compaction_start.elapsed()
            );
            last_printed = Instant::now();
        }
    }

    tracing::info!("Done copying RocksDB database");
    Ok(())
}
