//! Annotation of sequence variants.

use std::fs::File;
use std::io::{BufRead, Write};
use std::str::FromStr;
use std::time::Instant;

use clap::Parser;
use noodles::bgzf::Writer as BgzfWriter;
use noodles::vcf::header::{
    record::value::map::{info::Type, Info},
    Number,
};
use noodles::vcf::record::info::field::Value;
use noodles::vcf::{
    header::info::key::Key as InfoKey, header::info::key::Other as InfoKeyOther,
    header::record::value::map::Map, Header as VcfHeader, Writer as VcfWriter,
};
use noodles_util::variant::reader::Builder as VariantReaderBuilder;
use rocksdb::{DBWithThreadMode, ThreadMode};

use crate::db::create::seqvar_freqs::serialized::auto::Record as AutoRecord;
use crate::db::create::seqvar_freqs::serialized::vcf::Var as VcfVar;

/// Command line arguments for `annotate seqvars` sub command.
#[derive(Parser, Debug)]
#[command(about = "Annotate sequence variant VCF files", long_about = None)]
pub struct Args {
    /// Path to the mehari database folder.
    #[arg(long)]
    pub path_db: String,

    /// Path to the input VCF file.
    #[arg(long)]
    pub path_input_vcf: String,
    /// Path to the output VCF file.
    #[arg(long)]
    pub path_output_vcf: String,

    /// For debug purposes, maximal number of variants to annotate.
    #[arg(long)]
    pub max_var_count: Option<usize>,
}

pub mod keys {
    use std::str::FromStr;

    use noodles::vcf::{
        header::info::key::Key as InfoKey, header::info::key::Other as InfoKeyOther,
    };

    lazy_static::lazy_static! {
        pub static ref GNOMAD_EXOMES_AN: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_exomes_an").unwrap());
        pub static ref GNOMAD_EXOMES_HOM: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_exomes_hom").unwrap());
        pub static ref GNOMAD_EXOMES_HET: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_exomes_het").unwrap());
        pub static ref GNOMAD_EXOMES_HEMI: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_exomes_hemi").unwrap());

        pub static ref GNOMAD_GENOMES_AN: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_genomes_an").unwrap());
        pub static ref GNOMAD_GENOMES_HOM: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_genomes_hom").unwrap());
        pub static ref GNOMAD_GENOMES_HET: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_genomes_het").unwrap());
        pub static ref GNOMAD_GENOMES_HEMI: InfoKey = InfoKey::Other(InfoKeyOther::from_str("gnomad_genomes_hemi").unwrap());

        pub static ref HELIX_AN: InfoKey = InfoKey::Other(InfoKeyOther::from_str("helix_an").unwrap());
        pub static ref HELIX_HOM: InfoKey = InfoKey::Other(InfoKeyOther::from_str("helix_hom").unwrap());
        pub static ref HELIX_HET: InfoKey = InfoKey::Other(InfoKeyOther::from_str("helix_het").unwrap());
    }
}

fn build_header(header_in: &VcfHeader) -> VcfHeader {
    let mut header_out = header_in.clone();

    header_out.infos_mut().insert(
        keys::GNOMAD_EXOMES_AN.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of samples in gnomAD exomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_EXOMES_HOM.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of hom. alt. carriers in gnomAD exomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_EXOMES_HET.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of het. alt. carriers in gnomAD exomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_EXOMES_HEMI.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of hemi. alt. carriers in gnomAD exomes",
        ),
    );

    header_out.infos_mut().insert(
        keys::GNOMAD_GENOMES_AN.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of samples in gnomAD genomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_GENOMES_HOM.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of hom. alt. carriers in gnomAD genomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_GENOMES_HET.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of het. alt. carriers in gnomAD genomes",
        ),
    );
    header_out.infos_mut().insert(
        keys::GNOMAD_GENOMES_HEMI.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of hemi. alt. carriers in gnomAD genomes",
        ),
    );

    header_out.infos_mut().insert(
        keys::HELIX_AN.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of samples in HelixMtDb",
        ),
    );
    header_out.infos_mut().insert(
        keys::HELIX_HOM.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of hom. alt. carriers in HelixMtDb",
        ),
    );
    header_out.infos_mut().insert(
        keys::HELIX_HET.clone(),
        Map::<Info>::new(
            Number::Count(1),
            Type::Integer,
            "Number of het. alt. carriers in HelixMtDb",
        ),
    );
    header_out
}

fn run_with_writer<Inner: Write>(
    mut writer: VcfWriter<Inner>,
    args: &Args,
) -> Result<(), anyhow::Error> {
    // Open the database in read only mode.
    tracing::info!("Opening database");
    let options = rocksdb::Options::default();
    let db = rocksdb::DB::open_cf_for_read_only(
        &options,
        &args.path_db,
        ["meta", "autosomal", "gonosomal", "mitochondrial"],
        false,
    )?;

    let cf_autosomal = db.cf_handle("autosomal").unwrap();
    let cf_gonosomal = db.cf_handle("gonosomal").unwrap();
    let cf_mtdna = db.cf_handle("mitochondrial").unwrap();

    tracing::info!("Annotating VCF ...");
    let start = Instant::now();

    let mut reader = VariantReaderBuilder::default().build_from_path(&args.path_input_vcf)?;
    let header_in = reader.read_header()?;
    let header_out = build_header(&header_in);

    let mut total_written = 0usize;

    writer.write_header(&header_out)?;
    let mut records = reader.records(&header_in);
    loop {
        if let Some(record) = records.next() {
            let mut vcf_record = record?;
            let vcf_var = VcfVar::from_vcf(&vcf_record);
            let key: Vec<u8> = vcf_var.into();

            if let Some(freq) = db.get_cf(cf_autosomal, key)? {
                let auto_record = AutoRecord::from_buf(&freq);

                vcf_record.info_mut().insert(
                    keys::GNOMAD_EXOMES_AN.clone(),
                    Some(Value::Integer(auto_record.gnomad_exomes.an as i32)),
                );
                vcf_record.info_mut().insert(
                    keys::GNOMAD_EXOMES_HOM.clone(),
                    Some(Value::Integer(auto_record.gnomad_exomes.ac_hom as i32)),
                );
                vcf_record.info_mut().insert(
                    keys::GNOMAD_EXOMES_HET.clone(),
                    Some(Value::Integer(auto_record.gnomad_exomes.ac_het as i32)),
                );

                vcf_record.info_mut().insert(
                    keys::GNOMAD_GENOMES_AN.clone(),
                    Some(Value::Integer(auto_record.gnomad_genomes.an as i32)),
                );
                vcf_record.info_mut().insert(
                    keys::GNOMAD_GENOMES_HOM.clone(),
                    Some(Value::Integer(auto_record.gnomad_genomes.ac_hom as i32)),
                );
                vcf_record.info_mut().insert(
                    keys::GNOMAD_GENOMES_HET.clone(),
                    Some(Value::Integer(auto_record.gnomad_genomes.ac_het as i32)),
                );

                println!("{:?}", &vcf_record);
            };

            writer.write_record(&vcf_record)?;
        } else {
            break; // all done
        }

        total_written += 1;
        if let Some(max_var_count) = args.max_var_count {
            if total_written >= max_var_count {
                tracing::warn!(
                    "Stopping after {} records as requested by --max-var-count",
                    total_written
                );
                break;
            }
        }
    }
    tracing::info!(
        "... annotated {} records in {:?}",
        total_written,
        start.elapsed()
    );
    Ok(())
}

/// Main entry point for `annotate seqvars` sub command.
pub fn run(_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    if args.path_output_vcf.ends_with(".vcf.gz") || args.path_output_vcf.ends_with(".vcf.bgzf") {
        let writer = VcfWriter::new(File::create(&args.path_output_vcf).map(BgzfWriter::new)?);
        run_with_writer(writer, args)?;
    } else {
        let writer = VcfWriter::new(File::create(&args.path_output_vcf)?);
        run_with_writer(writer, args)?;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use clap_verbosity_flag::Verbosity;
    use pretty_assertions::assert_eq;
    use temp_testdir::TempDir;

    use super::{run, Args};

    #[test]
    fn smoke_test() -> Result<(), anyhow::Error> {
        let temp = TempDir::default();
        let path_out = temp.join("output.vcf");

        let args_common = crate::common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            path_db: String::from(
                "tests/data/db/create/seqvar_freqs/db-rs1263393206/seqvars/freqs",
            ),
            path_input_vcf: String::from(
                "tests/data/db/create/seqvar_freqs/db-rs1263393206/input.vcf",
            ),
            path_output_vcf: path_out.into_os_string().into_string().unwrap(),
            max_var_count: None,
        };

        run(&args_common, &args)?;

        let actual = std::fs::read_to_string(args.path_output_vcf)?;
        let expected = std::fs::read_to_string(
            "tests/data/db/create/seqvar_freqs/db-rs1263393206/output.vcf",
        )?;
        assert_eq!(&expected, &actual);

        Ok(())
    }
}
