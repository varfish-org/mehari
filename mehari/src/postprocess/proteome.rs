use crate::annotate::seqvars::ann::{ANN_AA_SEQ_ALT, ANN_COMPOUND_IDS, ANN_COMPOUND_VARIANTS};
use crate::common::noodles::NoodlesVariantReader;
use anyhow::{Context, Error};
use clap::Args as ClapArgs;
use noodles::fasta;
use noodles::vcf;
use noodles::vcf::variant::record_buf::info::field::Value;
use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;

#[derive(Debug, ClapArgs)]
pub struct Args {
    /// Path to the Mehari-annotated VCF file.
    #[arg(long, required = true)]
    pub path_input_vcf: String,

    /// Path to the output FASTA file.
    #[arg(long, required = true)]
    pub path_output_fasta: String,
}

pub async fn run(_common: &crate::common::Args, args: &Args) -> Result<(), Error> {
    tracing::info!("Opening VCF to extract proteome: {}", args.path_input_vcf);

    let mut reader = crate::common::noodles::open_variant_reader(&args.path_input_vcf).await?;
    let header = reader.read_header().await?;

    let ann_header = header
        .infos()
        .get("ANN")
        .context("VCF does not contain an INFO/ANN header. Was it processed by Mehari?")?;

    let description = ann_header.description();

    let format_string = description
        .split('\'')
        .nth(1)
        .context("Could not parse ANN format from VCF header description.")?;

    let columns: Vec<&str> = format_string.split(" | ").map(|s| s.trim()).collect();

    let aa_seq_idx = columns.iter().position(|&c| c == ANN_AA_SEQ_ALT)
        .context(format!("VCF was not annotated with {ANN_AA_SEQ_ALT}. Please run 'annotate seqvars' with '--report-protein-sequence alternative' (or both) first."))?;
    let feature_id_idx = columns.iter().position(|&c| c == "Feature_ID").unwrap_or(6);
    let gene_name_idx = columns.iter().position(|&c| c == "Gene_Name").unwrap_or(3);
    let hgvs_p_idx = columns.iter().position(|&c| c == "HGVS.p").unwrap_or(10);

    let group_id_idx = columns
        .iter()
        .position(|&c| c == ANN_COMPOUND_IDS)
        .unwrap_or(usize::MAX);
    let group_vars_idx = columns
        .iter()
        .position(|&c| c == ANN_COMPOUND_VARIANTS)
        .unwrap_or(usize::MAX);

    let out_file = File::create(&args.path_output_fasta)?;
    let mut fasta_writer = fasta::io::Writer::new(BufWriter::new(out_file));
    let mut seen_headers = HashSet::new();
    let mut sequences_written = 0;

    let mut records = reader.records(&header).await;
    use futures::TryStreamExt;

    while let Some(record) = records.try_next().await? {
        if let Some(Some(Value::Array(
            vcf::variant::record_buf::info::field::value::Array::String(ann_array),
        ))) = record.info().get("ANN")
        {
            for ann_str in ann_array.iter().flatten() {
                let fields: Vec<&str> = ann_str.split('|').collect();

                if fields.len() > aa_seq_idx {
                    let aa_seq = fields[aa_seq_idx].trim();

                    if !aa_seq.is_empty() && aa_seq != "." {
                        let feature_id = fields.get(feature_id_idx).unwrap_or(&"UnknownFeature");
                        let gene_name = fields.get(gene_name_idx).unwrap_or(&"UnknownGene");
                        let hgvs_p = fields.get(hgvs_p_idx).unwrap_or(&"");

                        let group_id = if group_id_idx < fields.len() {
                            fields[group_id_idx].trim()
                        } else {
                            "Single"
                        };
                        let group_vars = if group_vars_idx < fields.len() {
                            fields[group_vars_idx].trim()
                        } else {
                            "unknown_vars"
                        };

                        let header_string = format!(
                            "{}|{} variant={} group={} vars={}",
                            gene_name, feature_id, hgvs_p, group_id, group_vars
                        );

                        if seen_headers.insert(header_string.clone()) {
                            let definition = fasta::record::Definition::new(
                                format!("{}|{}", gene_name, feature_id),
                                Some(
                                    format!(
                                        "variant={} group={} vars={}",
                                        hgvs_p, group_id, group_vars
                                    )
                                    .into(),
                                ),
                            );

                            let sequence =
                                fasta::record::Sequence::from(aa_seq.as_bytes().to_vec());
                            let fasta_record = fasta::Record::new(definition, sequence);

                            fasta_writer.write_record(&fasta_record)?;
                            sequences_written += 1;
                        }
                    }
                }
            }
        }
    }

    tracing::info!(
        "Successfully extracted {} unique mutated protein sequences to {}",
        sequences_written,
        args.path_output_fasta
    );
    Ok(())
}
