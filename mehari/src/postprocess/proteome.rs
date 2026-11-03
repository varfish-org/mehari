use crate::annotate::seqvars::consequence::terms::{
    ANN_AA_SEQ_ALT, ANN_COMPOUND_IDS, ANN_COMPOUND_VARIANTS, ANN_KEY, AnnField,
};
use crate::common::noodles::NoodlesVariantReader;
use anyhow::{Context, Error};
use clap::{Args as ClapArgs, ValueEnum};
use noodles::fasta;
use noodles::vcf;
use noodles::vcf::variant::record_buf::info::field::Value;
use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::BufWriter;

#[derive(Debug, Clone, ValueEnum, Default, PartialEq)]
pub enum DedupMode {
    /// Output every valid protein sequence.
    #[default]
    None,

    /// Output exactly one sequence per transcript (the one with the most variants).
    Transcript,

    /// Output all compound events, plus any "orphaned" single variants.
    Maximal,
}

#[derive(Debug, ClapArgs)]
pub struct Args {
    /// Path to the Mehari-annotated VCF file.
    #[arg(long, required = true)]
    pub path_input_vcf: String,

    /// Path to the output FASTA file.
    #[arg(long, required = true)]
    pub path_output_fasta: String,

    /// Deduplication strategy for the FASTA output.
    #[arg(long, value_enum, default_value_t = DedupMode::None)]
    pub deduplicate: DedupMode,

    /// INFO key to read the annotations from.
    #[arg(long, default_value = ANN_KEY)]
    pub ann_key: String,
}

struct ParsedEvent {
    is_compound: bool,
    var_ids: Vec<String>,
    header_string: String,
    record: fasta::Record,
}

pub async fn run(_common: &crate::common::Args, args: &Args) -> Result<(), Error> {
    tracing::info!("Opening VCF to extract proteome: {}", args.path_input_vcf);

    let mut reader = crate::common::noodles::open_variant_reader(&args.path_input_vcf).await?;
    let header = reader.read_header().await?;

    let ann_key = &args.ann_key;
    let columns = AnnField::ann_field_names_from_header(&header, ann_key).with_context(|| {
        format!("VCF has no field names in its INFO/{ann_key} header. Was it processed by Mehari?")
    })?;

    let aa_seq_idx = columns.iter().position(|c| c == ANN_AA_SEQ_ALT)
        .context(format!("VCF was not annotated with {ANN_AA_SEQ_ALT}. Please run 'annotate seqvars' with '--plugin aa-sequence-alt' first."))?;
    let feature_id_idx = columns
        .iter()
        .position(|c| c == "Feature_ID")
        .with_context(|| format!("Required column 'Feature_ID' not found in {ann_key} header"))?;
    let gene_name_idx = columns
        .iter()
        .position(|c| c == "Gene_Name")
        .with_context(|| format!("Required column 'Gene_Name' not found in {ann_key} header"))?;
    let hgvs_p_idx = columns
        .iter()
        .position(|c| c == "HGVS.p")
        .with_context(|| format!("Required column 'HGVS.p' not found in {ann_key} header"))?;

    let group_id_idx = columns
        .iter()
        .position(|c| c == ANN_COMPOUND_IDS)
        .unwrap_or(usize::MAX);
    let group_vars_idx = columns
        .iter()
        .position(|c| c == ANN_COMPOUND_VARIANTS)
        .unwrap_or(usize::MAX);

    let out_file = File::create(&args.path_output_fasta)?;
    let mut fasta_writer = fasta::io::Writer::new(BufWriter::new(out_file));

    let mut seen_headers = HashSet::new();
    let mut transcript_buffer: BTreeMap<String, Vec<ParsedEvent>> = BTreeMap::new();
    let mut sequences_written = 0;

    let mut records = reader.records(&header).await;
    use futures::TryStreamExt;

    while let Some(record) = records.try_next().await? {
        let native_var = crate::annotate::seqvars::from_vcf_allele(&record, 0);
        let native_var_id = format!(
            "{}:{}:{}:{}",
            native_var.chrom, native_var.pos, native_var.reference, native_var.alternative
        );

        if let Some(Some(Value::Array(
            vcf::variant::record_buf::info::field::value::Array::String(ann_array),
        ))) = record.info().get(ann_key)
        {
            for ann_str in ann_array.iter().flatten() {
                let fields: Vec<&str> = ann_str.split('|').collect();

                if fields.len() > aa_seq_idx {
                    let aa_seq = fields[aa_seq_idx].trim();

                    if !aa_seq.is_empty() && aa_seq != "." {
                        let feature_id = fields
                            .get(feature_id_idx)
                            .unwrap_or(&"UnknownFeature")
                            .to_string();
                        let gene_name = fields.get(gene_name_idx).unwrap_or(&"UnknownGene");
                        let hgvs_p = fields.get(hgvs_p_idx).unwrap_or(&"");

                        let group_id = if group_id_idx < fields.len()
                            && !fields[group_id_idx].trim().is_empty()
                        {
                            fields[group_id_idx].trim()
                        } else {
                            "Single"
                        };
                        let group_vars = if group_vars_idx < fields.len()
                            && !fields[group_vars_idx].trim().is_empty()
                        {
                            fields[group_vars_idx].trim()
                        } else {
                            &native_var_id
                        };

                        let is_compound = group_id != "Single";
                        let var_ids: Vec<String> =
                            group_vars.split(',').map(|s| s.to_string()).collect();

                        let header_string = format!(
                            "{}|{} variant={} group={} vars={}",
                            gene_name, feature_id, hgvs_p, group_id, group_vars
                        );
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

                        let sequence = fasta::record::Sequence::from(aa_seq.as_bytes().to_vec());
                        let fasta_record = fasta::Record::new(definition, sequence);

                        if args.deduplicate == DedupMode::None {
                            if seen_headers.insert(header_string) {
                                fasta_writer.write_record(&fasta_record)?;
                                sequences_written += 1;
                            }
                        } else {
                            transcript_buffer
                                .entry(feature_id)
                                .or_default()
                                .push(ParsedEvent {
                                    is_compound,
                                    var_ids,
                                    header_string,
                                    record: fasta_record,
                                });
                        }
                    }
                }
            }
        }
    }

    if args.deduplicate != DedupMode::None {
        for (_, mut events) in transcript_buffer {
            match args.deduplicate {
                DedupMode::Transcript => {
                    events.sort_by_key(|a| std::cmp::Reverse(a.var_ids.len()));

                    if let Some(best_event) = events.into_iter().next()
                        && seen_headers.insert(best_event.header_string)
                    {
                        fasta_writer.write_record(&best_event.record)?;
                        sequences_written += 1;
                    }
                }
                DedupMode::Maximal => {
                    let mut incorporated_variants = HashSet::new();
                    for event in &events {
                        if event.is_compound {
                            incorporated_variants.extend(event.var_ids.clone());
                        }
                    }

                    for event in events {
                        let should_write = if event.is_compound {
                            true
                        } else {
                            let single_id = &event.var_ids[0];
                            !incorporated_variants.contains(single_id)
                        };

                        if should_write && seen_headers.insert(event.header_string) {
                            fasta_writer.write_record(&event.record)?;
                            sequences_written += 1;
                        }
                    }
                }
                DedupMode::None => unreachable!(),
            }
        }
    }

    tracing::info!(
        "Successfully extracted {} mutated protein sequences to {} (Mode: {:?})",
        sequences_written,
        args.path_output_fasta,
        args.deduplicate
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use temp_testdir::TempDir;

    /// Only `MEHARI` has the protein sequence, so reading `ANN` finds none.
    #[rstest]
    #[case::chosen_key("MEHARI", true)]
    #[case::other_field("ANN", false)]
    #[tokio::test]
    async fn reads_annotations_from_ann_key(
        #[case] ann_key: &str,
        #[case] finds_sequence: bool,
    ) -> Result<(), Error> {
        let temp = TempDir::default();
        let path_in = temp.join("input.vcf");
        let path_out = temp.join("output.fasta");
        std::fs::write(
            &path_in,
            format!(
                "##fileformat=VCFv4.2\n\
                 ##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: \
                 'Allele | Annotation'\">\n\
                 ##INFO=<ID=MEHARI,Number=.,Type=String,Description=\"Functional annotations: \
                 'Feature_ID | Gene_Name | HGVS.p | {ANN_AA_SEQ_ALT}'\">\n\
                 ##contig=<ID=17,length=81195210>\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 17\t1000\t.\tA\tT\t.\t.\tANN=T|other_tool;MEHARI=NM_1|GENE|p.Ala2Thr|MTK\n"
            ),
        )?;

        let result = run(
            &crate::common::Args {
                verbose: Default::default(),
            },
            &Args {
                path_input_vcf: path_in.display().to_string(),
                path_output_fasta: path_out.display().to_string(),
                deduplicate: DedupMode::None,
                ann_key: ann_key.to_string(),
            },
        )
        .await;

        assert_eq!(result.is_ok(), finds_sequence, "{result:?}");
        if finds_sequence {
            assert_eq!(
                std::fs::read_to_string(&path_out)?,
                ">GENE|NM_1 variant=p.Ala2Thr group=Single vars=17:1000:A:T\nMTK\n"
            );
        }
        Ok(())
    }
}
