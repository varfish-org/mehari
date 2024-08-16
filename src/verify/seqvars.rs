//! Verification of the sequence variant consequence prediction.

use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
    sync::Arc,
    time::Instant,
};

use crate::annotate::seqvars::{
    csq::{ConfigBuilder as ConsequencePredictorConfigBuilder, ConsequencePredictor, VcfVariant},
    load_tx_db, path_component,
    provider::{ConfigBuilder as MehariProviderConfigBuilder, Provider as MehariProvider},
    TranscriptPickMode, TranscriptPickType,
};
use biocommons_bioutils::assemblies::Assembly;
use clap::Parser;
use noodles::core::{Position, Region};
use quick_cache::unsync::Cache;

/// Command line arguments for `verify seqvars` sub command.
#[derive(Parser, Debug)]
#[command(about = "Compare variant effect predictions to VEP ones", long_about = None)]
pub struct Args {
    /// Path to the mehari database folder.
    #[arg(long)]
    pub path_db: String,

    /// Path to the input TSV file.
    #[arg(long)]
    pub path_input_tsv: String,
    /// Path to the reference FASTA file.

    #[arg(long)]
    pub path_reference_fasta: String,
    /// Path to output TSV file.

    #[arg(long)]
    pub path_output_tsv: String,

    /// Whether to report only the worst consequence for each picked transcript.
    #[arg(long, default_value_t = false)]
    pub report_most_severe_consequence_only: bool,

    /// Which kind of transcript to pick / restrict to. Default is not to pick at all.
    ///
    /// Depending on `--pick-transcript-mode`, if multiple transcripts match the selection,
    /// either the first one is kept or all are kept.
    #[arg(long)]
    pub pick_transcript: Vec<TranscriptPickType>,

    /// Determines how to handle multiple transcripts. Default is to keep all.
    ///
    /// When transcript picking is enabled via `--pick-transcript`,
    /// either keep the first one found or keep all that match.
    #[arg(long, default_value = "all")]
    pub pick_transcript_mode: TranscriptPickMode,

    /// For debug purposes, maximal number of variants to annotate.
    #[arg(long)]
    pub max_var_count: Option<usize>,
}

/// Guess genome release from VEP TSV file.
fn guess_assembly(path_input_tsv: &str) -> Result<Assembly, anyhow::Error> {
    tracing::info!("Guessing assembly from {}...", &path_input_tsv);
    let mut result = None;
    let lines = BufReader::new(File::open(path_input_tsv)?).lines();
    for line in lines {
        let line = line?;
        if line.starts_with("## assembly version") {
            let token = line
                .split_whitespace()
                .last()
                .expect("problem splitting 'assembly version' line");
            if token.starts_with("GRCh37") {
                result = Some(Assembly::Grch37);
                break;
            } else if token.starts_with("GRCh38") {
                result = Some(Assembly::Grch38);
                break;
            } else {
                anyhow::bail!("unknown genome release: {}", token);
            }
        } else if !line.starts_with("##") {
            break;
        }
    }

    if let Some(assembly) = result {
        tracing::info!("... guessed assembly to be: {:?}", assembly);
        Ok(assembly)
    } else {
        anyhow::bail!("could not guess assembly {}", path_input_tsv);
    }
}

/// Reading of VEP TSV files.
pub mod vep_tsv {
    use serde::Deserialize;

    /// Structure for deserializing VEP TSV line.
    #[derive(Debug, Deserialize)]
    pub struct VepRecord {
        pub uploaded_variation: String,
        pub location: String,
        pub allele: String,
        pub gene: String,
        pub feature: String,
        pub feature_type: String,
        pub consequence: String,
        pub cdna_position: String,
        pub cds_position: String,
        pub protein_position: String,
        pub amino_acids: String,
        pub codons: String,
        pub existing_variation: String,
        pub extra: String,
    }
}

/// Run the verification command.
pub fn run(_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    // Guess assembly from VEP TSV file.
    let assembly = guess_assembly(&args.path_input_tsv)?;

    // Output TSV file.
    let mut output_tsv = std::fs::File::create(&args.path_output_tsv)?;

    // Open the reference FASTA through a noodles FAI reader.
    tracing::info!(
        "Opening reference FASTA file: {}",
        &args.path_reference_fasta
    );
    let fai_index = noodles::fasta::fai::read(format!("{}.fai", args.path_reference_fasta))?;
    let mut fai_reader = noodles::fasta::indexed_reader::Builder::default()
        .set_index(fai_index)
        .build_from_path(&args.path_reference_fasta)?;

    // Read the serialized transcripts.
    tracing::info!("Opening transcript database");
    let tx_db = load_tx_db(format!(
        "{}/{}/txs.bin.zst",
        &args.path_db,
        path_component(assembly)
    ))?;

    let provider = Arc::new(MehariProvider::new(
        tx_db,
        assembly,
        MehariProviderConfigBuilder::default()
            .pick_transcript(args.pick_transcript.clone())
            .pick_transcript_mode(args.pick_transcript_mode)
            .build()?,
    ));

    let predictor = ConsequencePredictor::new(
        provider,
        assembly,
        ConsequencePredictorConfigBuilder::default()
            .report_most_severe_consequence_only(args.report_most_severe_consequence_only)
            .build()?,
    );

    // LRU caches used below to avoid re-reading from FAI and prediction.
    let mut ref_cache = Cache::new(100);
    let mut pred_cache = Cache::new(100);

    // Print header.
    writeln!(
        &mut output_tsv,
        "result\tlocation\tvep_feature\tvep_consequences\tmehari_feature\tmehari_consequence"
    )?;

    // Read through the VEP TSV file and compare to the mehari predictions.
    tracing::info!("Processing input TSV file ...");
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .comment(Some(b'#'))
        .from_path(&args.path_input_tsv)?;
    let mut prev = Instant::now();
    for (i, record) in reader.deserialize().enumerate() {
        let record: vep_tsv::VepRecord = record?;
        let record = {
            let mut record = record;
            record.consequence = record.consequence.replace(",NMD_transcript_variant", "");

            record
        };

        if prev.elapsed().as_secs() >= 5 {
            tracing::info!("at {:?}", &record.location);
            prev = Instant::now();
        }

        // Exract `USED_REF=` field from `record.extra`.
        let used_ref = record
            .extra
            .split(';')
            .find_map(|x| x.strip_prefix("USED_REF="));

        // Skip if not on an ENST transcript.
        if !record.feature.starts_with("ENST") {
            continue;
        }

        // Split the location into contig and position.
        let mut tokens = record.location.split(':');
        let contig = tokens.next().expect("problem splitting location");
        let pos = tokens.next().expect("problem splitting location");
        // Extract start and end position (if any) or fall back to start == end if no range.
        let mut pos_tokens = pos.split('-');
        let start = Position::try_from(
            pos_tokens
                .next()
                .expect("problem splitting position")
                .parse::<usize>()?,
        )?;
        let end = Position::try_from(
            pos_tokens
                .next()
                .map(|x| x.parse::<usize>())
                .transpose()?
                .unwrap_or(start.into()),
        )?;

        // The VEP TSV encoding is a bit cumbersome.
        //
        // In the case of deletions, the `record.allele` field contains `"-"`.  In the case of
        // insertions, `used_ref` will be `Some("-")`.  In both cases, we need to shift the
        // start position one to the left and possibly adjust the end position as well.
        let (is_del, is_ins, start, end) = if record.allele == "-" {
            // Is deletion, need to shift start position to the left by one.
            // Further down, we will need to use the first base of the reference
            // allele as the alternate allele. TODO
            (true, false, Position::new(start.get() - 1).unwrap(), end)
        } else if used_ref == Some("-") {
            // Further down, we will need to expand the alternate allele to the
            // left by the single reference allele base.
            (false, true, start, start)
        } else {
            (false, false, start, end)
        };

        // Extract the reference allele from the FASTA file (load from cache if possible).
        let reference_allele = if let Some(reference_allele) =
            ref_cache.get(&(record.location.clone(), is_del, is_ins))
        {
            reference_allele
        } else {
            let reference_allele = std::str::from_utf8(
                fai_reader
                    .query(&Region::new(contig, start..=end))?
                    .sequence()
                    .as_ref(),
            )?
            .to_owned();
            ref_cache.insert((record.location.clone(), is_del, is_ins), reference_allele);
            ref_cache
                .get(&(record.location.clone(), is_del, is_ins))
                .unwrap()
        };

        // Determine alternate allele, deletions and insertions need special handling.
        let alt_allele = if is_del {
            reference_allele.chars().next().unwrap().to_string()
        } else if is_ins {
            format!("{}{}", &reference_allele, &record.allele)
        } else {
            record.allele.clone()
        };

        // Extract prediction from mehari (load from cache if possible).
        let anns = if let Some(pred_cache) = pred_cache.get(&record.location) {
            pred_cache
        } else {
            let vcf_var = VcfVariant {
                chromosome: contig.to_owned(),
                position: start.get() as i32,
                reference: reference_allele.to_owned(),
                alternative: alt_allele.clone(),
            };
            let anns = predictor.predict(&vcf_var)?;
            pred_cache.insert(record.location.clone(), anns);
            pred_cache.get(&record.location).unwrap()
        };

        if let Some(anns) = anns {
            let mut found_any = false;
            for ann in anns {
                if ann.feature_id.starts_with(&record.feature) {
                    let ann_consequences = ann
                        .consequences
                        .iter()
                        .map(|s| s.to_string())
                        .collect::<Vec<_>>()
                        .join(",");
                    let result = if ann_consequences == record.consequence {
                        "OK"
                    } else {
                        "mismatch"
                    };
                    writeln!(
                        &mut output_tsv,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        result,
                        &record.location,
                        &record.feature,
                        &record.consequence,
                        &ann.feature_id,
                        &ann_consequences
                    )?;

                    found_any = true;
                    break;
                }
            }
            if !found_any {
                // There are annotations for this variant but there is no match.
                writeln!(
                    &mut output_tsv,
                    "mehari_no_match\t{}\t{}\t{}\t.\t.",
                    &record.location, &record.feature, &record.consequence
                )?;
            }
        } else {
            // There are no annotations for this variant.
            writeln!(
                &mut output_tsv,
                "mehari_no_tx\t{}\t{}\t{}\t.\t.",
                &record.location, &record.feature, &record.consequence
            )?;
        }

        // Break after `args.max_var_count` if provided.
        if let Some(max_var_count) = args.max_var_count {
            if i + 1 >= max_var_count {
                break;
            }
        }
    }
    tracing::info!("... done processing input TSV");

    Ok(())
}
