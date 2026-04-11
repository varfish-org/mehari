use clap::Parser;
use std::path::PathBuf;

/// Command line arguments for `db create txs` sub command.
#[derive(Parser, Debug)]
#[command(about = "Construct mehari transcripts and sequence database", long_about = None)]
pub struct Args {
    /// Targeted genome assembly to extract transcripts for.
    #[arg(long)]
    pub assembly: String,

    /// Version of the genome assembly, e.g. "GRCh37.p13".
    #[arg(long)]
    pub assembly_version: Option<String>,

    /// Paths to the transcript annotations to import (cdot JSON or arbitrary GFF3).
    #[arg(long, required = true)]
    pub annotation: Vec<PathBuf>,

    /// Version of annotation data (if applicable).
    #[arg(long)]
    pub annotation_version: Option<String>,

    /// Source of the transcripts. For example "RefSeq" or "Ensembl".
    #[arg(long)]
    pub transcript_source: String,

    /// Version of the transcript source. E.g. "112" for Ensembl.
    #[arg(long, required_if_eq("transcript_source", "ensembl"))]
    pub transcript_source_version: Option<String>,

    /// Path to the seqrepo instance directory to use.
    #[arg(long, required_unless_present = "transcript_sequences")]
    pub seqrepo: Option<PathBuf>,

    /// Path to FASTA file(s) containing transcript sequences (alternative to seqrepo).
    #[arg(long, required_unless_present = "seqrepo")]
    pub transcript_sequences: Option<PathBuf>,

    /// Path to TSV file for label transfer of transcripts.  Columns are
    /// transcript id (without version), (unused) gene symbol, and label.
    #[arg(long)]
    pub mane_transcripts: Option<PathBuf>,

    /// Disable rigorous filtering (useful for custom annotations).
    #[arg(long, default_value = "false")]
    pub disable_filters: bool,

    /// Number of threads to use for steps supporting parallel processing.
    #[arg(long, default_value = "1")]
    pub threads: usize,

    /// ZSTD compression level to use.
    #[arg(long, default_value = "19")]
    pub compression_level: i32,

    /// Path to output protobuf file to write to.
    #[arg(long)]
    pub output: PathBuf,
}
