use thiserror::Error;

#[derive(Debug, Error)]
pub enum MehariError {
    #[error(transparent)]
    Io(#[from] IoError),

    #[error(transparent)]
    Seqvars(#[from] SeqvarsError),

    #[error(transparent)]
    Strucvars(#[from] StrucvarsError),

    #[error(transparent)]
    Generic(#[from] anyhow::Error),
}
#[derive(Debug, Error)]
pub enum IoError {
    #[error("I/O error: {0}")]
    StdIo(#[from] std::io::Error),
    #[error("Failed to open file {path}: {msg}")]
    FileOpen { path: String, msg: String },
    #[error("Compression/Decompression error: {0}")]
    Compression(String),
}

#[derive(Debug, Error)]
pub enum SeqvarsError {
    #[error("Phased variant group validation failed: {0}")]
    GroupValidation(#[from] GroupValidationError),

    #[error("HGVS projection failed: {0}")]
    HgvsProjection(String),

    #[error("Could not determine chromosome accession")]
    UnknownChromosomeAccession,
}

#[derive(Debug, Error)]
pub enum GroupValidationError {
    #[error("All phased variants must be on the same chromosome.")]
    DifferentChromosomes,
    #[error("Overlapping variants detected at pos {0} and {1}.")]
    OverlappingVariants(i32, i32),
    #[error("Variant cluster spans too large a genomic region (>10kb).")]
    ClusterTooLarge,
}

#[derive(Debug, Error)]
pub enum StrucvarsError {
    #[error("Invalid PE orientation: {0}")]
    InvalidPeOrientation(String),
    #[error("Invalid SV type: {0}")]
    InvalidSvType(String),
    #[error("Missing required SV field: {0}")]
    MissingSvField(String),
}
