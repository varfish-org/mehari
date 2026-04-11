use thiserror::Error;

#[derive(Debug, Error)]
pub enum MehariError {
    #[error(transparent)]
    Seqvars(#[from] SeqvarsError),

    #[error(transparent)]
    Strucvars(#[from] StrucvarsError),
}

#[derive(Debug, Error)]
pub enum SeqvarsError {
    #[error("Phased variant group validation failed: {0}")]
    GroupValidation(#[from] GroupValidationError),

    #[error("HGVS projection failed: {0}")]
    HgvsProjection(String),

    #[error("Could not determine chromosome accession")]
    UnknownChromosomeAccession,

    #[error("Invalid coordinates start={0}, end={1} (max is 512M)")]
    InvalidCoordinates(i32, i32),

    #[error(transparent)]
    HgvsMapper(#[from] hgvs::mapper::Error),

    #[error(transparent)]
    HgvsNormalizer(#[from] hgvs::normalizer::Error),

    #[error(transparent)]
    HgvsData(#[from] hgvs::data::error::Error),

    #[error("Provider error: {0}")]
    Provider(String),
}

#[derive(Debug, Error)]
pub enum GroupValidationError {
    #[error("All phased variants must be on the same chromosome.")]
    DifferentChromosomes,

    #[error("Overlapping variants detected at pos {0} and {1}.")]
    OverlappingVariants(i32, i32),

    #[error("Variant cluster spans too large a genomic region (>10kb).")]
    ClusterTooLarge,

    #[error("Multiple insertion-only variants at the same normalized position.")]
    MultipleInsertionsSamePosition,
}

#[derive(Debug, Error)]
pub enum StrucvarsError {
    #[error("Invalid PE orientation: {0}")]
    InvalidPeOrientation(String),

    #[error("Invalid SV type: {0}")]
    InvalidSvType(String),

    #[error("Invalid SV sub type: {0}")]
    InvalidSvSubType(String),

    #[error("Missing required SV field: {0}")]
    MissingSvField(String),
}
