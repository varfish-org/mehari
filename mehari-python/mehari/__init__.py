import typing
from pathlib import Path

import polars as pl
import pyarrow as pa

from ._mehari import SeqvarsAnnotator as _SeqvarsAnnotator
from ._mehari import (
    consequence_variants,
    putative_impact_variants,
    feature_biotype_variants,
)
from ._mehari import build_transcript_db as _build_transcript_db

ConsequenceEnum = pl.Enum(consequence_variants())
ImpactEnum = pl.Enum(putative_impact_variants())
FeatureBiotypeEnum = pl.Enum(feature_biotype_variants())
FeatureBiotypeType = typing.Literal["Coding", "Noncoding"]


class VariantDict(typing.TypedDict):
    chromosome: str
    position: int
    reference: str
    alternative: str


class RankDict(typing.TypedDict):
    ord: int
    total: int


class PosDict(typing.TypedDict):
    ord: int
    total: int


class AnnotationDict(typing.TypedDict):
    allele: str
    consequences: list[str]
    putative_impact: str
    gene_symbol: str
    gene_id: str
    feature_type: str
    feature_id: str
    feature_biotype: list[FeatureBiotypeType]
    feature_tags: list[str]
    rank: RankDict | None
    hgvs_g: str | None
    hgvs_n: str | None
    hgvs_c: str | None
    hgvs_p: str | None
    cdna_pos: PosDict | None
    cds_pos: PosDict | None
    protein_pos: PosDict | None
    distance: int | None
    strand: int
    messages: list[str] | None
    custom_fields: dict[str, str | None] | None


class AnnotationResultDict(typing.TypedDict):
    annotation: list[AnnotationDict]


class SeqvarsAnnotator:
    def __init__(
        self,
        transcript_db_paths: str | list[str],
        reference_path: str | None = None,
        *,
        report_cdna_sequence: typing.Literal[
            "none", "reference", "alternative", "both"
        ] = "none",
        report_protein_sequence: typing.Literal[
            "none", "reference", "alternative", "both"
        ] = "none",
    ):
        if isinstance(transcript_db_paths, str):
            transcript_db_paths = [transcript_db_paths]

        self._annotator = _SeqvarsAnnotator(
            transcript_db_paths,
            reference_path,
            report_cdna_sequence,
            report_protein_sequence,
        )

    def _process_dataframe(self, df: pl.DataFrame) -> pl.DataFrame:
        """Internal helper to stream Arrow RecordBatches to Rust."""
        pa_table = df.to_arrow()
        batches = pa_table.to_batches()

        if not batches:
            batches = [pa.RecordBatch.from_pylist([], schema=pa_table.schema)]

        result_batches = [self._annotator.annotate_batch(b) for b in batches]
        res_df = pl.from_arrow(pa.Table.from_batches(result_batches))
        return df.with_columns(
            res_df.get_column("annotation")
            .list.eval(
                pl.element().struct.with_fields(
                    pl.element().struct.field("putative_impact").cast(ImpactEnum),
                    pl.element()
                    .struct.field("consequences")
                    .cast(pl.List(ConsequenceEnum)),
                    pl.element()
                    .struct.field("feature_biotype")
                    .cast(pl.List(FeatureBiotypeEnum)),
                )
            )
            .alias("annotation")
        )

    @typing.overload
    def annotate(self, data: str) -> AnnotationResultDict: ...

    @typing.overload
    def annotate(self, data: pl.DataFrame) -> pl.DataFrame: ...

    @typing.overload
    def annotate(self, data: pl.LazyFrame) -> pl.LazyFrame: ...

    @typing.overload
    def annotate(
        self, *, chromosome: str, position: int, reference: str, alternative: str
    ) -> AnnotationResultDict: ...

    def annotate(
        self,
        data: str | pl.DataFrame | pl.LazyFrame | None = None,
        *,
        chromosome: str | None = None,
        position: int | None = None,
        reference: str | None = None,
        alternative: str | None = None,
    ) -> typing.Any:
        """
        Annotates genetic variants.
        Supports Polars DataFrames/LazyFrames, 'chr:pos:ref:alt' strings, or explicit kwargs.
        """

        if isinstance(data, pl.DataFrame):
            return self._process_dataframe(data)

        if isinstance(data, pl.LazyFrame):
            in_schema = data.collect_schema()
            empty_df = pl.DataFrame(schema=in_schema)
            out_schema = self._process_dataframe(empty_df).schema
            return data.map_batches(self._process_dataframe, schema=out_schema)

        if isinstance(data, str):
            try:
                c, p, r, a = data.split(":")
                return self._annotator.annotate(c, int(p), r, a)
            except ValueError:
                raise ValueError(
                    f"Invalid format '{data}'. Expected 'chr:pos:ref:alt' (e.g., '17:41197701:G:C')"
                )

        if all(v is not None for v in (chromosome, position, reference, alternative)):
            return self._annotator.annotate(
                str(chromosome), int(position), str(reference), str(alternative)
            )

        raise ValueError(
            "Invalid input. Provide a Polars DataFrame/LazyFrame, a variant string ('chr:pos:ref:alt'), "
            "or explicit kwargs (chromosome=..., position=..., reference=..., alternative=...)."
        )

    def annotate_multiple(
        self, variants: typing.Iterable[str | VariantDict]
    ) -> AnnotationResultDict:
        """
        Annotate multiple phased variants together as a single compound event.
        All variants must be on the same chromosome and must not overlap.

        Args:
            variants: An iterable of 'chr:pos:ref:alt' strings OR dictionaries strictly matching
                      the VariantDict type ({'chromosome': str, 'position': int, 'reference': str, 'alternative': str}).
        """
        parsed_variants = []

        for var in variants:
            if isinstance(var, str):
                try:
                    c, p, r, a = var.split(":")
                    parsed_variants.append((c, int(p), r, a))
                except ValueError as e:
                    raise ValueError(
                        f"Invalid format '{var}'. Expected 'chr:pos:ref:alt'"
                    ) from e
            elif isinstance(var, dict):
                try:
                    parsed_variants.append(
                        (
                            str(var["chromosome"]),
                            int(var["position"]),
                            str(var["reference"]),
                            str(var["alternative"]),
                        )
                    )
                except KeyError as e:
                    raise ValueError(
                        f"Variant dict missing required key: {e}. "
                        "Expected keys: 'chromosome', 'position', 'reference', 'alternative'."
                    ) from e
            else:
                raise TypeError("Variants must be strings or VariantDict dictionaries.")

        return self._annotator.annotate_multiple(parsed_variants)


def _resolve_path(p: str | Path | None) -> str | None:
    """Expands '~' and resolves to an absolute path for Rust's PathBuf."""
    if p is None:
        return None
    return str(Path(p).expanduser().resolve())


def build_transcript_db(
    assembly: str,
    annotation: list[str | Path],
    output: str | Path,
    transcript_source: str,
    assembly_version: str | None = None,
    annotation_version: str | None = None,
    transcript_source_version: str | None = None,
    seqrepo: str | Path | None = None,
    transcript_sequences: str | Path | None = None,
    mane_transcripts: str | Path | None = None,
    disable_filters: bool = False,
    threads: int = 1,
    compression_level: int = 19,
) -> None:
    """
    Construct a mehari transcripts and sequence database (.bin.zst) from input annotations.

    Args:
        assembly: Targeted genome assembly (e.g., "GRCh38").
        annotation: List of paths to transcript annotations (cdot JSON or GFF3).
        output: Path to write the compressed database output (.bin.zst).
        transcript_source: Source of the transcripts (e.g., "RefSeq" or "Ensembl").
        assembly_version: Version of the genome assembly (e.g., "GRCh38.p13").
        annotation_version: Version of the annotation data.
        transcript_source_version: Version of the source (e.g., "112" for Ensembl). Required if source is Ensembl.
        seqrepo: Path to a seqrepo instance directory.
        transcript_sequences: Path to FASTA file(s) containing transcript sequences (alternative to seqrepo).
        mane_transcripts: Path to TSV file for MANE label transfer.
        disable_filters: Disable rigorous quality filtering (useful for custom/novel annotations).
        threads: Number of threads to use for parallel processing.
        compression_level: ZSTD compression level (default: 19).

    Raises:
        ValueError: If neither seqrepo nor transcript_sequences are provided, or if Ensembl version is missing.
        RuntimeError: If the database building fails during processing.
    """

    _build_transcript_db(
        assembly=str(assembly),
        annotation=[_resolve_path(p) for p in annotation],
        output=_resolve_path(output),
        transcript_source=str(transcript_source),
        assembly_version=str(assembly_version) if assembly_version else None,
        annotation_version=str(annotation_version) if annotation_version else None,
        transcript_source_version=str(transcript_source_version)
        if transcript_source_version
        else None,
        seqrepo=_resolve_path(seqrepo),
        transcript_sequences=_resolve_path(transcript_sequences),
        mane_transcripts=_resolve_path(mane_transcripts),
        disable_filters=bool(disable_filters),
        threads=int(threads),
        compression_level=int(compression_level),
    )
