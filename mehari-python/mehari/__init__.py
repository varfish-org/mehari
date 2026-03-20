import typing
import polars as pl
import pyarrow as pa
from ._mehari import SeqvarsAnnotator as _SeqvarsAnnotator
from ._mehari import consequence_variants, putative_impact_variants

ConsequenceEnum = pl.Enum(consequence_variants())
ImpactEnum = pl.Enum(putative_impact_variants())


class SeqvarsAnnotator:
    def __init__(self, tx_dbs: str | list[str], ref_path: str | None = None):
        if isinstance(tx_dbs, str):
            tx_dbs = [tx_dbs]
        self._annotator = _SeqvarsAnnotator(tx_dbs, ref_path)

    def _process_dataframe(self, df: pl.DataFrame) -> pl.DataFrame:
        """Internal helper to stream Arrow RecordBatches to Rust."""
        if df.is_empty():
            schema = pa.schema([
                ("chromosome", pa.string()),
                ("position", pa.int32()),
                ("reference", pa.string()),
                ("alternative", pa.string()),
            ])
            dummy_batch = pa.RecordBatch.from_pylist([], schema=schema)
            result_batch = self._annotator.annotate_batch(dummy_batch)
            return pl.from_arrow(pa.Table.from_batches([result_batch]))

        # Normal execution for populated dataframes
        pa_table = df.to_arrow()
        result_batches = [self._annotator.annotate_batch(b) for b in pa_table.to_batches()]
        return pl.from_arrow(pa.Table.from_batches(result_batches))

    @typing.overload
    def annotate(self, data: str) -> dict[str, typing.Any]: ...

    @typing.overload
    def annotate(self, data: pl.DataFrame) -> pl.DataFrame: ...

    @typing.overload
    def annotate(self, data: pl.LazyFrame) -> pl.LazyFrame: ...

    @typing.overload
    def annotate(
        self, *, chromosome: str, position: int, reference: str, alternative: str
    ) -> dict[str, typing.Any]: ...

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
