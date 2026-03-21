# mehari

Python bindings for the `mehari` Rust library.

## Features

* **Single variants:** Annotate a single variant using a format string (`chr:pos:ref:alt`) or keyword arguments.
* **DataFrames:** Process batches of variants by passing a `polars.DataFrame`.
* **LazyFrames:** Support for `polars.LazyFrame` to process large datasets (like Parquet files) without loading everything into memory.

## Usage

Initialize `SeqvarsAnnotator` with your transcript database (see [`mehari-data-tx`](https://github.com/varfish-org/mehari-data-tx/releases)) and a reference genome (FASTA, uncompressed, with index).

```python
from mehari import SeqvarsAnnotator

annotator = SeqvarsAnnotator(
    transcript_db_paths=["path/to/txs.bin.zst"],
    reference_path="path/to/reference.fa"
)
```

To annotate a single variant either use colon separated format string or keyword arguments:
```python
result1 = annotator.annotate("17:41197701:G:C")
result2 = annotator.annotate(chromosome="3", position=193332511, reference="G", alternative="T")
```

To annotate a batch of variants, pass a `polars.DataFrame` or `polars.LazyFrame`.
```python
import polars as pl

df = pl.DataFrame(
    {
        "chromosome": ["17", "3"],
        "position": [41197701, 193332511],
        "reference": ["G", "G"],
        "alternative": ["C", "T"],
    },
    schema={
        "chromosome": pl.Categorical, "position": pl.Int32,
        "reference": pl.String, "alternative": pl.String
    }
)

annotated_df = annotator.annotate(df)
```

## Schemas and types

### Enums
Mehari exports its internal enums to Python so you can use them for filtering or comparisons:
```python
from mehari import ConsequenceEnum, ImpactEnum
```

### DataFrame Schema

When annotating a DataFrame or LazyFrame, mehari appends an "annotation" column.
This column is a polars `List(Struct)` with the following fields:
- `allele`: `String`
- `consequences`: `List(ConsequenceEnum)`
- `putative_impact`: `ImpactEnum`
- `gene_symbol`: `String`
- `gene_id`: `String`
- `feature_type`: `String`
- `feature_id`: `String`
- `feature_biotype`: `List(String)`
- `feature_tags`: `List(String)`
- `rank`: `Struct(ord: Int32, total: Int32)`
- `cdna_pos`: `Struct(ord: Int32, total: Int32)`
- `cds_pos`: `Struct(ord: Int32, total: Int32)`
- `protein_pos`: `Struct(ord: Int32, total: Int32)`
- `hgvs_g`: `String`
- `hgvs_n`: `String`
- `hgvs_c`: `String`
- `hgvs_p`: `String`
- `distance`: `Int32`
- `strand`: `Int32`
- `messages`: `List(String)`
