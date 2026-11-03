# mehari

Python bindings for the [`mehari`](https://github.com/varfish-org/mehari) Rust library.

## Features

* **Single variants:** Annotate a single variant using a format string (`chr:pos:ref:alt`) or keyword arguments.
* **Multiple variants (_experimental_):** Evaluate the compound effect of multiple variants.
* **DataFrames:** Process batches of variants by passing a `polars.DataFrame`.
* **LazyFrames:** Support for `polars.LazyFrame` to process large datasets (like Parquet files) without loading
  everything into memory.

## Usage

Initialize `SeqvarsAnnotator` with your transcript database (see [
`mehari-data-tx`](https://github.com/varfish-org/mehari-data-tx/releases)) and a reference genome (FASTA, uncompressed,
with index).

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

A symbolic allele such as `<DEL>` takes the INFO fields of its VCF record as keyword arguments: `end` (`END`), `sv_length` (`SVLEN`), `sv_type` (`SVTYPE`, which a custom allele such as `<FOO>` counts as), `sv_claim` (`SVCLAIM`), `cn` (`CN`), `cipos` (`CIPOS`), `ciend` (`CIEND`) and `imprecise` (`IMPRECISE`).
`cn` is the copy number of the allele as in VCF 4.4, rounded; for an older file, convert its total copy number first. A `<CNV>` with `cn=1` is no variant and gets no annotations.
`cipos` and `ciend` are pairs `(min, max)`, as a tuple or a list.
A breakend needs only its ALT, and takes `imprecise` and `cipos` for the uncertainty of its position.
A variant that cannot be annotated, e.g. a `<DEL>` without `end`, raises a `ValueError` with the reason.
The result lists the gene fusions of a structural variant or breakend in `fusions`, e.g. `["OPA1::BRCA1"]`, as the `FUSION` INFO field of `mehari annotate seqvars` does.

```python
result3 = annotator.annotate(
    "17:41197710:N:<DEL>", end=41197760, imprecise=True, cipos=(-5, 5)
)
result4 = annotator.annotate(
    "17:41215930:G:G]13:32900000]", imprecise=True, cipos=(-500, 500)
)
```

To annotate multiple phased variants together as a single compound event (Experimental):
> **Note:** Mehari does not infer phasing.
> When using `annotate_multiple`, mehari assumes all provided variants are on the same chromosome, exist on the same
> haplotype, and do not overlap.

```python
result1 = annotator.annotate_multiple(["1:37799635:TA:A", "1:37799639:C:CG"])

result2 = annotator.annotate_multiple([
    {"chromosome": "1", "position": 37799635, "reference": "TA", "alternative": "A"},
    {"chromosome": "1", "position": 37799639, "reference": "C", "alternative": "CG"}
])
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

For symbolic alleles and breakends, the DataFrame may hold the columns `end`, `sv_length`, `sv_type`, `sv_claim`, `cn`, `cipos`, `ciend` and `imprecise`; a null value stands for a missing field.
`cipos` and `ciend` are list columns of two integers each, `pl.List` or `pl.Array`, e.g. `[-500, 500]`.
As keyword arguments next to a DataFrame, they raise a `ValueError`.
A row that cannot be annotated gets no annotations and a warning with the reason, and the other rows are annotated.

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
- `hgnc_id`: `String`
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

It also appends a "fusions" column, a `List(String)` of the gene fusions of each structural variant or breakend.


## Building a transcript database
To build a transcript database, you can use the `build_transcript_db` function:

```python
from mehari import build_transcript_db
build_transcript_db(
    assembly="grch38",
    annotation=["grch38.gff.gz"],
    transcript_sequences="grch38.fasta",
    transcript_source="ensembl",
    output="grch38.bin.zst"
)
```

This may take a while (several minutes for GRCh38 + Ensembl).
Most of the time goes into compressing the database at the default `compression_level=19`.
To show a progress bar for each step, pass a tqdm-compatible class as `progress`, for example `progress=tqdm.auto.tqdm`.
To follow the steps in between, enable [logging](#logging).

## Logging

mehari writes its log messages to Python's `logging` module, under the logger `mehari`.
By default, Python shows only warnings and errors.
To also see progress messages, enable the `INFO` level:

```python
import logging

logging.basicConfig(level=logging.INFO)
```

Configure logging before you create a `SeqvarsAnnotator` or call `build_transcript_db`.
mehari caches the log levels and rereads them only when one of these two starts.