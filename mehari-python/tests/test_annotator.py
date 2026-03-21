import os

import polars as pl
import pytest

from mehari import SeqvarsAnnotator

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(TEST_DIR, "..", ".."))
DB_DIR = os.path.join(
    PROJECT_ROOT, "mehari", "tests", "data", "annotate", "db", "grch37"
)


@pytest.fixture(scope="module")
def annotator():
    """
    Initializes the Mehari annotator once for the entire test module.
    """
    tx_db_path = os.path.join(DB_DIR, "txs.bin.zst")
    assert os.path.exists(tx_db_path), f"Transcript DB not found at {tx_db_path}"

    return SeqvarsAnnotator(transcript_db_paths=tx_db_path, reference_path=None)


@pytest.fixture
def sample_variants() -> pl.DataFrame:
    """Provides a standard eager Polars DataFrame for batch testing."""
    return pl.DataFrame(
        {
            "chromosome": ["17", "3"],
            "position": [41197701, 193332511],
            "reference": ["G", "G"],
            "alternative": ["C", "T"],
        },
        schema={
            "chromosome": pl.String,
            "position": pl.Int32,
            "reference": pl.String,
            "alternative": pl.String,
        },
    )


def test_annotate_string_format(annotator):
    """
    Tests the single variant string routing ('chr:pos:ref:alt').
    Tests coordinate: 17:41197701:G:C (exonic BRCA1)
    """
    result = annotator.annotate("17:41197701:G:C")

    assert "annotation" in result
    annotations = result["annotation"]
    assert len(annotations) > 0

    nm_007294 = next(
        (c for c in annotations if c.get("feature_id") == "NM_007294.4"), None
    )
    assert nm_007294 is not None
    assert nm_007294["distance"] == 0
    assert nm_007294["strand"] == -1
    assert nm_007294["putative_impact"] == "MODERATE"
    assert "missense_variant" in nm_007294["consequences"]


def test_annotate_kwargs(annotator):
    """
    Tests the explicit kwargs routing.
    Tests coordinate: 3:193332511:G:T (-1bp intronic OPA1)
    """
    result = annotator.annotate(
        chromosome="3", position=193332511, reference="G", alternative="T"
    )

    annotations = result["annotation"]
    assert len(annotations) > 0

    nm_130837 = next(
        (c for c in annotations if c.get("feature_id") == "NM_130837.3"), None
    )
    assert nm_130837 is not None
    assert nm_130837["distance"] == -1
    assert "splice_acceptor_variant" in nm_130837["consequences"]
    assert "coding_transcript_intron_variant" in nm_130837["consequences"]


def test_annotate_eager_dataframe(annotator, sample_variants):
    """
    Tests the batch processing routing using an eager Polars DataFrame.
    """
    result_df = annotator.annotate(sample_variants)

    assert isinstance(result_df, pl.DataFrame)
    assert "annotation" in result_df.columns
    assert len(result_df) == 2

    brca1_row = result_df.row(0, named=True)
    brca1_annotations = brca1_row["annotation"]

    nm_007294 = next(
        (c for c in brca1_annotations if c.get("feature_id") == "NM_007294.4"), None
    )
    assert nm_007294 is not None
    assert nm_007294["putative_impact"] == "MODERATE"
    assert "missense_variant" in nm_007294["consequences"]


def test_annotate_lazy_streaming(annotator, sample_variants):
    """
    Tests the batch processing routing using a LazyFrame to ensure schema mapping
    is configured correctly for out-of-core streaming.
    """
    lazy_df = sample_variants.lazy()

    # This should not compute anything yet, just build the query plan
    result_lazy = annotator.annotate(lazy_df)

    # Verify it is still lazy
    assert isinstance(result_lazy, pl.LazyFrame)

    # Trigger computation
    result_df = result_lazy.collect()

    # Verify OPA1 logic (Row 1)
    opa1_row = result_df.row(1, named=True)
    opa1_annotations = opa1_row["annotation"]

    nm_130837 = next(
        (c for c in opa1_annotations if c.get("feature_id") == "NM_130837.3"), None
    )
    assert nm_130837 is not None
    assert "splice_acceptor_variant" in nm_130837["consequences"]
