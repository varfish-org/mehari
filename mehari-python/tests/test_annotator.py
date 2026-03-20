import os

import mehari
import polars as pl
import pytest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(TEST_DIR, "..", ".."))
DB_DIR = os.path.join(PROJECT_ROOT, "mehari", "tests", "data", "annotate", "db", "grch37")


@pytest.fixture(scope="module")
def annotator():
    """
    Initializes the Mehari annotator once for the entire test module
    to avoid reloading the .zst databases for every test.
    """
    tx_db_path = os.path.join(DB_DIR, "txs.bin.zst")

    assert os.path.exists(tx_db_path), f"Transcript DB not found at {tx_db_path}"

    return mehari.SeqvarsAnnotator(
        transcript_db_paths=[tx_db_path], reference_path=None
    )


def test_annotate_brca1_exonic(annotator):
    """
    Port of `annotate_snv_brca1_one_variant` from csq.rs.
    Tests coordinate: 17:41197701:G:C (exonic)
    """
    result = annotator.annotate(
        chrom="17", position=41197701, reference="G", alternative="C"
    )

    assert "consequences" in result
    csqs = result["consequences"]

    assert len(csqs) > 0

    nm_007294 = next((c for c in csqs if c.get("feature_id") == "NM_007294.4"), None)
    assert nm_007294 is not None

    assert nm_007294["distance"] == 0

    assert nm_007294["strand"] == -1  # BRCA1 is on the minus strand
    assert nm_007294["putative_impact"] == "moderate"
    assert "missense_variant" in nm_007294["consequences"]


def test_annotate_opa1_intronic(annotator):
    """
    Port of `annotate_snv_opa1_csq` from csq.rs.
    Tests coordinate: 3:193332511:G:T (-1bp intronic)
    """
    result = annotator.annotate(
        chrom="3", position=193332511, reference="G", alternative="T"
    )

    csqs = result["consequences"]
    assert len(csqs) > 0

    nm_130837 = next((c for c in csqs if c.get("feature_id") == "NM_130837.3"), None)
    assert nm_130837 is not None

    assert nm_130837["distance"] == -1

    assert "splice_acceptor_variant" in nm_130837["consequences"]
    assert "coding_transcript_intron_variant" in nm_130837["consequences"]


def test_annotate_batch_polars(annotator):
    df = pl.DataFrame(
        {
            "chromosome": ["17", "3"],
            "position": [41197701, 193332511],
            "reference": ["G", "G"],
            "alternative": ["C", "T"],
        }
    ).with_columns(pl.col("position").cast(pl.Int32))

    arrow_batch = df.to_arrow().to_batches()[0]
    result_batch = annotator.annotate_batch(arrow_batch)
    result_df = pl.from_arrow(result_batch)

    assert "consequences" in result_df.columns
    assert len(result_df) == 2

    brca1_csqs = result_df["consequences"][0]
    assert len(brca1_csqs) > 0
    nm_007294 = next(
        (c for c in brca1_csqs if c.get("feature_id") == "NM_007294.4"), None
    )
    assert nm_007294 is not None
    assert nm_007294["putative_impact"] == "moderate"
    assert "missense_variant" in nm_007294["consequences"]

    opa1_csqs = result_df["consequences"][1]
    assert len(opa1_csqs) > 0
    nm_130837 = next(
        (c for c in opa1_csqs if c.get("feature_id") == "NM_130837.3"), None
    )
    assert nm_130837 is not None
    assert "splice_acceptor_variant" in nm_130837["consequences"]
