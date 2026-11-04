import logging
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
    assert nm_007294["hgnc_id"] == "HGNC:1100"
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


def brca1(annotations: list[dict]) -> dict:
    """The annotation of BRCA1 NM_007294.4."""
    return next(c for c in annotations if c.get("feature_id") == "NM_007294.4")


def test_annotate_symbolic_allele(annotator):
    """
    A symbolic allele takes its end from `end`, as from END in a VCF: here the deletion of 50
    bases in the last exon of BRCA1.
    """
    result = annotator.annotate("17:41197710:N:<DEL>", end=41197760)

    nm_007294 = brca1(result["annotation"])
    assert nm_007294["hgvs_g"] == "g.41197711_41197760del"
    assert "frameshift_variant" in nm_007294["consequences"]


def sv_dataframe(variant: str, sv_info: dict) -> pl.DataFrame:
    """The variant 'chr:pos:ref:alt' as a DataFrame, with a column per INFO field."""
    chromosome, position, reference, alternative = variant.split(":", 3)
    columns = {
        "end": None,
        "sv_length": None,
        "sv_claim": None,
        "sv_type": None,
        "cn": None,
        "cipos": None,
        "ciend": None,
    } | sv_info
    return pl.DataFrame(
        {
            "chromosome": [chromosome],
            "position": [int(position)],
            "reference": [reference],
            "alternative": [alternative],
        }
        | {key: [value] for key, value in columns.items()},
        schema_overrides={
            "position": pl.Int32,
            "end": pl.Int32,
            "sv_length": pl.Int32,
            "sv_claim": pl.String,
            "sv_type": pl.String,
            "cn": pl.Float64,
            "cipos": pl.List(pl.Int32),
            "ciend": pl.List(pl.Int32),
        },
    )


def annotate_all_forms(annotator, variant: str, sv_info: dict) -> list[list[dict]]:
    """
    The annotations of the variant from the string form, kwargs, the DataFrame form and the
    LazyFrame form.
    """
    chromosome, position, reference, alternative = variant.split(":", 3)
    df = sv_dataframe(variant, sv_info)
    return [
        annotator.annotate(variant, **sv_info)["annotation"],
        annotator.annotate(
            chromosome=chromosome,
            position=int(position),
            reference=reference,
            alternative=alternative,
            **sv_info,
        )["annotation"],
        annotator.annotate(df).get_column("annotation").to_list()[0],
        annotator.annotate(df.lazy()).collect().get_column("annotation").to_list()[0],
    ]


@pytest.mark.parametrize(
    "variant, sv_info, message",
    [
        ("17:41197710:N:<DEL>", {}, "has no end past its position"),
        ("17:41197710:N:<DEL>", {"end": 41197700}, "has no end past its position"),
        ("17:41197710:N:<FOO>", {"end": 41197760}, "is not supported"),
        ("17:41197710:N:G]13:x]", {}, "invalid mate position"),
        ("17:41197710:N:<DEL>", {"end": 41197760, "sv_claim": "X"}, "SVCLAIM"),
        ("17:41197710:N:<FOO>", {"end": 41197760, "sv_type": "FOO"}, "SVTYPE"),
        ("17:41197710:N:<CNV>", {"end": 41197760, "cn": -1.0}, "CN"),
        ("17:41197710:N:<DEL>", {"end": 41197760, "cipos": [-500]}, "cipos"),
        ("17:41197710:N:<DEL>", {"end": 41197760, "ciend": [500, -500]}, "CIEND"),
    ],
)
def test_annotate_rejects_bad_structural_variant(
    annotator, variant, sv_info, message, caplog
):
    """
    A symbolic allele that cannot be annotated, a malformed breakend, an unknown claim or
    type, a negative copy number and a malformed confidence interval raise a ValueError with
    the reason, from the string form and kwargs. In a DataFrame, the row gets no annotations
    and a warning with the reason, as `annotate seqvars` carries on past such records.
    """
    with pytest.raises(ValueError, match=message):
        annotator.annotate(variant, **sv_info)

    chromosome, position, reference, alternative = variant.split(":", 3)
    with pytest.raises(ValueError, match=message):
        annotator.annotate(
            chromosome=chromosome,
            position=int(position),
            reference=reference,
            alternative=alternative,
            **sv_info,
        )

    with caplog.at_level(logging.WARNING, logger="mehari"):
        result = annotator.annotate(sv_dataframe(variant, sv_info))
    assert result.get_column("annotation").to_list() == [[]]
    assert message in caplog.text


def test_annotate_dataframe_skips_bad_rows(annotator, caplog):
    """A row that cannot be annotated leaves the other rows of its DataFrame annotated."""
    df = pl.DataFrame(
        {
            "chromosome": ["17", "17"],
            "position": [41197701, 41197710],
            "reference": ["G", "N"],
            "alternative": ["C", "<DEL>"],
        },
        schema_overrides={"position": pl.Int32},
    )
    with caplog.at_level(logging.WARNING, logger="mehari"):
        rows = annotator.annotate(df).get_column("annotation").to_list()

    assert "missense_variant" in brca1(rows[0])["consequences"]
    assert rows[1] == []
    assert "has no end past its position" in caplog.text


@pytest.mark.parametrize(
    "variant, sv_info, hgvs_g, imprecise",
    [
        (
            "17:41197710:N:<FOO>",
            {"end": 41197760, "sv_type": "DEL"},
            "g.41197711_41197760del",
            False,
        ),
        (
            "17:41197710:N:<CNV>",
            {"end": 41197760, "cn": 0.4},
            "g.(?_41197711)_(41197760_?)del",
            True,
        ),
        ("17:41197710:N:<CNV>", {"end": 41197760, "cn": 1}, None, False),
        (
            "17:41197710:N:<DEL>",
            {"end": 41197760, "cipos": (-5, 5), "imprecise": True},
            "g.(41197706_41197716)_(41197755_41197765)del",
            True,
        ),
        (
            "17:41197710:N:<DEL>",
            {"end": 41197760, "imprecise": True},
            "g.(?_41197711)_(41197760_?)del",
            True,
        ),
        (
            "17:41197710:N:<DEL>",
            {"end": 41197760, "sv_type": "INVDUP"},
            "g.41197711_41197760del",
            False,
        ),
        (
            "17:41215930:G:G]13:32900000]",
            {"cipos": [-500, 500], "imprecise": True},
            "g.(41215431_41216431)_qterdelins[NC_000013.10:g.pter_32900000inv]",
            True,
        ),
        (
            "17:41215930:G:G]13:32900000]",
            {"sv_type": "TRA"},
            "g.41215931_qterdelins[NC_000013.10:g.pter_32900000inv]",
            False,
        ),
    ],
)
def test_annotate_info_fields(annotator, variant, sv_info, hgvs_g, imprecise):
    """
    A custom allele counts as its sv_type, which other alleles ignore. A <CNV> with 0 copies
    (cn, rounded) is an imprecise deletion and one with 1 copy no variant, and imprecise makes
    the breakpoints imprecise, within cipos where given, also of a breakend; from the string
    form, kwargs, the DataFrame form and the LazyFrame form.
    """
    for annotations in annotate_all_forms(annotator, variant, sv_info):
        if hgvs_g is None:
            assert annotations == []
            continue
        nm_007294 = brca1(annotations)
        assert nm_007294["hgvs_g"] == hgvs_g
        messages = nm_007294["messages"] or []
        assert ("INFO_IMPRECISE_BREAKPOINTS" in messages) == imprecise


def test_annotate_dataframe_rejects_sv_kwargs(annotator, sample_variants):
    """
    A DataFrame holds the INFO fields of symbolic alleles in its columns, so keyword arguments
    for them are an error rather than ignored.
    """
    with pytest.raises(ValueError, match="columns"):
        annotator.annotate(sample_variants, end=41197760)
    with pytest.raises(ValueError, match="columns"):
        annotator.annotate(sample_variants.lazy(), imprecise=True)


def test_annotate_breakend(annotator):
    """
    A breakend needs no more than its ALT. It keeps the bases in front of 17:41215930, so
    BRCA1 on the minus strand loses its 5' end.
    """
    result = annotator.annotate(
        chromosome="17", position=41215930, reference="G", alternative="G]13:32900000]"
    )

    nm_007294 = brca1(result["annotation"])
    assert (
        nm_007294["hgvs_g"] == "g.41215931_qterdelins[NC_000013.10:g.pter_32900000inv]"
    )
    assert "feature_truncation" in nm_007294["consequences"]


@pytest.mark.parametrize("alternative", ["*", "<*>", "<NON_REF>", "."])
def test_annotate_no_variant_of_its_own(annotator, alternative):
    """
    `*`, `<*>`, `<NON_REF>` and ALT `.` are no variants of their own, so they have no
    annotations, as in `annotate seqvars`.
    """
    result = annotator.annotate(f"17:41197701:G:{alternative}")

    assert result["annotation"] == []


def test_annotate_dataframe_with_structural_variants(annotator):
    """
    The optional columns end, sv_length, sv_claim and imprecise hold the INFO fields of
    symbolic alleles; null stands for a missing field.
    """
    df = pl.DataFrame(
        {
            "chromosome": ["17", "17", "17", "17"],
            "position": [41197701, 41197710, 41197710, 41215930],
            "reference": ["G", "N", "N", "G"],
            "alternative": ["C", "<DEL>", "<DEL>", "G]13:32900000]"],
            "end": [None, 41197760, None, None],
            "sv_length": [None, None, -50, None],
            "sv_claim": [None, None, "D", None],
            "imprecise": [None, True, None, None],
        },
        schema_overrides={"position": pl.Int32, "end": pl.Int32, "sv_length": pl.Int32},
    )
    rows = annotator.annotate(df).get_column("annotation").to_list()

    assert "missense_variant" in brca1(rows[0])["consequences"]
    for row in rows[1:3]:
        nm_007294 = brca1(row)
        assert nm_007294["hgvs_g"] == "g.(?_41197711)_(41197760_?)del"
        assert "INFO_IMPRECISE_BREAKPOINTS" in nm_007294["messages"]
    assert "feature_truncation" in brca1(rows[3])["consequences"]


def test_annotate_multiple_with_breakend(annotator):
    """
    A breakend's ALT holds colons. A group with one is not evaluated together, as in
    `annotate seqvars`, so it gives no annotations.
    """
    result = annotator.annotate_multiple(
        ["17:41215930:G:G]13:32900000]", "17:41197701:G:C"]
    )

    assert result["annotation"] == []


def test_annotate_batch_forwards_rust_warnings(annotator, sample_variants, caplog):
    """
    Batch annotation logs from rayon worker threads. These log records must reach
    Python logging without deadlocking on the GIL.
    """
    unknown_chromosome = sample_variants.with_columns(chromosome=pl.lit("FOO"))

    with caplog.at_level(logging.WARNING, logger="mehari"):
        result_df = annotator.annotate(unknown_chromosome)

    assert result_df.get_column("annotation").to_list() == [[], []]
    assert "Could not determine chromosome accession" in caplog.text


def test_annotate_reports_protein_sequences(sample_variants):
    """
    Sequence reporting fills the custom_fields struct, in column order.
    """
    annotator = SeqvarsAnnotator(
        transcript_db_paths=os.path.join(DB_DIR, "txs.bin.zst"),
        plugins=["aa-sequence-ref", "aa-sequence-alt"],
    )
    result_df = annotator.annotate(sample_variants)

    nm_007294 = next(
        c
        for c in result_df.row(0, named=True)["annotation"]
        if c.get("feature_id") == "NM_007294.4"
    )
    custom_fields = nm_007294["custom_fields"]
    assert list(custom_fields) == ["aa_sequence_ref", "aa_sequence_alt"]
    assert custom_fields["aa_sequence_ref"] != custom_fields["aa_sequence_alt"]


def test_annotate_with_nmd_plugin():
    """
    The NMD plugin predicts NMD for a stop-gained variant in BRCA1 exon 10 (c.3020).
    """
    annotator = SeqvarsAnnotator(
        transcript_db_paths=os.path.join(DB_DIR, "txs.bin.zst"),
        plugins=["nmd"],
    )
    result = annotator.annotate("17:41244528:G:T")

    nm_007294 = next(
        c for c in result["annotation"] if c.get("feature_id") == "NM_007294.4"
    )
    assert "stop_gained" in nm_007294["consequences"]
    assert nm_007294["custom_fields"] == {"nmd": "triggering"}


def test_unknown_plugin_is_rejected():
    with pytest.raises(ValueError, match="Invalid plugin"):
        SeqvarsAnnotator(
            transcript_db_paths=os.path.join(DB_DIR, "txs.bin.zst"),
            plugins=["no-such-plugin"],
        )
