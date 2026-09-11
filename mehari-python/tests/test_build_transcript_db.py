import logging
import os

from mehari import build_transcript_db

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(TEST_DIR, "..", ".."))
TXS_DIR = os.path.join(PROJECT_ROOT, "mehari", "tests", "data", "db", "create", "txs")


def _build(output):
    build_transcript_db(
        assembly="grch37",
        annotation=[
            os.path.join(TXS_DIR, "cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json")
        ],
        output=output,
        transcript_source="refseq",
        seqrepo=os.path.join(TXS_DIR, "latest"),
    )


def test_build_transcript_db_logs_progress(tmp_path, caplog):
    """
    The build runs on a rayon thread pool. Its progress messages must reach Python
    logging, also if logging was configured only after an earlier build.
    """
    _build(tmp_path / "first.bin.zst")

    with caplog.at_level(logging.INFO, logger="mehari"):
        _build(tmp_path / "second.bin.zst")

    assert (tmp_path / "second.bin.zst").exists()
    assert "Loading annotations" in caplog.text
