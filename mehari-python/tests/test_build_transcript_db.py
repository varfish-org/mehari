import logging
import os
import sys

from mehari import build_transcript_db

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(TEST_DIR, "..", ".."))
TXS_DIR = os.path.join(PROJECT_ROOT, "mehari", "tests", "data", "db", "create", "txs")


def _build(output, **kwargs):
    build_transcript_db(
        assembly="grch37",
        annotation=[
            os.path.join(TXS_DIR, "cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json")
        ],
        output=output,
        transcript_source="refseq",
        seqrepo=os.path.join(TXS_DIR, "latest"),
        **kwargs,
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


def test_build_transcript_db_reports_progress(tmp_path):
    """
    A tqdm-compatible class gets one bar per step. Each bar reaches its total and
    is closed.
    """
    bars = []

    class RecordingBar:
        def __init__(self, **kwargs):
            self.kwargs = kwargs
            self.n = 0
            self.closed = False
            bars.append(self)

        def update(self, n):
            self.n += n

        def close(self):
            self.closed = True

    _build(tmp_path / "txs.bin.zst", progress=RecordingBar)

    assert [bar.kwargs["desc"] for bar in bars] == [
        "Loading cdot-0.2.22.refseq.grch37_grch38.brca1_opa1.json",
        "Fetching transcript sequences",
        "Writing txs.bin.zst",
    ]
    for bar in bars:
        assert bar.n == bar.kwargs["total"]
        assert bar.closed


def test_build_transcript_db_survives_broken_progress_bar(tmp_path, monkeypatch):
    """
    An exception in the progress bar goes to sys.unraisablehook. The build still
    succeeds.
    """
    unraisable = []
    monkeypatch.setattr(sys, "unraisablehook", unraisable.append)

    class BrokenBar:
        def __init__(self, **kwargs):
            pass

        def update(self, n):
            raise RuntimeError("broken bar")

        def close(self):
            pass

    _build(tmp_path / "txs.bin.zst", progress=BrokenBar)

    assert (tmp_path / "txs.bin.zst").exists()
    assert unraisable
    assert all(isinstance(u.exc_value, RuntimeError) for u in unraisable)
