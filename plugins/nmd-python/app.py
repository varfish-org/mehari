"""mehari's built-in `nmd` plugin, written in Python against `mehari/wit/plugin.wit`."""

from componentize_py_types import Err
from wit_world import exports
from wit_world.imports.types import Annotation, Prediction, Setup

# Premature stop codons starting within this many coding bases escape NMD by re-initiation.
START_PROXIMAL_BASES = 100
# Premature stop codons starting at most this many bases before the last exon-exon junction
# escape NMD.
LAST_JUNCTION_BASES = 50


class TranscriptPlugin(exports.TranscriptPlugin):
    def init(self, setup: Setup) -> None:
        for key, _ in setup.options():
            raise Err(f"unknown option {key!r}")
        setup.add_string_column("nmd")

    def annotate(self, prediction: Prediction) -> None:
        for annotation in prediction.annotations():
            value = predict(annotation)
            if value is not None:
                annotation.set_string("nmd", value)


def predict(annotation: Annotation) -> str | None:
    """One value for every stop-gained or frameshift annotation."""
    if not {"stop_gained", "frameshift_variant"} & set(annotation.consequences()):
        return None
    return classify(annotation) or "unknown"


def classify(annotation: Annotation) -> str | None:
    """`None` if the premature stop codon (PTC) or the exon structure is unknown."""
    exons = annotation.exons()
    cds_start = annotation.reference_cds_start()
    reference = annotation.reference_transcript()
    alternative = annotation.alternative_transcript()
    protein = annotation.alternative_protein()
    if cds_start is None or reference is None or alternative is None or protein is None:
        return None
    stop = protein.find("*")
    if stop == -1:
        return "no_stop_codon"
    cdna_pos = annotation.cdna_pos()
    if not exons or cdna_pos is None:
        return None

    # 1-based positions on the mutant transcript.
    ptc = cds_start + 3 * stop
    last_exon_start = exons[-1].tx_start
    # Insertions and deletions before the last junction move it.
    last_junction = last_exon_start
    if cdna_pos.ord < last_exon_start:
        last_junction += len(alternative) - len(reference)

    ptc_c = ptc - cds_start + 1
    if last_junction - ptc <= LAST_JUNCTION_BASES or ptc_c <= START_PROXIMAL_BASES:
        return "escaping"
    return "triggering"
