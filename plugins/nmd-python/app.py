"""mehari's built-in `nmd` plugin, written in Python against `mehari/wit/plugin.wit`."""

import wit_world
from componentize_py_types import Err
from wit_world.imports.types import Annotations, VcfVariant

# Premature stop codons starting within this many coding bases escape NMD by re-initiation.
START_PROXIMAL_BASES = 100
# Premature stop codons starting at most this many bases before the last exon-exon junction
# escape NMD.
LAST_JUNCTION_BASES = 50


class WitWorld(wit_world.WitWorld):
    def init(self, options: list[tuple[str, str]]) -> list[str]:
        if options:
            raise Err(f"unknown option {options[0][0]!r}")
        return ["nmd"]

    def annotate(
        self, variants: list[VcfVariant], annotations: Annotations
    ) -> list[list[str | None]]:
        return [[predict(annotations, index)] for index in range(annotations.len())]


def predict(annotations: Annotations, index: int) -> str | None:
    """One value for every stop-gained or frameshift annotation."""
    if not {"stop_gained", "frameshift_variant"} & set(annotations.consequences(index)):
        return None
    return classify(annotations, index) or "unknown"


def classify(annotations: Annotations, index: int) -> str | None:
    """`None` if the premature stop codon (PTC) or the exon structure is unknown."""
    transcript = annotations.transcript(index)
    reference = annotations.reference(index)
    alternative = annotations.alternative(index)
    if transcript is None or reference is None or alternative is None:
        return None
    stop = alternative.protein.find("*")
    if stop == -1:
        return "no_stop_codon"
    cdna_pos = annotations.cdna_pos(index)
    if not transcript.exons or cdna_pos is None:
        return None

    # 1-based positions on the mutant transcript.
    ptc = reference.cds_start + 3 * stop
    last_exon_start = transcript.exons[-1].tx_start
    # Insertions and deletions before the last junction move it.
    last_junction = last_exon_start
    if cdna_pos.ord < last_exon_start:
        last_junction += len(alternative.transcript) - len(reference.transcript)

    ptc_c = ptc - reference.cds_start + 1
    if last_junction - ptc <= LAST_JUNCTION_BASES or ptc_c <= START_PROXIMAL_BASES:
        return "escaping"
    return "triggering"
