"""mehari's built-in `nmd` plugin, written in Python against `mehari/wit/plugin.wit`."""

import wit_world
from wit_world.imports.types import Annotation, Contexts, VcfVariant

# Premature stop codons starting within this many coding bases escape NMD by re-initiation.
START_PROXIMAL_BASES = 100
# Premature stop codons starting at most this many bases before the last exon-exon junction
# escape NMD.
LAST_JUNCTION_BASES = 50


class WitWorld(wit_world.WitWorld):
    def columns(self) -> list[str]:
        return ["nmd"]

    def annotate(
        self,
        variants: list[VcfVariant],
        annotations: list[Annotation],
        contexts: Contexts,
    ) -> list[list[str | None]]:
        return [[predict(contexts, index, ann)] for index, ann in enumerate(annotations)]


def predict(contexts: Contexts, index: int, ann: Annotation) -> str | None:
    """One value for every stop-gained or frameshift annotation."""
    if not {"stop_gained", "frameshift_variant"} & set(ann.consequences):
        return None
    return classify(contexts, index, ann) or "unknown"


def classify(contexts: Contexts, index: int, ann: Annotation) -> str | None:
    """`None` if the premature stop codon (PTC) or the exon structure is unknown."""
    transcript = contexts.transcript(index)
    reference = contexts.reference(index)
    alternative = contexts.alternative(index)
    if transcript is None or reference is None or alternative is None:
        return None
    stop = alternative.protein.find("*")
    if stop == -1:
        return "no_stop_codon"
    if not transcript.exons or ann.cdna_pos is None:
        return None

    # 1-based positions on the mutant transcript.
    ptc = reference.cds_start + 3 * stop
    last_exon_start = transcript.exons[-1].tx_start
    # Insertions and deletions before the last junction move it.
    last_junction = last_exon_start
    if ann.cdna_pos < last_exon_start:
        last_junction += len(alternative.transcript) - len(reference.transcript)

    ptc_c = ptc - reference.cds_start + 1
    if last_junction - ptc <= LAST_JUNCTION_BASES or ptc_c <= START_PROXIMAL_BASES:
        return "escaping"
    return "triggering"
