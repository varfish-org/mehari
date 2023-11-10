#!/usr/bin/env python
"""Helper script to convert ClinCNV files to VCF format for import."""

import csv
import enum
from typing import Annotated
from pathlib import Path
import logging

from bioutils import assemblies
import pydantic
import logzero
from logzero import logger
import typer
import vcfpy


@enum.unique
class GenomeRelease(enum.Enum):
    """Enumeration for genome releases"""

    #: GRCh37 release
    GRCH37 = "GRCh37"
    #: GRCh38 release
    GRCH38 = "GRCh38"


#: Canonical contig names.
CANONICAL_CONTIGS = [
    *[str(i) for i in range(1, 23)],
    "X", "Y", "MT"
]

def get_clincnv_version(path: Path) -> str:
    """Get ClinCNV version from file."""
    with path.open("rt") as inputf:
        for line in inputf:
            if line.startswith("##ClinCNV version:"):
                return line.strip().split(": ")[1].strip()
            elif not line.startswith("#"):
                break
    raise RuntimeError("Could not determine ClinCNV version")


def create_header(
    sample_name: str, genome_release: GenomeRelease, clincnv_version: str
) -> vcfpy.Header:
    """Create header with the given values."""
    header = vcfpy.Header(samples=vcfpy.SamplesInfos([sample_name]))
    header.add_line(vcfpy.HeaderLine(key="fileformat", value="VCFv4.2"))
    header.add_line(vcfpy.HeaderLine(key="source", value=f"ClinCNV {clincnv_version}"))
    # FILTER
    header.add_filter_line({"ID": "PASS", "Description": "All filters passed"})
    header.add_filter_line(
        {
            "ID": "LowQual",
            "Description": "Loglikelyhood (reported as GQ) is less than 20q",
        }
    )
    # INFO
    header.add_info_line(
        {
            "ID": "SVLEN",
            "Number": ".",
            "Type": "Integer",
            "Description": "Difference in length between REF and ALT alleles",
        }
    )
    header.add_info_line(
        {
            "ID": "END",
            "Number": 1,
            "Type": "Integer",
            "Description": "End position of the variant described in this record",
        }
    )
    header.add_info_line(
        {
            "ID": "SVTYPE",
            "Number": 1,
            "Type": "String",
            "Description": "Type of structural variant",
        }
    )
    # FORMAT
    header.add_format_line(
        {
            "ID": "CN",
            "Number": 1,
            "Type": "Integer",
            "Description": "Segment most-likely copy-number call",
        }
    )
    header.add_format_line(
        {
            "ID": "GT",
            "Number": 1,
            "Type": "String",
            "Description": "Segment genotype 0 or 1",
        }
    )
    header.add_format_line(
        {"ID": "NP", "Number": 1, "Type": "Integer", "Description": "Number of regions"}
    )
    header.add_format_line(
        {
            "ID": "GQ",
            "Number": 1,
            "Type": "Integer",
            "Description": "Loglikelyhood of call",
        }
    )
    # CONTIG
    if genome_release == GenomeRelease.GRCH37:
        # we need p13 as we want chrMT
        assembly_info = assemblies.get_assembly("GRCh37.p12")
    else:
        assembly_info = assemblies.get_assembly("GRCh38")
    for seq in assembly_info["sequences"]:
        if seq["name"].replace("chr", "") in CANONICAL_CONTIGS:
            header.add_contig_line(
                {
                    "ID": seq["name"] if genome_release == GenomeRelease.GRCH37 else f"chr{seq['name']}",
                    "length": seq["length"],
                    "assembly": genome_release.value,
                }
            )

    return header


class ClinCnvRecord(pydantic.BaseModel):
    chr: str
    start: int
    end: int
    cn_change: int = pydantic.Field(validation_alias="CN_change")
    loglikelihood: int
    no_of_regions: int
    length_kb: str = pydantic.Field(validation_alias="length_KB")
    potential_af: str = pydantic.Field(validation_alias="potential_AF")
    genes: str
    qvalue: str
    overlap_af_genomes_imgag: str = pydantic.Field(validation_alias="overlap af_genomes_imgag")
    cn_pathogenic: str
    dosage_sensitive_disease_genes: str
    clinvar_cnvs: str
    omim: str
    gene_info: str
    ngsd_pathogenic_cnvs: str



def run_processing(
    path_in: Path,
    path_out: Path,
    header: vcfpy.Header,
):
    """Run the actual processing"""
    logger.info("    - skipping header")
    with path_in.open("rt") as inputf:
        # skip over ## header
        while True:
            pos = inputf.tell()
            line = inputf.readline()
            if line is None:
                break
            if not line.startswith("##"):
                break
        inputf.seek(pos)
        inputf.read(1)  # skip comment from #chr... header

        reader = csv.DictReader(inputf, delimiter="\t")
        with vcfpy.Writer.from_path(path_out, header) as writer:
            for raw_record in reader:
                clincnv_record = ClinCnvRecord(**raw_record)
                record = vcfpy.Record(
                    CHROM=clincnv_record.chr,
                    POS=clincnv_record.start,
                    ID=[],
                    QUAL=None,
                    REF="N",
                    ALT=[vcfpy.SymbolicAllele("DEL" if clincnv_record.cn_change < 2 else "DUP")],
                    FILTER=[] if clincnv_record.loglikelihood >= 20 else ["LowQual"],
                    INFO={
                        "END": clincnv_record.end,
                        "SVTYPE": "DEL" if clincnv_record.cn_change <2 else "DUP",
                        "SVLEN": [clincnv_record.end - clincnv_record.start + 1],
                    },
                    FORMAT=["GT", "CN", "GQ", "NP"],
                    calls=[vcfpy.Call(sample=header.samples.names[0], data={
                        "GT": "1",
                        "CN": clincnv_record.cn_change,
                        "GQ": clincnv_record.loglikelihood,
                        "NP": clincnv_record.no_of_regions,
                    })]
                )

                writer.write_record(record)


def main(
    path_in: Annotated[Path, typer.Option(help="Path to input TSV file")],
    path_out: Annotated[Path, typer.Option(help="Path to output VCF")],
    sample_name: Annotated[str, typer.Option(help="Sample name to use in VCF header")],
    genome_release: Annotated[
        GenomeRelease,
        typer.Option(
            help="Genome release to use for VCF header",
            show_default=True,
        ),
    ],
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable verbose output")
    ] = False,
):
    """Convert ClinCNV TSV file to VCF format.

    Notes:

    Loglikelyhood is written out as in GQ tag.
    """
    # Setup logging
    if verbose:  # pragma: no cover
        level = logging.DEBUG
    else:
        # Remove module name and line number if not running in debug mode.s
        formatter = logzero.LogFormatter(
            fmt="%(color)s[%(levelname)1.1s %(asctime)s]%(end_color)s %(message)s"
        )
        logzero.formatter(formatter)
        level = logging.INFO
    logzero.loglevel(level=level)

    logger.info("Starting conversion...")
    logger.info("  - getting ClinCNV version")
    clincnv_version = get_clincnv_version(path_in)
    logger.info("  - creating header")
    header = create_header(sample_name, genome_release, clincnv_version)
    logger.info("  - processing records")
    run_processing(path_in, path_out, header)
    logger.info("... done with conversion")

    logger.info("All done. Have a nice day!")


if __name__ == "__main__":
    typer.run(main)
