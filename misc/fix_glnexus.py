#!/usr/bin/env python3
"""
Fix GLNexus output for input to `mehari annotate seqvars`.

First, note that multiallelic sites must be corrected by a call to
`bcftools norm` already as an input to this script.

What we will do is splitting the `FORMAT/RNC` field with a comma
as noodles-vcf expects a list of characters for rather than a
string of length 2 for `##FORMAT=<ID=RNC,Number=2,Type=Character>`.
"""

from pathlib import Path
import sys
from typing import Annotated

import typer
import vcfpy


def main(
    path_in: Annotated[Path, typer.Option(help="Path to input VCF")],
    path_out: Annotated[Path, typer.Option(help="Path to output VCF")],
    quiet: Annotated[bool, typer.Option(help="Disable verbose output")] = False,
):
    """
    Fix GLNexus output to be compatible with `noodles-vcf`.

    cf: https://github.com/varfish-org/mehari/issues/356
    """
    if not quiet:
        print(f"Opening input file {path_in}", file=sys.stderr)
    reader = vcfpy.Reader.from_path(path_in)
    header = reader.header.copy()
    for line in header.lines:
        if line.key == "FORMAT" and line.mapping.get("ID") == "GQ":
            line.mapping["Number"] = 1
            line.mapping["Type"] = "Integer"
    if not quiet:
        print(f"Opening output file {path_out}", file=sys.stderr)
    writer = vcfpy.Writer.from_path(path_out, header)
    if not quiet:
        print("Processing records...", file=sys.stderr)
    with reader, writer:
        for idx, record in enumerate(reader):
            if idx % 10_000 == 0:
                print(
                    f"  at {idx} records {record.CHROM}:{record.POS}", file=sys.stderr
                )
            for call in record.calls:
                if "RNC" in call.data:
                    if (
                        isinstance(call.data["RNC"], list)
                        and len(call.data["RNC"]) == 1
                    ):
                        call.data["RNC"] = list(call.data["RNC"][0])
            writer.write_record(record)
    if not quiet:
        print("... done", file=sys.stderr)
        print("All done. Have a nice day!", file=sys.stderr)


if __name__ == "__main__":
    typer.run(main)
