#!/usr/bin/env python
"""Helper script to fix FreeBayes VCF file."""

from typing import Annotated
from pathlib import Path
import sys

import typer
import vcfpy


def main(
    path_in: Annotated[Path, typer.Option(help="Path to input VCF")],
    path_out: Annotated[Path, typer.Option(help="Path to output VCF")],
    quiet: Annotated[bool, typer.Option(help="Disable verbose output")] = False,
):
    """
    Fix FreeBayes VCF files to be compatible with the VCF4.2 standard.

    - Ensure the FORMAT=GQ field is an Integer.
    - If AD is missing, derive AD from DP and AO.
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
            if "AD" not in record.FORMAT and "AO" in record.FORMAT and "DP" in record.FORMAT:
                record.FORMAT.append("AD")
            for call in record.calls:
                if "GQ" in call.data:
                    call.data["GQ"] = int(call.data["GQ"])
                if "AD" not in call.data and "AO" in call.data and "DP" in call.data:
                    assert len(call.data["AO"]) == 1
                    call.data["AD"] = [call.data["DP"] - call.data["AO"][0], call.data["AO"][0]]
            writer.write_record(record)
    if not quiet:
        print("... done", file=sys.stderr)
        print("All done. Have a nice day!", file=sys.stderr)


if __name__ == "__main__":
    typer.run(main)
