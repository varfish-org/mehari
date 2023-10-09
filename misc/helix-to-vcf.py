#!/usr/bin/env python
"""Helper script to convert the HelixMtDb format to VCF.

The resulting file will look like this::

    ##fileformat=VCFv4.2
    ##contig=<ID=chrM,length=16569>
    ##FILTER=<ID=PASS,Description="Variant passes all filters">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Overall allele number (Number of alleles with non-missing genotype)">
    ##INFO=<ID=AC_hom,Number=1,Type=Integer,Description="Allele counts called as homoplasmic">
    ##INFO=<ID=AC_het,Number=1,Type=Integer,Description="Alelle counts called as heteroplasmic">
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
    chrM    5       .       A       C       .       PASS    AN=196554;AC_hom=1;AC_het=0
    chrM    10      .       T       C       .       PASS    AN=196554;AC_hom=7;AC_het=1
    chrM    11      .       C       T       .       PASS    AN=196554;AC_hom=0;AC_het=1
"""

import csv
import json
import sys

import vcfpy

#: Number of individiduals in HelixMtDb
AN = 196_554


def build_writer():
    """Build ``vcfpy`` writer for writing out the VCF file."""
    header = vcfpy.Header(samples=vcfpy.SamplesInfos(sample_names=[]))
    header.add_line(vcfpy.HeaderLine("fileformat", "VCFv4.2"))
    header.add_contig_line({"ID": "chrM", "length": 16569})
    header.add_filter_line({"ID": "PASS", "Description": "Variant passes all filters"})
    header.add_info_line(
        {
            "ID": "AN",
            "Number": 1,
            "Type": "Integer",
            "Description": "Overall allele number (Number of alleles with non-missing genotype)",
        }
    )
    header.add_info_line(
        {
            "ID": "AC_hom",
            "Number": 1,
            "Type": "Integer",
            "Description": "Allele counts called as homoplasmic",
        }
    )
    header.add_info_line(
        {
            "ID": "AC_het",
            "Number": 1,
            "Type": "Integer",
            "Description": "Alelle counts called as heteroplasmic",
        }
    )
    return vcfpy.Writer(sys.stdout, header)


def main():
    writer = build_writer()
    reader = csv.DictReader(sys.stdin, delimiter="\t")
    for row in reader:
        locus = row["locus"].split(":")
        chrom, pos = locus[0], int(locus[1])
        alleles = json.loads(row["alleles"])
        ac_hom = int(row["counts_hom"])
        ac_het = int(row["counts_het"])
        writer.write_record(
            vcfpy.Record(
                CHROM=chrom,
                POS=pos,
                ID=[],
                REF=alleles[0],
                ALT=[vcfpy.Substitution(vcfpy.SNV, alleles[1])],
                QUAL=None,
                FILTER=["PASS"],
                INFO={
                    "AN": AN,
                    "AC_hom": ac_hom,
                    "AC_het": ac_het,
                },
                FORMAT=None,
                calls=None,
            )
        )


if __name__ == "__main__":
    sys.exit(main())
