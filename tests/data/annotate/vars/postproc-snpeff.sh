#!/usr/bin/env bash

INPUT=tests/data/annotate/vars/clinvar.excerpt.snpeff.vcf.gz
OUTPUT=${INPUT%.vcf.gz}.tsv.gz

bcftools query -f "%CHROM-%POS-%REF-%ALT\t%ANN\n" $INPUT \
| awk -F'\t' -v OFS='\t' '
    {
        split($2, a, ",");
        for (i in a) {
            split(a[i], b, "|");
            if (b[4] == "BRCA1") {
                print $1, b[7], b[2];
            }
        }
    }' \
| sort -k1,2 -u \
| gzip -c \
>$OUTPUT
