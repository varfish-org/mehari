#!/usr/bin/env bash

# Post-process VEP results from the clinvar excerpt.
#
# The file will contain the variant in column 1, transcript in column 2, and
# the effect in column 3.  We will limit to to the transcripts in the JSON
# file.

JSON=tests/data/db/create/txs/cdot-0.2.12.refseq.grch37_grch38.brca1_opa1.json
INPUT=tests/data/annotate/vars/clinvar.excerpt.vep.vcf.gz
OUTPUT=${INPUT%.vcf.gz}.tsv

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" ERR EXIT

grep '"id": "NM_' $JSON \
| cut -b 14-24 \
> $TMPDIR/txs.txt

bcftools query -f "%CHROM-%POS-%REF-%ALT\t%CSQ\n" $INPUT \
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
| grep -f $TMPDIR/txs.txt \
| sort -k1,2 -u \
>$OUTPUT
