#!/usr/bin/bash

# Setup Logging -------------------------------------------------------------

log()
{
    >&2 echo $@
}

debug()
{
    [[ "${VERBOSE-0}" -ne 0 ]] && >&2 echo $@
}

set -euo pipefail

if [[ "${VERBOSE-0}" -ne 0 ]]; then
    set -x
fi

# Initialization ------------------------------------------------------------

if [[ "$#" -ne 0 ]]; then
    log "USAGE: bootstrap.sh"
    log ""
    log "Set VERBOSE=1 to increase verbosity."
    exit 1
fi

# path to the directory where the script resides.
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Destination directory.
DST=$SCRIPT_DIR

# URL of VCF to fetch.
URL=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz

# Regions to fetch
REGIONS="17:41,191,312-41,282,500"  # BRCA1 +/- 5kbp

# Download ClinVar excerpt --------------------------------------------------

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" ERR EXIT

rm -rf $DST/clinvar.excerpt.vcf*

tabix https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz 17:41,191,312-41,282,500 --print-header \
| bgzip -c \
> $TMPDIR/clinvar.excerpt.vcf.gz
tabix -f $TMPDIR/clinvar.excerpt.vcf.gz

bcftools annotate -x ID,INFO $TMPDIR/clinvar.excerpt.vcf.gz \
> $DST/clinvar.excerpt.vcf

bgzip -c $DST/clinvar.excerpt.vcf \
> $DST/clinvar.excerpt.vcf.gz

tabix -f $DST/clinvar.excerpt.vcf.gz
