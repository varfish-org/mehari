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

if [[ "$#" -ne 2 ]]; then
    log "USAGE: bootstrap.sh SEQREPO INSTANCE"
    log ""
    log "Set VERBOSE=1 to increase verbosity."
    exit 1
fi

# path to the directory where the script resides.
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# SeqRepo source directory.
SRC=$1

# SeqRepo instance.
INSTANCE=$2

# Destination directory.
DST=$SCRIPT_DIR

# Import SQLite database ----------------------------------------------------

rm -rf $DST/latest brca1.fasta

seqrepo --root-directory $DST init --instance-name latest

gene=brca1
tx="NM_007294.4 NM_007297.4 NM_007298.3 NM_007299.4 NM_007300.4"

seqrepo --root-directory $SRC export -i $INSTANCE $tx \
| sed -e 's/NCBI://g' -e 's/ refseq:.*//g' \
> $DST/$gene.fasta

seqrepo --root-directory $DST load --instance-name latest --namespace refseq $DST/$gene.fasta
