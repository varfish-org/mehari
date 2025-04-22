#!/usr/bin/bash

set -x
set -euo pipefail

# Interpreted environment variables.
#
#   PATH_DB         -- path to the database directory containing,
#                      e.g., `grch3{7,8}/*.zst`.
#                      default: /data/mehari/db
#   HTTP_HOST       -- host to listen on
#                      default: 0.0.0.0
#   HTTP_PORT       -- port
#                      default: 8080

PATH_DB=${PATH_DB-/data/mehari/db}
HTTP_HOST=${HTTP_HOST-0.0.0.0}
HTTP_PORT=${HTTP_PORT-8080}

PATH_TRANSCRIPTS_37=$PATH_DB/grch37/seqvars/txs.bin.zst
PATH_TRANSCRIPTS_38=$PATH_DB/grch38/seqvars/txs.bin.zst
PATH_FREQUENCIES_37=$PATH_DB/grch37/seqvars/freqs/rocksdb
PATH_FREQUENCIES_38=$PATH_DB/grch38/seqvars/freqs/rocksdb
PATH_CLINVAR_37=$PATH_DB/grch37/seqvars/clinvar/rocksdb
PATH_CLINVAR_38=$PATH_DB/grch38/seqvars/clinvar/rocksdb

first=${1-}

if [ "$first" == exec ]; then
  shift
  exec "$@"
else
  exec \
    mehari \
    server run \
      $(test -e "$PATH_TRANSCRIPTS_37" && echo --transcripts "$PATH_TRANSCRIPTS_37") \
      $(test -e "$PATH_TRANSCRIPTS_38" && echo --transcripts "$PATH_TRANSCRIPTS_38") \
      $(test -e "$PATH_FREQUENCIES_37" && echo --frequencies "$PATH_FREQUENCIES_37") \
      $(test -e "$PATH_FREQUENCIES_38" && echo --frequencies "$PATH_FREQUENCIES_38") \
      $(test -e "$PATH_CLINVAR_37" && echo --clinvar "$PATH_CLINVAR_37") \
      $(test -e "$PATH_CLINVAR_38" && echo --clinvar "$PATH_CLINVAR_38") \
      --listen-host "$HTTP_HOST" \
      --listen-port "$HTTP_PORT"
fi

exit $?
