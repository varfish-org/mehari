#!/usr/bin/bash

set -x
set -euo pipefail

# Interpreted environment variables.
#
#   PATH_DB         -- path to the database directory containing,
#                      e.g., `grch3{7,8}/*.zst`.
#                      default: /data/hpo
#   HTTP_HOST       -- host to listen on
#                      default: 0.0.0.0
#   HTTP_PORT       -- port
#                      default: 8080

PATH_DB=${PATH_DB-/data/mehari}
HTTP_HOST=${HTTP_HOST-0.0.0.0}
HTTP_PORT=${HTTP_PORT-8080}

first=${1-}

if [ "$first" == exec ]; then
  shift
  exec "$@"
else
  exec \
    mehari \
    server run \
      --path-db "$PATH_DB" \
      --listen-host "$HTTP_HOST" \
      --listen-port "$HTTP_PORT"
fi

exit $?
