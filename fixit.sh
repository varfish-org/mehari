#!/usr/bin/bash

set -euo pipefail
set -x

for s in $(find . -name 'bootstrap.sh' | sort); do
    bash $s
done
