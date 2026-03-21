#!/usr/bin/env bash

set -x

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd $SCRIPT_DIR

rm -rf 14kb.txt 14kb.txt.gz 14kb.txt.bgz

seq 3000 > 14kb.txt
bgzip -c 14kb.txt >14kb.txt.gz
ln -sr 14kb.txt.gz 14kb.txt.bgz
