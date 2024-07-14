#!/usr/bin/bash

# Will install into ~/.local/share/protoc, so make sure to add the following
# to your PATH: ~/.local/share/protoc/bin
#
# Will go into ./utils/var for cloning/building.

set -x
set -euo pipefail

RELEASE=${RELEASE-27.2}
ARCH=${ARCH-linux-x86_64}
PREFIX=${PREFIX-$HOME/.local/share/protoc}

wget -O /tmp/protoc-${RELEASE}-${ARCH}.zip \
    https://github.com/protocolbuffers/protobuf/releases/download/v${RELEASE}/protoc-${RELEASE}-${ARCH}.zip

mkdir -p $PREFIX
cd $PREFIX
unzip /tmp/protoc-${RELEASE}-${ARCH}.zip
