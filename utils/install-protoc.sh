#!/usr/bin/bash

# Will install into ~/.local/share/protoc, so make sure to add the following
# to your PATH: ~/.local/share/protoc/bin
#
# Will go into ./utils/var for cloning/building.

set -x
set -euo pipefail

CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX-$HOME/.local/share/protoc}

mkdir -p utils/var
cd utils/var

apt-get update
apt-get install -y git cmake build-essential

if [[ ! -e protobuf ]]; then
    git clone https://github.com/protocolbuffers/protobuf.git
fi
cd protobuf
git submodule update --init --recursive

cmake . -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX
make -j 8 install
