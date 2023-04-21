#!/usr/bin/bash

# Will install into ~/.local/share/flatbuffers, so make sure to add the following
# to your PATH: ~/.local/share/flatbuffers/bin
#
# Will go into ./utils/var for cloning/building.

mkdir -p utils/var
cd utils/var

sudo apt-get install g++ git

if [[ ! -e protobuf ]]; then
    git clone https://github.com/protocolbuffers/protobuf.git
fi
cd protobuf
git submodule update --init --recursive

cmake . CMAKE_INSTALL_PREFIX=$HOME/.local/share/protoc
make -j 8 install

mkdir -p ~/.local/share/protoc/bin
cp bazel-bin/protoc ~/.local/share/protoc/bin
