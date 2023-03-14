#!/usr/bin/bash

# Will install into ~/.local/share/flatbuffers, so make sure to add the following
# to your PATH: ~/.local/share/flatbuffers/bin
#
# Will go into ./utils/var for cloning/building.

mkdir -p utils/var
cd utils/var

git clone https://github.com/google/flatbuffers.git
cd flatbuffers
git checkout v22.12.06
cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=$HOME/.local/share/flatbuffers
make
./flattests
make install
flatc --version
