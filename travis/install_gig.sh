#!/usr/bin/env bash

set -e

wget https://github.com/Horta/gig/archive/v0.0.3.tar.gz
tar xzf v0.0.3.tar.gz

pushd gig-0.0.3
mkdir build && cd build
cmake .. && make && make test && sudo make install
popd
