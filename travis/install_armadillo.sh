#!/usr/bin/env bash

set -e

wget https://sourceforge.net/projects/arma/files/armadillo-7.600.2.tar.xz
tar xvfJ armadillo-7.600.2.tar.xz

pushd armadillo-7.600.2
cmake . && make && sudo make install
popd
