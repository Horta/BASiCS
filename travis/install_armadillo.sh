#!/usr/bin/env bash

set -e

wget https://sourceforge.net/projects/arma/files/armadillo-7.600.2.tar.xz
tar xvfJ armadillo-7.600.2.tar.xz

pushd armadillo-7.600.2
# perl -i -pe 's/ARMA_USE_WRAPPER true/ARMA_USE_WRAPPER false/gc' CMakeLists.txt
# perl -i -pe 's/armadillo PROPERTIES/armadillo PROPERTIES POSITION_INDEPENDENT_CODE ON/gc' CMakeLists.txt
cmake . && make && sudo make install
popd
