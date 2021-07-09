#!/usr/bin/env bash

mkdir -p build
cd build

export CFLAGS="${CFLAGS} -I${CONDA_PREFIX}/include"
export CXXFLAGS="${CXXFLAGS} -I${CONDA_PREFIX}/include"

# export CXXC='mpicc'

cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_CXX_COMPILER=$CXX

make VERBOSE=1
make install
