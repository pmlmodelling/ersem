#!/bin/bash
BRANCH="master"
CMAKE_FLAG=""
while getopts ":b:f:" flag; do
    case "${flag}" in
        b) BRANCH=${OPTARG};;
        f) CMAKE_FLAG=${OPTARG};;
    esac
done

CPU="$(nproc)"

echo "Cloning ERSEM"
git clone https://github.com/pmlmodelling/ersem.git

echo "Cloning FABM"
git clone https://github.com/fabm-model/fabm.git

echo "Cloning GOTM"
git clone https://github.com/gotm-model/code.git gotm
cd gotm && git checkout v6.0 && git submodule update --init --recursive && cd ..

echo "Checking out branch: $BRANCH"
cd ersem && git checkout $BRANCH && cd ..

echo "Building GOTM-FABM-ERSEM"
mkdir build && cd build
cmake ../gotm -DFABM_BASE=../fabm -DFABM_ERSEM_BASE=../ersem $CMAKE_FLAG
make install -j $CPU
