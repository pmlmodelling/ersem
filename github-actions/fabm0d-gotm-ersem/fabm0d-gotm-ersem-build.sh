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
git clone --recursive https://github.com/gotm-model/code.git gotm

echo "Checking out branch: $BRANCH"
cd ersem && git checkout $BRANCH && cd ..

echo "Building FABM-GOTM-ERSEM"
mkdir build && cd build
cmake ../fabm/src/drivers/0d -DGOTM_BASE=../gotm -DFABM_ERSEM_BASE=../ersem $CMAKE_FLAG
make install -j $CPU
