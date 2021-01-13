#!/bin/bash
BRANCH=${1:-master}

echo "Cloning ERSEM"
git clone https://github.com/pmlmodelling/ersem.git

echo "Cloning FABM"
git clone https://github.com/fabm-model/fabm.git

echo "Cloning GOTM"
git clone https://github.com/gotm-model/code.git gotm
cd gotm && git submodule update --init --recursive && cd ..

echo "Checking out branch: $BRANCH"
cd ersem && git checkout $BRANCH && cd ..

echo "Installing netCDF"
sudo apt install libnetcdff-dev

echo "Building FABM-GOTM-ERSEM"
mkdir build && cd build
cmake ../fabm/scr/drivers/0d -DGOTM_BASE=../gotm -DFABM_ERSEM_BASE=../ersem
make install
