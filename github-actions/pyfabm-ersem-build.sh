#!/bin/bash
BRANCH=${1:-master}

echo "Cloning ERSEM"
git clone https://github.com/pmlmodelling/ersem.git

echo "Cloning FABM"
git clone https://github.com/fabm-model/fabm.git

echo "Checking out branch: $BRANCH"
cd ersem && git checkout $BRANCH && cd ..

echo "Installing numpy and wheel"
python -m pip install wheel numpy

echo "Building FABM-ERSEM"
mkdir build && cd build
cmake ../fabm/src/drivers/python -DFABM_ERSEM_BASE=../ersem
make install
