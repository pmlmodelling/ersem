#!/bin/bash
BRANCH=${1:-master}

echo "Cloning ERSEM"
git clone https://github.com/pmlmodelling/ersem.git

echo "Cloning FABM"
git clone https://github.com/fabm-model/fabm.git

echo "Checking out branch: $BRANCH"
cd ersem && git checkout $BRANCH && cd ..

echo "Moving ERSEM config files for pyFABM"
cp github-actions/pyfabm-ersem/setup.cfg fabm

echo "Building PyFABM-ERSEM"
cd fabm
python -m pip install .
