#!/bin/bash
banner() {
    msg="# $* #"
    edge=$(echo "$msg" | sed 's/./#/g')
    echo -e "\n$edge"
    echo -e "$msg"
    echo -e "$edge\n"
}

BRANCH=${1:-master}

banner "Cloning ERSEM"
git clone https://github.com/pmlmodelling/ersem.git

banner "Cloning FABM"
git clone https://github.com/fabm-model/fabm.git

banner "Cloning GOTM"
git clone https://github.com/gotm-model/code.git gotm
cd gotm && git submodule update --init --recursive && cd ..

banner "Checking out branch: $BRANCH"
cd ersem && git checkout $BRANCH && cd ..

banner "Installing netCDF"
sudo apt install libnetcdff-dev

banner "Building GOTM-FABM-ERSEM"
mkdir build && cd build
cmake ../gotm -DFABM_BASE=../fabm -DFABM_ERSEM_BASE=../ersem
make install
