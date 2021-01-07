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

banner "Checking out branch: $BRANCH"
cd ersem && git checkout $BRANCH && cd ..

banner "Installing numpy and wheel"
python -m pip install wheel numpy

banner "Building FABM-ERSEM"
mkdir build && cd build
cmake ../fabm/src/drivers/python -DFABM_ERSEM_BASE=../ersem
make install
