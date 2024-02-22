#!/bin/bash

echo "Cloning config repo"
git clone https://github.com/pmlmodelling/ersem-setups.git

cd ersem-setups/0d-aquarium

echo "Running FABM with repo configuration"
~/local/fabm/0d/bin/fabm0d -y ../../ersem/testcases/fabm-ersem-15.06-L4-ben-docdyn-iop.yaml
