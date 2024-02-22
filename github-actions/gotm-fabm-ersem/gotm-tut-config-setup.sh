#!/bin/bash

echo "Cloning config repo"
git clone https://github.com/pmlmodelling/ersem-setups.git

cp ersem/testcases/fabm-ersem-15.06-L4-ben-docdyn-iop.yaml ersem-setups/L4/fabm.yaml
cd ersem-setups/L4

echo "Running GOTM with repo configuration"
~/local/gotm/bin/gotm
