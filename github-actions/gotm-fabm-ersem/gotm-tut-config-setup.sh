#!/bin/bash

echo "Cloning config repo"
git clone https://github.com/pmlmodelling/gotm-ersem-setups.git

cd gotm-ersem-setups/L4


echo "Running GOTM with repo configuration"
./gotm
