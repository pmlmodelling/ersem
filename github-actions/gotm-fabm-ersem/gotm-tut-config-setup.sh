#!/bin/bash

echo "Cloning config repo"
git clone https://github.com/pmlmodelling/ersem-setups.git

cd ersem-setups/L4

echo "Running GOTM with repo configuration"
~/local/gotm/bin/gotm
