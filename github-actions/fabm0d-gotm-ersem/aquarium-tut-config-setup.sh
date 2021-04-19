#!/bin/bash

echo "Cloning config repo"
git clone https://github.com/pmlmodelling/ersem-setups.git

cd ersem-setups/0d-aquarium

echo "Running FABM with repo configuration"
~/local/fabm/0d/bin/fabm0d -y ../L4/fabm.yaml
