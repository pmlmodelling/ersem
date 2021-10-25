#!/bin/bash

echo "Cloning config repo"
git clone https://github.com/pmlmodelling/ersem-setups.git

cd ersem-setups/L4

echo "Running GOTM with repo configuration"
for i in {1..5}
do
    ~/local/gotm/bin/gotm
    gprof -b -Q  ~/local/gotm/bin/gotm gmon.out | awk -v OFS=',' '{print $1,$2,$3,$4,$5,$6,$7}' | sed -e '1,3d' > gprof_$i.csv
done
