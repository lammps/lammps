#!/bin/bash -e

lmpbin=$1
if [ ! -f $lmpbin ]; then
    echo "LAMMPS binary '$lmpbin' is not a file"
    exit 1
fi

ref_out="plate_cap.csv"
if [ ! -f $ref_out ]; then
    echo "Generating reference data"
    python3 plate_cap.py > $ref_out
fi

echo "Running Lammps inputs"
rm -rf madelung.txt && touch madelung.txt
for file in in.*; do
    printf "\n$file\n" >> madelung.txt
    rm -f out.csv inv.csv vec.csv 
    $lmpbin -i $file &> /dev/null
    python3 eval.py $ref_out out.csv inv.csv vec.csv
done
cat madelung.txt
