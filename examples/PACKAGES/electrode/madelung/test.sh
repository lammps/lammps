#!/bin/bash -e

lmpbin=$1
if [ ! -f $lmpbin ]; then
    echo "LAMMPS binary '$lmpbin' is not a file"
    exit 1
fi

ref_out="plate_cap.csv"
ref_mix_out="plate_cap_eta_mix.csv"
if [ ! -f $ref_out ] || [ ! -f $ref_mix_out ]; then
    echo "Generating reference data"
    python3 plate_cap.py
fi

echo "Running Lammps inputs"
# w/o eta mixing
rm -rf madelung.txt && touch madelung.txt
for file in in.eta in.ewald-ew3dc in.ewald-ew2d in.pppm-ew3dc in.cg; do
    printf "\n$file\n" >> madelung.txt
    rm -f out.csv inv.csv vec.csv 
    $lmpbin -i $file &> /dev/null
    python3 eval.py $ref_out out.csv inv.csv vec.csv
done

# with eta mixing
for file in in.eta_mix in.eta_cg; do
    printf "\n$file\n" >> madelung.txt
    rm -f out.csv inv.csv vec.csv 
    $lmpbin -i $file &> /dev/null
    python3 eval.py $ref_mix_out out.csv inv.csv vec.csv
done
cat madelung.txt
