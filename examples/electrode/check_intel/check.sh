#!/bin/bash -e

lmpbin=$1
if [ ! -f $lmpbin ]; then
    echo "LAMMPS binary '$lmpbin' is not a file"
    exit 1
fi


rm -f *.csv
$lmpbin -i ref.in #> /dev/null
mv ref.csv ref_nointel.csv
$lmpbin -i ref.in -pk intel 0 omp 1 mode double -sf intel #> /dev/null
mv ref.csv ref_intel.csv

python3 check.py 
#rm *.csv log.lammps*
