#!/bin/bash -e

lmpbin=$1
if [ ! -f $lmpbin ]; then
    echo "LAMMPS binary '$lmpbin' is not a file"
    exit 1
fi


rm -f *.csv
$lmpbin < ref.in > /dev/null
$lmpbin < ffield.in > /dev/null
$lmpbin < ffield_flip.in > /dev/null
for printout in ref.csv ffield.csv ffield_flip.csv; do
    if [ ! -f $printout ]; then
        echo "$printout missing after running LAMMPS"
        exit 1
    fi
#    cat $printout
done
python3 check.py 
rm *.csv log.lammps*
