#!/bin/bash -e

lmpbin=$1
if [ ! -f $lmpbin ]; then
    echo "LAMMPS binary '$lmpbin' is not a file"
    exit 1
fi


rm -f *.csv
$lmpbin < onestep.in > /dev/null
$lmpbin < twostep.in > /dev/null
for matrix in twostep.csv onestep.csv; do
    if [ ! -f $matrix ]; then
        echo "$matrix missing after running LAMMPS"
        exit 1
    fi
done
python3 check.py 
rm *.csv log.lammps
