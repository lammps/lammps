#!/bin/bash -e

lmpbin=$1
if [ ! -f $lmpbin ]; then
    echo "LAMMPS binary '$lmpbin' is not a file"
    exit 1
fi


rm -f *.csv
$lmpbin -i ref.in > /dev/null
$lmpbin -i etypes.in > /dev/null
for printout in ref.csv etypes.csv; do
    if [ ! -f $printout ]; then
        echo "$printout missing after running LAMMPS"
        exit 1
    fi
    #cat $printout
done
python3 check.py 
rm *.csv log.lammps*
