#!/bin/bash -e

lmpbin=$1
if [ ! -f $lmpbin ]; then
    echo "LAMMPS binary '$lmpbin' is not a file"
    exit 1
fi

for file in in.*; do
    echo "$file"
    echo "1 proc"
    $lmpbin -i $file > /dev/null
    grep '  1.2  ' log.lammps
    echo "2 procs"
    mpirun -np 2 $lmpbin -i $file > /dev/null
    grep '  1.2  ' log.lammps
done
