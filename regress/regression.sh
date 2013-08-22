#!/bin/bash
cd /code/lammps-atc/regression
./benchmark.py 4 12 min dipole >& latest
fail=`grep -c FAIL latest`
addrs="rjones@sandia.gov jatempl@sandia.gov jzimmer@sandia.gov sjplimp@sandia.gov akohlmey@gmail.com" 
if [ $fail == 0 ] ; then
  mhmail  $addrs -subject "LAMMPS regression passes" <  latest
else
  mhmail  $addrs -subject "LAMMPS regression $fail tests failed" < latest
fi

