#!/bin/bash
cd /code/lammps-atc/regress
#./benchmark.py 4 12 min melt >& latest
fail=`grep -c FAIL latest`
addresses="rjones@sandia.gov jatempl@sandia.gov jzimmer@sandia.gov sjplimp@sandia.gov pscrozi@sandia.gov akohlmey@gmail.com" 
subject="\"LAMMPS regression $fail tests failed\""
if [ $fail == 0 ] ; then
  subject='"LAMMPS regression passed"'
fi
scp latest vikramarka.ca.sandia.gov:.
#echo $subject
#echo "mhmail $addresses -subject $subject < latest"
ssh vikramarka.ca.sandia.gov "mhmail $addresses -subject $subject < latest"

