#!/bin/bash

tempi=0.0
tempf=20.0

rm res_*.dat

# compute Lammps 
N=20
for (( i=0; i<$N; i++ ))
do
  temp="$(echo "$tempi+$i*($tempf-$tempi)/$N" | bc -l)"
  sed s/temperature/${temp}/g test-prec-spin.template > \
    test-prec-spin.in
  ./../../../../src/lmp_serial -in test-prec-spin.in 
  Hz="$(tail -n 1 average_spin | awk -F " " '{print $3}')"
  sz="$(tail -n 1 average_spin | awk -F " " '{print $5}')"
  en="$(tail -n 1 average_spin | awk -F " " '{print $6}')"
  echo $temp $Hz $sz $en >> res_lammps.dat
done

# compute Langevin
python3 -m langevin.py > res_langevin.dat

# plot results
python3 -m plot_precession.py res_lammps.dat res_langevin.dat
