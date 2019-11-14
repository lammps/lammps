#!/bin/bash

# set initial and final temperature (K)
tempi=0.0
tempf=2000.0

rm res_*.dat

# run Lammps calculation
N=20
for (( i=0; i<$N; i++ ))
do
  temp="$(echo "$tempi+$i*($tempf-$tempi)/$N" | bc -l)"
  sed s/temperature/${temp}/g bench-exchange-spin.template > \
    bench-exchange-spin.in
  ../../../../src/lmp_serial \
    -in bench-exchange-spin.in 
  Hz="$(tail -n 1 _av_spin | awk -F " " '{print $3}')"
  sm="$(tail -n 1 _av_spin | awk -F " " '{print $5}')"
  en="$(tail -n 1 _av_spin | awk -F " " '{print $6}')"
  echo $temp $Hz $sm $en >> res_lammps.dat
done

# run Langevin function calculation
python3 -m langevin-exchange.py > res_langevin.dat

# plot comparison
python3 -m plot_exchange.py res_lammps.dat res_langevin.dat
