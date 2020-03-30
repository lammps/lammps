#!/bin/bash

# clean old res
rm res_*.dat

# compute Lammps 
./../../../../src/lmp_serial \
  -in in.spin.iron-nve 
in="$(grep -n Step log.lammps | awk -F ':' '{print $1}')"
en="$(grep -n Loop log.lammps | awk -F ':' '{print $1}')"
in="$(echo "$in+1" | bc -l)"
en="$(echo "$en-$in" | bc -l)"
tail -n +$in log.lammps | head -n $en > res_lammps.dat

# plot results
python3 -m plot_nve.py res_lammps.dat res_llg.dat
