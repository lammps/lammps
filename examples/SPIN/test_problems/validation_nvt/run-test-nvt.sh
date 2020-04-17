#!/bin/bash

# clean old res
rm res_*.dat

# compute NVT Spin -> Lattice 
./../../../../src/lmp_serial -in in.spin.nvt_spin 
in="$(grep -n Step log.lammps | awk -F ':' '{print $1}')"
en="$(grep -n Loop log.lammps | awk -F ':' '{print $1}')"
in="$(echo "$in+1" | bc -l)"
en="$(echo "$en-$in" | bc -l)"
tail -n +$in log.lammps | head -n $en > res_nvt_spin.dat

# compute NVT Lattice -> Spin
./../../../../src/lmp_serial -in in.spin.nvt_lattice 
in="$(grep -n Step log.lammps | awk -F ':' '{print $1}')"
en="$(grep -n Loop log.lammps | awk -F ':' '{print $1}')"
in="$(echo "$in+1" | bc -l)"
en="$(echo "$en-$in" | bc -l)"
tail -n +$in log.lammps | head -n $en > res_nvt_lattice.dat

# plot results
python3 plot_nvt.py res_nvt_spin.dat res_nvt_lattice.dat
