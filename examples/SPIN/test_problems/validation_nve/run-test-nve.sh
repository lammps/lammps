#!/bin/bash

# clean old res
rm res_*.dat

# test standard Lammps 
../../../../src/lmp_serial -in in.spin.iron-nve 

# test spin/kk with Kokkos Lammps
# mpirun -np 1 ../../../../src/lmp_kokkos_mpi_only \
#   -k on -sf kk -in in.spin.iron-nve

# extract data from Lammps run
in="$(grep -n Step log.lammps | awk -F ':' '{print $1}')"
en="$(grep -n Loop log.lammps | awk -F ':' '{print $1}')"
in="$(echo "$in+1" | bc -l)"
en="$(echo "$en-$in" | bc -l)"
tail -n +$in log.lammps | head -n $en > res_lammps.dat

# plot results
python3 plot_nve.py res_lammps.dat res_llg.dat
