#!/bin/bash

# clean old res
rm res_*.dat

### compute NVT Spin -> Lattice 

# test standard Lammps 
./../../../../src/lmp_serial -in in.spin.nvt_spin 

# test spin/kk with Kokkos Lammps
# mpirun -np 1 ../../../../src/lmp_kokkos_mpi_only \
#   -k on -sf kk -in in.spin.nvt_spin

# extract data from Lammps run
in="$(grep -n Step log.lammps | awk -F ':' '{print $1}')"
en="$(grep -n Loop log.lammps | awk -F ':' '{print $1}')"
in="$(echo "$in+1" | bc -l)"
en="$(echo "$en-$in" | bc -l)"
tail -n +$in log.lammps | head -n $en > res_nvt_spin.dat

### compute NVT Lattice -> Spin

# test standard Lammps 
./../../../../src/lmp_serial -in in.spin.nvt_lattice 

# test spin/kk with Kokkos Lammps
# mpirun -np 1 ../../../../src/lmp_kokkos_mpi_only \
#   -k on -sf kk -in in.spin.nvt_lattice

# extract data from Lammps run
in="$(grep -n Step log.lammps | awk -F ':' '{print $1}')"
en="$(grep -n Loop log.lammps | awk -F ':' '{print $1}')"
in="$(echo "$in+1" | bc -l)"
en="$(echo "$en-$in" | bc -l)"
tail -n +$in log.lammps | head -n $en > res_nvt_lattice.dat

# plot results
python3 plot_nvt.py res_nvt_spin.dat res_nvt_lattice.dat
