#! /bin/bash -e

#run the LAMMPS simulation (needs a current LAMMPS version compiled with the user pair_style sw/table)
lmp < spce.in > spce.out

