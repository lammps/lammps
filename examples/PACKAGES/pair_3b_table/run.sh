#! /bin/bash -e

#run the LAMMPS simulation of the first part of the tutorial (needs a current LAMMPS version compiled with the user pair_style 3b/table)
lmp < spce.in > spce.out

#run the LAMMPS simulation of the second part of the tutorial (needs a current LAMMPS version compiled with the user pair_style 3b/table)
lmp < spce_2.in > spce_2.out

