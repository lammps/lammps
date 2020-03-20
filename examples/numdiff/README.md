# LAMMPS FIX NUMDIFF EXAMPLE

## Numerical Difference Fix

This directory contains the ingredients to run an NVE simulation using the numerical difference fix and calculate error in forces.

Example:
```
NP=4 #number of processors
mpirun -np $NP lmp_mpi -in.numdiff
```

## Required LAMMPS packages: MOLECULE package
