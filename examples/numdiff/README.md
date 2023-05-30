# LAMMPS NUMDIFF EXAMPLES FOR FORCES, VIRIAL, and BORN MATRIX

## Numerical Difference Fixes and Computes

This directory contains the input script for an NVE simulation with
fix numdiff, fix numdiff/virial, and compute born/matrix numdiff.
In each cases, results are compared to exact analytic expressions
and the small relative differences are reported.

Example:
```
NP=4 #number of processors
mpirun -np $NP lmp_mpi -in.numdiff
```

## Required LAMMPS packages: EXTRA-FIX, EXTRA-COMPUTE
