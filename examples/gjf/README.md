# LAMMPS GJF-2GJ THERMOSTAT EXAMPLE

## GJF-2GJ THERMOSTAT

This directory contains the ingredients to run an NVT simulation using the GJF-2GJ thermostat.

Example:
```
NP=4 #number of processors
mpirun -np $NP lmp_mpi -in.argon -out.argon
```

## Required LAMMPS packages: MOLECULE package
