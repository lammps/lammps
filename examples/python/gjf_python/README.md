# LAMMPS GJF-2GJ THERMOSTAT EXAMPLE W/ PYTHON

## GJF-2GJ THERMOSTAT

This directory contains a python script to run NVT simulations using the GJF-2GJ thermostat.
The script will vary the timestep and write thermodynamic output to screen.
This script has True/False options to change how you would like to dump/write your output.

Example:
```
NP=4 #number of processors
mpirun -np $NP python gjf.py
```

## Required LAMMPS packages: MOLECULE package
## LAMMPS COMPILE MODE: SHLIB
## LAMMPS OPTIONAL INSTALL: make install-python
## Required Python packages: mpi4py 
