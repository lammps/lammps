# LAMMPS LATTICE DYNAMICS COMMANDS

## DYNAMICAL MATRIX CALCULATOR

This directory contains the ingredients to calculate a dynamical matrix.  

Example:
```
NP=4 #number of processors
mpirun -np $NP lmp_mpi -in in.silicon -out out.silicon
```

To test out a different silicon example:
```
LMP_FILE=amorphous_silicon.lmp
cp lmp_bank/$LMP_FILE ./silicon_input_file.lmp
NP=4 #number of processors
mpirun -np $NP lmp_mpi -in in.silicon -out out.silicon
```

## Requires: MANYBODY and MOLECULE packages
