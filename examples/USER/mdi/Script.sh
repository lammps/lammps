#!/bin/bash
# sample launch scripts

# TCP, running LAMMPS on one proc

python driver.py -mdi "-name driver -role DRIVER -method TCP -port 8021" &
../../../src/lmp_mdi -mdi "-name LAMMPS -role ENGINE -method TCP -port 8021 -hostname localhost" -in lammps.in > lammps.out &
wait


# MPI, running LAMMPS on one proc

mpirun -n 1 python driver.py -mdi "-name driver -role DRIVER -method MPI" : -n 1 ../../../src/lmp_mdi -mdi "-name LAMMPS -role ENGINE -method MPI"
