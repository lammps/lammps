#!/bin/bash
# sample launch scripts


# TCP, running LAMMPS on one proc

python driver.py -mdi "-name driver -role DRIVER -method TCP -port 8021" &
../../../src/lmp_mdi -mdi "-name LAMMPS -role ENGINE -method TCP -port 8021 -hostname localhost" -in lammps.in > lammps.out &
wait


# TCP, running LAMMPS on two procs

python driver.py -mdi "-name driver -role DRIVER -method TCP -port 8021" &
mpiexec -n 2 ../../../src/lmp_mdi -mdi "-name LAMMPS -role ENGINE -method TCP -port 8021 -hostname localhost" -in lammps.in > lammps.out &
wait

# MPI, running LAMMPS on two procs

mpiexec -n 1 python driver.py -mdi "-name driver -role DRIVER -method MPI" : \
        -n 2 ../../../src/lmp_mdi -mdi "-name LAMMPS -role ENGINE -method MPI" -in lammps.in > lammps.out
