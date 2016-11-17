#!/bin/sh

export MPI_IMPLEMENTATION=openmpi
PROC="4x1"

mpirun -np 4 /project/normoliq/dstelter/grem-feature/lammps/src/lmp_openmpi.scc -p $PROC -in in.gREM-temper
