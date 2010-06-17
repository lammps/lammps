#!/bin/sh

# user config
lmp=/home/peter/Repository/lammps-awesome/src/lmp_openmpi-omp
job=melt
##################

for m in 2 4 8 ; do \
  mpirun -x OMP_NUM_THREADS=1 --mca mpi_paffinity_alone 1 -npernode ${m} ${lmp} -in in.${job}     -log log.${job}-${m}m-plain
  mpirun -x OMP_NUM_THREADS=1 --mca mpi_paffinity_alone 1 -npernode ${m} ${lmp} -in in.${job}-omp -log log.${job}-${m}m-1t
done
m=4; t=2; 
mpirun -x OMP_NUM_THREADS=${t} --mca mpi_paffinity_alone 0 -npernode ${m} ${lmp} -in in.${job}-omp -log log.${job}-${m}m-${t}t
m=2; t=2; 
mpirun -x OMP_NUM_THREADS=${t} --mca mpi_paffinity_alone 0 -npernode ${m} ${lmp} -in in.${job}-omp -log log.${job}-${m}m-${t}t
m=2; t=4; 
mpirun -x OMP_NUM_THREADS=${t} --mca mpi_paffinity_alone 0 -npernode ${m} ${lmp} -in in.${job}-omp -log log.${job}-${m}m-${t}t
m=1; t=2; 
mpirun -x OMP_NUM_THREADS=${t} --mca mpi_paffinity_alone 0 -npernode ${m} ${lmp} -in in.${job}-omp -log log.${job}-${m}m-${t}t
m=1; t=4; 
mpirun -x OMP_NUM_THREADS=${t} --mca mpi_paffinity_alone 0 -npernode ${m} ${lmp} -in in.${job}-omp -log log.${job}-${m}m-${t}t
m=1; t=8; 
mpirun -x OMP_NUM_THREADS=${t} --mca mpi_paffinity_alone 0 -npernode ${m} ${lmp} -in in.${job}-omp -log log.${job}-${m}m-${t}t
