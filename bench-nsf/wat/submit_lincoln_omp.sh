#!/bin/sh
#PBS -A epe 
#PBS -l nodes=16:ppn=8
#PBS -q lincoln 
#PBS -l walltime=1:00:00
#PBS -N wat-cpu
#PBS
#

cd $PBS_O_WORKDIR

mpirun -npernode 8 --mca mpi_paffinity_alone 0 -x OMP_NUM_THREADS=1 ./lmp_abe -in in.wat-omp -var num 10 -log log.wat-omp

mpirun -npernode 4 --mca mpi_paffinity_alone 0 -x OMP_NUM_THREADS=2 ./lmp_abe -in in.wat-omp -var num 10 -log log.wat-omp2

mpirun -npernode 2 --mca mpi_paffinity_alone 0 -x OMP_NUM_THREADS=4 ./lmp_abe -in in.wat-omp -var num 10 -log log.wat-omp3

