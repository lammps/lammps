#!/bin/sh
#PBS -A <project>
#PBS -l nodes=16:ppn=8
#PBS -q lincoln 
#PBS -l walltime=1:00:00
#PBS -N lj-gpu
#PBS
#

cd $PBS_O_WORKDIR
mpirun -npernode 2 --mca mpi_paffinity_alone 1 -x OMP_NUM_THREADS=1 ./lmp_lincoln -in in.lj-gpu -log log.lj-gpu

