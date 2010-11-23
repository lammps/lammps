#!/bin/sh
#PBS -A epe 
#PBS -l nodes=16:ppn=8
#PBS -q lincoln 
#PBS -l walltime=1:00:00
#PBS -N wat-gpu
#PBS
#

cd $PBS_O_WORKDIR

mpirun -npernode 2 --mca mpi_paffinity_alone 1 ./lmp_lincoln -in in.wat-gpu -var num 10 -log log.wat-gpu

