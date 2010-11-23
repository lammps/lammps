#!/bin/sh
#PBS -A <project>
#PBS -l nodes=16:ppn=8
#PBS -q lincoln 
#PBS -l walltime=1:00:00
#PBS -N sds-gpu
#PBS
#

cd $PBS_O_WORKDIR

mpirun -npernode 2 --mca mpi_paffinity_alone 1 ./lmp_lincoln -in in.sds-gpu -log log.sds-gpu

