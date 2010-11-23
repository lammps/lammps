#!/bin/sh 
#PBS -A TG-TRA100009
#PBS -q normal
#PBS -N wat-all
#PBS -l nodes=16:ppn=8
#PBS -j oe
#PBS -l walltime=03:00:00

set -v
set -x

cd $PBS_O_WORKDIR

mpirun=/usr/local/openmpi-1.3.2-gnu/bin/mpirun

$mpirun -npernode 7 -x OMP_NUM_THREADS=1 ./lmp_abe -in in.wat-cpu -var num 10 -log log.wat-cpu-7m1t
$mpirun -npernode 7 -x OMP_NUM_THREADS=1 ./lmp_abe -in in.wat-omp -var num 10 -log log.wat-omp-7m1t
$mpirun -npernode 4 -x OMP_NUM_THREADS=1 ./lmp_abe -in in.wat-cpu -var num 10 -log log.wat-cpu-4m1t
$mpirun -npernode 4 -x OMP_NUM_THREADS=2 ./lmp_abe -in in.wat-omp -var num 10 -log log.wat-omp-4m2t
$mpirun -npernode 2 -x OMP_NUM_THREADS=1 ./lmp_abe -in in.wat-cpu -var num 10 -log log.wat-cpu-2m1t
$mpirun -npernode 2 -x OMP_NUM_THREADS=4 ./lmp_abe -in in.wat-omp -var num 10 -log log.wat-omp-2m4t

