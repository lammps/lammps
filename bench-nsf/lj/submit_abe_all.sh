#!/bin/sh 
#PBS -A <project>
#PBS -q normal
#PBS -N lj-all
#PBS -l nodes=16:ppn=8
#PBS -j oe
#PBS -l walltime=02:00:00

set -v
set -x

cd $PBS_O_WORKDIR

mpirun=/usr/local/openmpi-1.3.2-gnu/bin/mpirun

$mpirun -npernode 7 -x OMP_NUM_THREADS=1 ./lmp_abe -in in.lj-cpu -var len 66.0 -log log.lj-cpu-7m1t
$mpirun -npernode 7 -x OMP_NUM_THREADS=1 ./lmp_abe -in in.lj-omp -var len 66.0 -log log.lj-omp-7m1t
$mpirun -npernode 4 -x OMP_NUM_THREADS=1 ./lmp_abe -in in.lj-cpu -var len 66.0 -log log.lj-cpu-4m1t
$mpirun -npernode 4 -x OMP_NUM_THREADS=2 ./lmp_abe -in in.lj-omp -var len 66.0 -log log.lj-omp-4m2t
$mpirun -npernode 2 -x OMP_NUM_THREADS=1 ./lmp_abe -in in.lj-cpu -var len 66.0 -log log.lj-cpu-2m1t
$mpirun -npernode 2 -x OMP_NUM_THREADS=4 ./lmp_abe -in in.lj-omp -var len 66.0 -log log.lj-omp-2m4t

