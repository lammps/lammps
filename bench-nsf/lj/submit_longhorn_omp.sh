#!/bin/sh
#$ -A <project>
#$ -V
#$ -cwd
#$ -N lj-omp
#$ -q normal
#$ -P gpgpu
#$ -pe 4way 128
#$ -l h_rt=00:30:00

OMP_NUM_THREADS=2
export OMP_NUM_THREADS
ibrun tacc_affinity ./lmp_longhorn-omp -in in.lj-omp -var len 66.0 -log log.lj-omp

