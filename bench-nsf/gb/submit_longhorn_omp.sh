#!/bin/sh
#$ -A TG-DMR100081
#$ -V
#$ -cwd
#$ -N gb-omp
#$ -q normal
#$ -P gpgpu
#$ -pe 2way 128
#$ -l h_rt=00:30:00

OMP_NUM_THREADS=4
export OMP_NUM_THREADS
ibrun tacc_affinity ./lmp_longhorn-omp -in in.gb-omp -var len 40.0 -log log.gb-omp

