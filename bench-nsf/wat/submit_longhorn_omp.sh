#!/bin/sh
#$ -A TG-DMR100081
#$ -V
#$ -cwd
#$ -N wat-omp
#$ -q normal
#$ -P gpgpu
#$ -pe 4way 128
#$ -l h_rt=00:30:00

OMP_NUM_THREADS=2
export OMP_NUM_THREADS
ibrun tacc_affinity ./lmp_longhorn-omp -in in.wat-omp -var num 10 -log log.wat-omp

