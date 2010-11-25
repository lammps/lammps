#!/bin/sh
#$ -A TG-TRA100009
#$ -V
#$ -cwd
#$ -N sds-omp
#$ -q development
#$ -pe 4way 128
#$ -l h_rt=00:30:00

OMP_NUM_THREADS=4
export OMP_NUM_THREADS
ibrun tacc_affinity ./lmp_ranger -in in.sds-omp -log log.sds-omp

