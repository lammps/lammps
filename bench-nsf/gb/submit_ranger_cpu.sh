#!/bin/sh
#$ -A TG-TRA100009
#$ -V
#$ -cwd
#$ -N gb-cpu
#$ -q development
#$ -pe 16way 128
#$ -l h_rt=00:30:00

OMP_NUM_THREADS=1
export OMP_NUM_THREADS=1
ibrun tacc_affinity ./lmp_ranger -in in.gb-cpu -var len 40.0 -log log.gb-cpu

