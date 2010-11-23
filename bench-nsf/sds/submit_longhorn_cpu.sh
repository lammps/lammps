#!/bin/sh
#$ -A TG-DMR100081
#$ -V
#$ -cwd
#$ -N sds-cpu
#$ -q normal
#$ -P gpgpu
#$ -pe 8way 128
#$ -l h_rt=00:30:00

ibrun tacc_affinity ./lmp_longhorn-omp -in in.sds-cpu -log log.sds-cpu

