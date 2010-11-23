#!/bin/sh
#$ -A <project>
#$ -V
#$ -cwd
#$ -N wat-cpu
#$ -q normal
#$ -P gpgpu
#$ -pe 8way 128
#$ -l h_rt=00:30:00

ibrun tacc_affinity ./lmp_longhorn-omp -in in.wat-cpu -var num 10 -log log.wat-cpu

