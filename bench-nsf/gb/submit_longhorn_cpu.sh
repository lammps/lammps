#!/bin/sh
#$ -A TG-DMR100081
#$ -V
#$ -cwd
#$ -N gb-cpu
#$ -q normal
#$ -P gpgpu
#$ -pe 8way 128
#$ -l h_rt=00:30:00

ibrun tacc_affinity ./lmp_longhorn -in in.gb-cpu -var len 40.0 -log log.gb-cpu

