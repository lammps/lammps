#!/bin/sh
#$ -A TG-DMR100081
#$ -V
#$ -cwd
#$ -N lj-cpu
#$ -q normal
#$ -P gpgpu
#$ -pe 8way 128
#$ -l h_rt=00:30:00

ibrun tacc_affinity ./lmp_longhorn -in in.lj-cpu -var len 66.0 -log log.lj-cpu

