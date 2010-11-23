#!/bin/sh
#$ -A TG-DMR100081
#$ -V
#$ -cwd
#$ -N gb-gpu
#$ -q normal
#$ -P gpgpu
#$ -pe 2way 128
#$ -l h_rt=00:30:00

ibrun ./lmp_longhorn -in in.gb-gpu -var len 40.0 -log log.gb-gpu

