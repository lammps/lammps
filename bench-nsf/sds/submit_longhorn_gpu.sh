#!/bin/sh
#$ -A TG-DMR100081
#$ -V
#$ -cwd
#$ -N sds-gpu
#$ -q normal
#$ -P gpgpu
#$ -pe 2way 128
#$ -l h_rt=00:30:00

ibrun ./lmp_longhorn -in in.sds-gpu -log log.sds-gpu

