#!/bin/bash
#SBATCH --time=00:05:00 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=langevin
#SBATCH --partition=low
#SBATCH --output=vtest1-new.out
#SBATCH --error=vtest1-new.err

srun lmp.exe -in in.gjf.vfull