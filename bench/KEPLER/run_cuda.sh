#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -N 1 lmp_cuda_double -c on -sf cuda -pk cuda 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.cuda.double.128K.1

mpirun -N 2 lmp_cuda_double -c on -sf cuda -pk cuda 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.cuda.double.128K.2

mpirun -N 1 lmp_cuda_mixed -c on -sf cuda -pk cuda 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.cuda.mixed.128K.1

mpirun -N 2 lmp_cuda_mixed -c on -sf cuda -pk cuda 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.cuda.mixed.128K.2

mpirun -N 1 lmp_cuda_single -c on -sf cuda -pk cuda 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.cuda.single.128K.1

mpirun -N 2 lmp_cuda_single -c on -sf cuda -pk cuda 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.cuda.single.128K.2
