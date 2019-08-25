#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -np 1 lmp_omp -sf omp -pk omp 16 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.omp.128K.1.16

mpirun -np 2 lmp_omp -sf omp -pk omp 8 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.omp.128K.2.8

mpirun -np 4 lmp_omp -sf omp -pk omp 4 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.omp.128K.4.4

mpirun -np 8 lmp_omp -sf omp -pk omp 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.omp.128K.8.2

mpirun -np 16 lmp_omp -sf omp -pk omp 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.omp.128K.16.1
