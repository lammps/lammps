#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 1 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.1

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 2 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.2

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 3 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.3

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 4 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.4

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 5 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.5

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 6 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.6

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 7 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.7

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 8 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.8

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 9 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.9

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 10 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.10

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 11 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.11

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 12 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.12

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 13 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.13

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 14 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.14

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 15 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.15

mpirun -np 1 lmp_kokkos_cuda -k on g 1 t 16 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.1.16

mpirun -np 2 lmp_kokkos_cuda -k on g 2 t 1 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.2.1

mpirun -np 2 lmp_kokkos_cuda -k on g 2 t 2 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.2.2

mpirun -np 2 lmp_kokkos_cuda -k on g 2 t 3 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.2.3

mpirun -np 2 lmp_kokkos_cuda -k on g 2 t 4 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.2.4

mpirun -np 2 lmp_kokkos_cuda -k on g 2 t 5 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.2.5

mpirun -np 2 lmp_kokkos_cuda -k on g 2 t 6 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.2.6

mpirun -np 2 lmp_kokkos_cuda -k on g 2 t 7 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.2.7

mpirun -np 2 lmp_kokkos_cuda -k on g 2 t 8 -sf kk -pk kokkos binsize 2.8 comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.cuda.128K.2.8
