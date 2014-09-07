#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -bind-to none -N 1 ./lmp_omp -v x 8 -v y 8 -v z 8 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.2K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 8 -v y 8 -v z 16 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.4K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 8 -v y 16 -v z 16 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.8K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 16 -v y 16 -v z 16 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.16K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 16 -v y 16 -v z 32 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.32K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 16 -v y 32 -v z 32 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.64K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 32 -v y 32 -v z 32 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.128K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 32 -v y 32 -v z 64 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.256K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 32 -v y 64 -v z 64 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.512K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 64 -v y 64 -v z 64 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.1024K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 64 -v y 64 -v z 128 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.2048K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 64 -v y 128 -v z 128 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.4096K.1.16
mpirun -bind-to none -N 1 ./lmp_omp -v x 128 -v y 128 -v z 128 -v t 100 -v h 16 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.8192K.1.16

mpirun -bind-to socket -N 2 ./lmp_omp -v x 8 -v y 8 -v z 8 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.2K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 8 -v y 8 -v z 16 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.4K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 8 -v y 16 -v z 16 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.8K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 16 -v y 16 -v z 16 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.16K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 16 -v y 16 -v z 32 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.32K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 16 -v y 32 -v z 32 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.64K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 32 -v y 32 -v z 32 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.128K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 32 -v y 32 -v z 64 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.256K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 32 -v y 64 -v z 64 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.512K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 64 -v y 64 -v z 64 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.1024K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 64 -v y 64 -v z 128 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.2048K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 64 -v y 128 -v z 128 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.4096K.2.8
mpirun -bind-to socket -N 2 ./lmp_omp -v x 128 -v y 128 -v z 128 -v t 100 -v h 8 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.8192K.2.8

mpirun -bind-to socket -N 4 ./lmp_omp -v x 8 -v y 8 -v z 8 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.2K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 8 -v y 8 -v z 16 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.4K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 8 -v y 16 -v z 16 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.8K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 16 -v y 16 -v z 16 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.16K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 16 -v y 16 -v z 32 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.32K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 16 -v y 32 -v z 32 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.64K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 32 -v y 32 -v z 32 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.128K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 32 -v y 32 -v z 64 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.256K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 32 -v y 64 -v z 64 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.512K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 64 -v y 64 -v z 64 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.1024K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 64 -v y 64 -v z 128 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.2048K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 64 -v y 128 -v z 128 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.4096K.4.4
mpirun -bind-to socket -N 4 ./lmp_omp -v x 128 -v y 128 -v z 128 -v t 100 -v h 4 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.8192K.4.4

mpirun -bind-to socket -N 8 ./lmp_omp -v x 8 -v y 8 -v z 8 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.2K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 8 -v y 8 -v z 16 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.4K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 8 -v y 16 -v z 16 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.8K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 16 -v y 16 -v z 16 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.16K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 16 -v y 16 -v z 32 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.32K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 16 -v y 32 -v z 32 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.64K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 32 -v y 32 -v z 32 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.128K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 32 -v y 32 -v z 64 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.256K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 32 -v y 64 -v z 64 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.512K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 64 -v y 64 -v z 64 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.1024K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 64 -v y 64 -v z 128 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.2048K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 64 -v y 128 -v z 128 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.4096K.8.2
mpirun -bind-to socket -N 8 ./lmp_omp -v x 128 -v y 128 -v z 128 -v t 100 -v h 2 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.8192K.8.2

mpirun -bind-to core -N 16 ./lmp_omp -v x 8 -v y 8 -v z 8 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.2K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 8 -v y 8 -v z 16 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.4K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 8 -v y 16 -v z 16 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.8K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 16 -v y 16 -v z 16 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.16K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 16 -v y 16 -v z 32 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.32K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 16 -v y 32 -v z 32 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.64K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 32 -v y 32 -v z 32 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.128K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 32 -v y 32 -v z 64 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.256K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 32 -v y 64 -v z 64 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.512K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 64 -v y 64 -v z 64 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.1024K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 64 -v y 64 -v z 128 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.2048K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 64 -v y 128 -v z 128 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.4096K.16.1
mpirun -bind-to core -N 16 ./lmp_omp -v x 128 -v y 128 -v z 128 -v t 100 -v h 1 -sf omp < in.lj.omp
mv log.lammps log.28Jun14.lj.omp.8192K.16.1
