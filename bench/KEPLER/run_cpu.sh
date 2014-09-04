#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -N 1 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.1
mpirun -N 1 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.1
mpirun -N 1 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.1
mpirun -N 1 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.1
mpirun -N 1 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.1
mpirun -N 1 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.1
mpirun -N 1 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.1
mpirun -N 1 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.1
mpirun -N 1 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.1
mpirun -N 1 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.1
mpirun -N 1 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.1
mpirun -N 1 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.1
mpirun -N 1 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.1

mpirun -N 2 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.2
mpirun -N 2 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.2
mpirun -N 2 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.2
mpirun -N 2 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.2
mpirun -N 2 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.2
mpirun -N 2 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.2
mpirun -N 2 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.2
mpirun -N 2 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.2
mpirun -N 2 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.2
mpirun -N 2 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.2
mpirun -N 2 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.2
mpirun -N 2 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.2
mpirun -N 2 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.2

mpirun -N 3 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.3
mpirun -N 3 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.3
mpirun -N 3 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.3
mpirun -N 3 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.3
mpirun -N 3 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.3
mpirun -N 3 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.3
mpirun -N 3 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.3
mpirun -N 3 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.3
mpirun -N 3 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.3
mpirun -N 3 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.3
mpirun -N 3 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.3
mpirun -N 3 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.3
mpirun -N 3 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.3

mpirun -N 4 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.4
mpirun -N 4 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.4
mpirun -N 4 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.4
mpirun -N 4 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.4
mpirun -N 4 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.4
mpirun -N 4 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.4
mpirun -N 4 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.4
mpirun -N 4 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.4
mpirun -N 4 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.4
mpirun -N 4 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.4
mpirun -N 4 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.4
mpirun -N 4 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.4
mpirun -N 4 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.4

mpirun -N 5 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.5
mpirun -N 5 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.5
mpirun -N 5 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.5
mpirun -N 5 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.5
mpirun -N 5 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.5
mpirun -N 5 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.5
mpirun -N 5 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.5
mpirun -N 5 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.5
mpirun -N 5 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.5
mpirun -N 5 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.5
mpirun -N 5 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.5
mpirun -N 5 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.5
mpirun -N 5 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.5

mpirun -N 6 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.6
mpirun -N 6 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.6
mpirun -N 6 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.6
mpirun -N 6 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.6
mpirun -N 6 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.6
mpirun -N 6 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.6
mpirun -N 6 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.6
mpirun -N 6 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.6
mpirun -N 6 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.6
mpirun -N 6 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.6
mpirun -N 6 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.6
mpirun -N 6 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.6
mpirun -N 6 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.6

mpirun -N 7 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.7
mpirun -N 7 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.7
mpirun -N 7 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.7
mpirun -N 7 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.7
mpirun -N 7 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.7
mpirun -N 7 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.7
mpirun -N 7 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.7
mpirun -N 7 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.7
mpirun -N 7 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.7
mpirun -N 7 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.7
mpirun -N 7 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.7
mpirun -N 7 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.7
mpirun -N 7 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.7

mpirun -N 8 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.8
mpirun -N 8 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.8
mpirun -N 8 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.8
mpirun -N 8 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.8
mpirun -N 8 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.8
mpirun -N 8 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.8
mpirun -N 8 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.8
mpirun -N 8 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.8
mpirun -N 8 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.8
mpirun -N 8 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.8
mpirun -N 8 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.8
mpirun -N 8 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.8
mpirun -N 8 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.8

mpirun -N 9 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.9
mpirun -N 9 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.9
mpirun -N 9 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.9
mpirun -N 9 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.9
mpirun -N 9 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.9
mpirun -N 9 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.9
mpirun -N 9 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.9
mpirun -N 9 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.9
mpirun -N 9 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.9
mpirun -N 9 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.9
mpirun -N 9 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.9
mpirun -N 9 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.9
mpirun -N 9 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.9

mpirun -N 10 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.10
mpirun -N 10 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.10
mpirun -N 10 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.10
mpirun -N 10 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.10
mpirun -N 10 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.10
mpirun -N 10 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.10
mpirun -N 10 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.10
mpirun -N 10 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.10
mpirun -N 10 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.10
mpirun -N 10 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.10
mpirun -N 10 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.10
mpirun -N 10 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.10
mpirun -N 10 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.10

mpirun -N 11 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.11
mpirun -N 11 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.11
mpirun -N 11 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.11
mpirun -N 11 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.11
mpirun -N 11 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.11
mpirun -N 11 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.11
mpirun -N 11 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.11
mpirun -N 11 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.11
mpirun -N 11 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.11
mpirun -N 11 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.11
mpirun -N 11 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.11
mpirun -N 11 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.11
mpirun -N 11 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.11

mpirun -N 12 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.12
mpirun -N 12 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.12
mpirun -N 12 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.12
mpirun -N 12 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.12
mpirun -N 12 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.12
mpirun -N 12 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.12
mpirun -N 12 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.12
mpirun -N 12 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.12
mpirun -N 12 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.12
mpirun -N 12 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.12
mpirun -N 12 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.12
mpirun -N 12 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.12
mpirun -N 12 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.12

mpirun -N 13 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.13
mpirun -N 13 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.13
mpirun -N 13 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.13
mpirun -N 13 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.13
mpirun -N 13 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.13
mpirun -N 13 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.13
mpirun -N 13 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.13
mpirun -N 13 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.13
mpirun -N 13 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.13
mpirun -N 13 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.13
mpirun -N 13 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.13
mpirun -N 13 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.13
mpirun -N 13 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.13

mpirun -N 14 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.14
mpirun -N 14 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.14
mpirun -N 14 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.14
mpirun -N 14 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.14
mpirun -N 14 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.14
mpirun -N 14 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.14
mpirun -N 14 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.14
mpirun -N 14 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.14
mpirun -N 14 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.14
mpirun -N 14 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.14
mpirun -N 14 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.14
mpirun -N 14 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.14
mpirun -N 14 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.14

mpirun -N 15 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.15
mpirun -N 15 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.15
mpirun -N 15 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.15
mpirun -N 15 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.15
mpirun -N 15 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.15
mpirun -N 15 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.15
mpirun -N 15 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.15
mpirun -N 15 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.15
mpirun -N 15 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.15
mpirun -N 15 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.15
mpirun -N 15 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.15
mpirun -N 15 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.15
mpirun -N 15 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.15

mpirun -N 16 ./lmp_cpu -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2K.16
mpirun -N 16 ./lmp_cpu -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4K.16
mpirun -N 16 ./lmp_cpu -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8K.16
mpirun -N 16 ./lmp_cpu -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.16K.16
mpirun -N 16 ./lmp_cpu -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.32K.16
mpirun -N 16 ./lmp_cpu -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.64K.16
mpirun -N 16 ./lmp_cpu -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.128K.16
mpirun -N 16 ./lmp_cpu -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.256K.16
mpirun -N 16 ./lmp_cpu -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.512K.16
mpirun -N 16 ./lmp_cpu -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.1024K.16
mpirun -N 16 ./lmp_cpu -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.2048K.16
mpirun -N 16 ./lmp_cpu -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.4096K.16
mpirun -N 16 ./lmp_cpu -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cpu
mv log.lammps log.28Jun14.lj.cpu.8192K.16
