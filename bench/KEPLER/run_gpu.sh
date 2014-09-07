#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.1
mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.1
mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.1
mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.1
mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.1
mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.1
mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.1
mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.1
mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.1
mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.1
mpirun -N 1 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.1

mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.2

mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.3
mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.3
mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.3
mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.3
mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.3
mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.3
mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.3
mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.3
mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.3
mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.3
mpirun -N 3 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.3

mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.4

mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.5
mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.5
mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.5
mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.5
mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.5
mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.5
mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.5
mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.5
mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.5
mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.5
mpirun -N 5 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.5

mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.6

mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.7
mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.7
mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.7
mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.7
mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.7
mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.7
mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.7
mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.7
mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.7
mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.7
mpirun -N 7 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.7

mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.8

mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.9
mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.9
mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.9
mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.9
mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.9
mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.9
mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.9
mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.9
mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.9
mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.9
mpirun -N 9 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.9

mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.10

mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.11
mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.11
mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.11
mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.11
mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.11
mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.11
mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.11
mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.11
mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.11
mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.11
mpirun -N 11 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.11

mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.12

mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.13
mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.13
mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.13
mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.13
mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.13
mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.13
mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.13
mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.13
mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.13
mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.13
mpirun -N 13 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.13

mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.14

mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.15
mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.15
mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.15
mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.15
mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.15
mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.15
mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.15
mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.15
mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.15
mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.15
mpirun -N 15 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.15

mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.1.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.1.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.1.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.1.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.1.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.1.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.1.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.1.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.1.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.1.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.1.16

mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.2.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.2.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.2.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.2.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.2.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.2.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.2.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.2.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.2.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.2.2
mpirun -N 2 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.2.2

mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.2.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.2.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.2.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.2.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.2.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.2.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.2.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.2.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.2.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.2.4
mpirun -N 4 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.2.4

mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.2.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.2.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.2.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.2.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.2.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.2.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.2.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.2.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.2.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.2.6
mpirun -N 6 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.2.6

mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.2.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.2.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.2.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.2.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.2.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.2.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.2.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.2.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.2.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.2.8
mpirun -N 8 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.2.8

mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.2.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.2.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.2.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.2.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.2.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.2.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.2.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.2.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.2.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.2.10
mpirun -N 10 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.2.10

mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.2.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.2.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.2.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.2.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.2.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.2.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.2.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.2.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.2.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.2.12
mpirun -N 12 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.2.12

mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.2.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.2.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.2.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.2.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.2.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.2.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.2.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.2.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.2.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.2.14
mpirun -N 14 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.2.14

mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2K.2.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.4K.2.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.8K.2.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.16K.2.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.32K.2.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.64K.2.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.128K.2.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.256K.2.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.512K.2.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.1024K.2.16
mpirun -N 16 ./lmp_gpu_double -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.double.2048K.2.16

mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.1
mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.1
mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.1
mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.1
mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.1
mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.1
mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.1
mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.1
mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.1
mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.1
mpirun -N 1 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.1

mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.2

mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.3
mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.3
mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.3
mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.3
mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.3
mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.3
mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.3
mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.3
mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.3
mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.3
mpirun -N 3 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.3

mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.4

mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.5
mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.5
mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.5
mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.5
mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.5
mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.5
mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.5
mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.5
mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.5
mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.5
mpirun -N 5 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.5

mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.6

mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.7
mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.7
mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.7
mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.7
mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.7
mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.7
mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.7
mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.7
mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.7
mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.7
mpirun -N 7 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.7

mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.8

mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.9
mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.9
mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.9
mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.9
mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.9
mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.9
mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.9
mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.9
mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.9
mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.9
mpirun -N 9 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.9

mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.10

mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.11
mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.11
mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.11
mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.11
mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.11
mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.11
mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.11
mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.11
mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.11
mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.11
mpirun -N 11 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.11

mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.12

mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.13
mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.13
mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.13
mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.13
mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.13
mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.13
mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.13
mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.13
mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.13
mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.13
mpirun -N 13 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.13

mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.14

mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.15
mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.15
mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.15
mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.15
mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.15
mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.15
mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.15
mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.15
mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.15
mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.15
mpirun -N 15 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.15

mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.1.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.1.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.1.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.1.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.1.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.1.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.1.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.1.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.1.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.1.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.1.16

mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.2.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.2.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.2.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.2.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.2.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.2.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.2.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.2.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.2.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.2.2
mpirun -N 2 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.2.2

mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.2.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.2.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.2.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.2.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.2.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.2.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.2.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.2.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.2.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.2.4
mpirun -N 4 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.2.4

mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.2.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.2.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.2.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.2.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.2.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.2.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.2.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.2.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.2.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.2.6
mpirun -N 6 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.2.6

mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.2.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.2.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.2.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.2.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.2.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.2.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.2.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.2.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.2.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.2.8
mpirun -N 8 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.2.8

mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.2.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.2.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.2.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.2.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.2.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.2.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.2.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.2.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.2.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.2.10
mpirun -N 10 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.2.10

mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.2.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.2.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.2.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.2.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.2.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.2.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.2.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.2.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.2.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.2.12
mpirun -N 12 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.2.12

mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.2.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.2.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.2.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.2.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.2.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.2.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.2.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.2.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.2.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.2.14
mpirun -N 14 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.2.14

mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2K.2.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.4K.2.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.8K.2.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.16K.2.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.32K.2.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.64K.2.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.128K.2.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.256K.2.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.512K.2.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.1024K.2.16
mpirun -N 16 ./lmp_gpu_mixed -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.mixed.2048K.2.16

mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.1
mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.1
mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.1
mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.1
mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.1
mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.1
mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.1
mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.1
mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.1
mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.1
mpirun -N 1 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.1

mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.2

mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.3
mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.3
mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.3
mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.3
mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.3
mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.3
mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.3
mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.3
mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.3
mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.3
mpirun -N 3 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.3

mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.4

mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.5
mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.5
mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.5
mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.5
mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.5
mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.5
mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.5
mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.5
mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.5
mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.5
mpirun -N 5 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.5

mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.6

mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.7
mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.7
mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.7
mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.7
mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.7
mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.7
mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.7
mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.7
mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.7
mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.7
mpirun -N 7 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.7

mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.8

mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.9
mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.9
mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.9
mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.9
mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.9
mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.9
mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.9
mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.9
mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.9
mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.9
mpirun -N 9 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.9

mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.10

mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.11
mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.11
mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.11
mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.11
mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.11
mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.11
mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.11
mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.11
mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.11
mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.11
mpirun -N 11 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.11

mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.12

mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.13
mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.13
mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.13
mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.13
mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.13
mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.13
mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.13
mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.13
mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.13
mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.13
mpirun -N 13 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.13

mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.14

mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.15
mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.15
mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.15
mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.15
mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.15
mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.15
mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.15
mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.15
mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.15
mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.15
mpirun -N 15 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.15

mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.1.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.1.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.1.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.1.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.1.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.1.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.1.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.1.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.1.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.1.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.1.16

mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.2.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.2.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.2.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.2.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.2.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.2.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.2.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.2.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.2.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.2.2
mpirun -N 2 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.2.2

mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.2.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.2.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.2.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.2.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.2.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.2.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.2.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.2.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.2.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.2.4
mpirun -N 4 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.2.4

mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.2.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.2.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.2.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.2.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.2.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.2.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.2.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.2.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.2.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.2.6
mpirun -N 6 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.2.6

mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.2.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.2.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.2.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.2.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.2.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.2.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.2.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.2.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.2.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.2.8
mpirun -N 8 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.2.8

mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.2.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.2.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.2.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.2.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.2.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.2.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.2.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.2.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.2.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.2.10
mpirun -N 10 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.2.10

mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.2.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.2.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.2.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.2.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.2.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.2.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.2.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.2.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.2.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.2.12
mpirun -N 12 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.2.12

mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.2.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.2.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.2.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.2.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.2.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.2.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.2.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.2.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.2.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.2.14
mpirun -N 14 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.2.14

mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2K.2.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.4K.2.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.8K.2.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.16K.2.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.32K.2.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.64K.2.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.128K.2.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.256K.2.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.512K.2.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.1024K.2.16
mpirun -N 16 ./lmp_gpu_single -sf gpu -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.gpu
mv log.lammps log.28Jun14.lj.gpu.single.2048K.2.16
