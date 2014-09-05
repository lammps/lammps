#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.1
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 1 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.1

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.2
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 2 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.2

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.3
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 3 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.3

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.4
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 4 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.4

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.5
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 5 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.5

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.6
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 6 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.6

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.7
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 7 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.7

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.8
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 8 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.8

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.9
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 9 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.9

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.10
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 10 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.10

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.11
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 11 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.11

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.12
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 12 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.12

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.13
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 13 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.13

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.14
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 14 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.14

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.15
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 15 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.15

mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.1.16
mpirun -np 1 ./lmp_kokkos_cuda -k on g 1 t 16 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.1.16

mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.2.1
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 1 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.2.1

mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.2.2
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 2 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.2.2

mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.2.3
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 3 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.2.3

mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.2.4
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 4 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.2.4

mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.2.5
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 5 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.2.5

mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.2.6
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 6 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.2.6

mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.2.7
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 7 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.2.7

mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.16K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.32K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.64K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.128K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.256K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.512K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.1024K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.2048K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.4096K.2.8
mpirun -np 2 ./lmp_kokkos_cuda -k on g 2 t 8 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.cuda
mv log.lammps log.28Jun14.lj.kokkos.cuda.8192K.2.8
