#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.2K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.4K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.8K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.16K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.32K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.64K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.128K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.256K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.512K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.1024K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.2048K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.4096K.1
mpirun -N 1 ./lmp_cuda_double -c on -sf cuda -v g 1 -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.8192K.1

mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.2K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.4K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.8K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.16K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.32K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.64K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.128K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.256K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.512K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.1024K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.2048K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.4096K.2
mpirun -N 2 ./lmp_cuda_double -c on -sf cuda -v g 2 -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.double.8192K.2

mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.2K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.4K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.8K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.16K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.32K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.64K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.128K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.256K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.512K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.1024K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.2048K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.4096K.1
mpirun -N 1 ./lmp_cuda_mixed -c on -sf cuda -v g 1 -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.8192K.1

mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.2K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.4K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.8K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.16K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.32K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.64K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.128K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.256K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.512K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.1024K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.2048K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.4096K.2
mpirun -N 2 ./lmp_cuda_mixed -c on -sf cuda -v g 2 -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.mixed.8192K.2

mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.2K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.4K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.8K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.16K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.32K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.64K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.128K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.256K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.512K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.1024K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.2048K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.4096K.1
mpirun -N 1 ./lmp_cuda_single -c on -sf cuda -v g 1 -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.8192K.1

mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.2K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.4K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.8K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.16K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.32K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.64K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.128K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.256K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.512K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.1024K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.2048K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.4096K.2
mpirun -N 2 ./lmp_cuda_single -c on -sf cuda -v g 2 -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.cuda
mv log.lammps log.28Jun14.lj.cuda.single.8192K.2
