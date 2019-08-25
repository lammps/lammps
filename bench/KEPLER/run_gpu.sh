#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -np 1 lmp_gpu_single -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.1.1

mpirun -np 2 lmp_gpu_single -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.2.1

mpirun -np 2 lmp_gpu_single -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.2.2

mpirun -np 4 lmp_gpu_single -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.4.1

mpirun -np 4 lmp_gpu_single -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.4.2

mpirun -np 6 lmp_gpu_single -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.6.1

mpirun -np 6 lmp_gpu_single -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.6.2

mpirun -np 8 lmp_gpu_single -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.8.1

mpirun -np 8 lmp_gpu_single -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.8.2

mpirun -np 10 lmp_gpu_single -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.10.1

mpirun -np 10 lmp_gpu_single -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.10.2

mpirun -np 12 lmp_gpu_single -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.12.1

mpirun -np 12 lmp_gpu_single -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.12.2

mpirun -np 14 lmp_gpu_single -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.14.1

mpirun -np 14 lmp_gpu_single -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.14.2

mpirun -np 16 lmp_gpu_single -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.16.1

mpirun -np 16 lmp_gpu_single -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.single.128K.16.2

mpirun -np 1 lmp_gpu_mixed -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.1.1

mpirun -np 2 lmp_gpu_mixed -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.2.1

mpirun -np 2 lmp_gpu_mixed -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.2.2

mpirun -np 4 lmp_gpu_mixed -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.4.1

mpirun -np 4 lmp_gpu_mixed -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.4.2

mpirun -np 6 lmp_gpu_mixed -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.6.1

mpirun -np 6 lmp_gpu_mixed -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.6.2

mpirun -np 8 lmp_gpu_mixed -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.8.1

mpirun -np 8 lmp_gpu_mixed -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.8.2

mpirun -np 10 lmp_gpu_mixed -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.10.1

mpirun -np 10 lmp_gpu_mixed -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.10.2

mpirun -np 12 lmp_gpu_mixed -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.12.1

mpirun -np 12 lmp_gpu_mixed -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.12.2

mpirun -np 14 lmp_gpu_mixed -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.14.1

mpirun -np 14 lmp_gpu_mixed -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.14.2

mpirun -np 16 lmp_gpu_mixed -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.16.1

mpirun -np 16 lmp_gpu_mixed -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.mixed.128K.16.2

mpirun -np 1 lmp_gpu_double -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.1.1

mpirun -np 2 lmp_gpu_double -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.2.1

mpirun -np 2 lmp_gpu_double -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.2.2

mpirun -np 4 lmp_gpu_double -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.4.1

mpirun -np 4 lmp_gpu_double -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.4.2

mpirun -np 6 lmp_gpu_double -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.6.1

mpirun -np 6 lmp_gpu_double -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.6.2

mpirun -np 8 lmp_gpu_double -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.8.1

mpirun -np 8 lmp_gpu_double -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.8.2

mpirun -np 10 lmp_gpu_double -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.10.1

mpirun -np 10 lmp_gpu_double -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.10.2

mpirun -np 12 lmp_gpu_double -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.12.1

mpirun -np 12 lmp_gpu_double -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.12.2

mpirun -np 14 lmp_gpu_double -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.14.1

mpirun -np 14 lmp_gpu_double -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.14.2

mpirun -np 16 lmp_gpu_double -sf gpu -pk gpu 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.16.1

mpirun -np 16 lmp_gpu_double -sf gpu -pk gpu 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.gpu.double.128K.16.2
