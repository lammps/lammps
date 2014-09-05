#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -np 16 lmp_opt -sf opt -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.2K.16
mpirun -np 16 lmp_opt -sf opt -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.4K.16
mpirun -np 16 lmp_opt -sf opt -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.8K.16
mpirun -np 16 lmp_opt -sf opt -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.16K.16
mpirun -np 16 lmp_opt -sf opt -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.32K.16
mpirun -np 16 lmp_opt -sf opt -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.64K.16
mpirun -np 16 lmp_opt -sf opt -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.128K.16
mpirun -np 16 lmp_opt -sf opt -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.256K.16
mpirun -np 16 lmp_opt -sf opt -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.512K.16
mpirun -np 16 lmp_opt -sf opt -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.1024K.16
mpirun -np 16 lmp_opt -sf opt -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.2048K.16
mpirun -np 16 lmp_opt -sf opt -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.4096K.16
mpirun -np 16 lmp_opt -sf opt -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.opt
mv log.lammps log.28Jun14.lj.opt.8192K.16
