#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.2K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.4K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.8K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.16K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.32K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.64K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.128K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.256K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.512K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.1024K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.2048K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.4096K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 1 -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.single.8192K.16

mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.2K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.4K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.8K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.16K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.32K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.64K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.128K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.256K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.512K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.1024K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.2048K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.4096K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 3 -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.mixed.8192K.16

mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.2K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.4K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.8K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.16K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.32K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.64K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.128K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.256K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.512K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.1024K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.2048K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.4096K.16
mpirun -np 16 lmp_intel_cpu -sf intel -v a 2 -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.intel.cpu
mv log.lammps log.28Jun14.lj.intel.cpu.double.8192K.16
