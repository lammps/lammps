#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -np 1 lmp_intel_cpu -sf intel -pk intel 1 prec single -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.single.128K.1

mpirun -np 2 lmp_intel_cpu -sf intel -pk intel 1 prec single -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.single.128K.2

mpirun -np 4 lmp_intel_cpu -sf intel -pk intel 1 prec single -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.single.128K.4

mpirun -np 6 lmp_intel_cpu -sf intel -pk intel 1 prec single -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.single.128K.6

mpirun -np 8 lmp_intel_cpu -sf intel -pk intel 1 prec single -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.single.128K.8

mpirun -np 10 lmp_intel_cpu -sf intel -pk intel 1 prec single -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.single.128K.10

mpirun -np 12 lmp_intel_cpu -sf intel -pk intel 1 prec single -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.single.128K.12

mpirun -np 14 lmp_intel_cpu -sf intel -pk intel 1 prec single -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.single.128K.14

mpirun -np 16 lmp_intel_cpu -sf intel -pk intel 1 prec single -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.single.128K.16

mpirun -np 1 lmp_intel_cpu -sf intel -pk intel 1 prec mixed -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.mixed.128K.1

mpirun -np 2 lmp_intel_cpu -sf intel -pk intel 1 prec mixed -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.mixed.128K.2

mpirun -np 4 lmp_intel_cpu -sf intel -pk intel 1 prec mixed -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.mixed.128K.4

mpirun -np 6 lmp_intel_cpu -sf intel -pk intel 1 prec mixed -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.mixed.128K.6

mpirun -np 8 lmp_intel_cpu -sf intel -pk intel 1 prec mixed -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.mixed.128K.8

mpirun -np 10 lmp_intel_cpu -sf intel -pk intel 1 prec mixed -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.mixed.128K.10

mpirun -np 12 lmp_intel_cpu -sf intel -pk intel 1 prec mixed -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.mixed.128K.12

mpirun -np 14 lmp_intel_cpu -sf intel -pk intel 1 prec mixed -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.mixed.128K.14

mpirun -np 16 lmp_intel_cpu -sf intel -pk intel 1 prec mixed -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.mixed.128K.16

mpirun -np 1 lmp_intel_cpu -sf intel -pk intel 1 prec double -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.double.128K.1

mpirun -np 2 lmp_intel_cpu -sf intel -pk intel 1 prec double -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.double.128K.2

mpirun -np 4 lmp_intel_cpu -sf intel -pk intel 1 prec double -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.double.128K.4

mpirun -np 6 lmp_intel_cpu -sf intel -pk intel 1 prec double -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.double.128K.6

mpirun -np 8 lmp_intel_cpu -sf intel -pk intel 1 prec double -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.double.128K.8

mpirun -np 10 lmp_intel_cpu -sf intel -pk intel 1 prec double -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.double.128K.10

mpirun -np 12 lmp_intel_cpu -sf intel -pk intel 1 prec double -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.double.128K.12

mpirun -np 14 lmp_intel_cpu -sf intel -pk intel 1 prec double -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.double.128K.14

mpirun -np 16 lmp_intel_cpu -sf intel -pk intel 1 prec double -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.intel.cpu.double.128K.16
