#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.2K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.4K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.8K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.16K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.32K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.64K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.128K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.256K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.512K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.1024K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.2048K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.4096K.1.8
mpirun -np 1 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.8192K.1.8

mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.2K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.4K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.8K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.16K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.32K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.64K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.128K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.256K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.512K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.1024K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.2048K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.4096K.2.8
mpirun -np 2 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 8 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.8192K.2.8

mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.2K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.4K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.8K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.16K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.32K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.64K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.128K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.256K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.512K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.1024K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.2048K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.4096K.4.4
mpirun -np 4 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 4 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.8192K.4.4

mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.2K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.4K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.8K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.16K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.32K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.64K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.128K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.256K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.512K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.1024K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.2048K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.4096K.8.2
mpirun -np 8 -bind-to socket -map-by socket -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 2 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.omp
mv log.lammps log.28Jun14.lj.kokkos.omp.8192K.8.2

mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 8 -v y 8 -v z 8 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.2K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 8 -v y 8 -v z 16 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.4K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 8 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.8K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 16 -v y 16 -v z 16 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.16K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 16 -v y 16 -v z 32 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.32K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 16 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.64K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.128K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 32 -v y 32 -v z 64 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.256K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 32 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.512K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 64 -v y 64 -v z 64 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.1024K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 64 -v y 64 -v z 128 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.2048K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 64 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.4096K.16.1
mpirun -np 16 -bind-to core -map-by core -x KMP_AFFINITY=scatter ./lmp_kokkos_omp -k on t 1 -sf kk -v c device -v x 128 -v y 128 -v z 128 -v t 100 < in.lj.kokkos.omphalf
mv log.lammps log.28Jun14.lj.kokkos.omp.8192K.16.1
