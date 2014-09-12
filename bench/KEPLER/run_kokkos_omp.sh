#!/bin/bash
#SBATCH -N 1 --time=12:00:00

mpirun -np full -bind-to socket -map-by socket -x KMP_AFFINITY=scatter lmp_kokkos_omp -k on t 16 -sf kk -pk kokkos neigh full newton off comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.omp.128K.1.16

mpirun -np full -bind-to socket -map-by socket -x KMP_AFFINITY=scatter lmp_kokkos_omp -k on t 8 -sf kk -pk kokkos neigh full newton off comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.omp.128K.2.8

mpirun -np full -bind-to socket -map-by socket -x KMP_AFFINITY=scatter lmp_kokkos_omp -k on t 4 -sf kk -pk kokkos neigh full newton off comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.omp.128K.4.4

mpirun -np full -bind-to socket -map-by socket -x KMP_AFFINITY=scatter lmp_kokkos_omp -k on t 2 -sf kk -pk kokkos neigh full newton off comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.omp.128K.8.2

mpirun -np half -bind-to socket -map-by socket -x KMP_AFFINITY=scatter lmp_kokkos_omp -k on t 1 -sf kk -pk kokkos neigh half newton on comm device -v x 32 -v y 32 -v z 32 -v t 100 < in.lj
mv log.lammps log.10Sep14.lj.kokkos.omp.128K.16.1
