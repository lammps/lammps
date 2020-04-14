#!/bin/bash

DATE=14Apr20

# bfo
cd bfo/
../../../src/lmp_serial -in in.spin.bfo
cp log.lammps log.${DATE}.spin.bfo.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.bfo
cp log.lammps log.${DATE}.spin.bfo.g++.4
rm log.lammps log.cite dump*.lammpstrj
cd ..

# fcc cobalt
cd cobalt_fcc/
../../../src/lmp_serial -in in.spin.cobalt_fcc
cp log.lammps log.${DATE}.spin.cobalt_fcc.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.cobalt_fcc
cp log.lammps log.${DATE}.spin.cobalt_fcc.g++.4
rm log.lammps log.cite dump*.lammpstrj
cd ..

# hcp cobalt
cd cobalt_hcp/
../../../src/lmp_serial -in in.spin.cobalt_hcp
cp log.lammps log.${DATE}.spin.cobalt_hcp.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.cobalt_hcp
cp log.lammps log.${DATE}.spin.cobalt_hcp.g++.4
rm log.lammps log.cite dump*.lammpstrj
cd ..

# dipole spin
cd dipole_spin/
../../../src/lmp_serial -in in.spin.iron_dipole_cut
cp log.lammps log.${DATE}.spin.iron_dipole_cut.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.iron_dipole_cut
cp log.lammps log.${DATE}.spin.iron_dipole_cut.g++.4
../../../src/lmp_serial -in in.spin.iron_dipole_ewald
cp log.lammps log.${DATE}.spin.iron_dipole_ewald.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.iron_dipole_ewald
cp log.lammps log.${DATE}.spin.iron_dipole_ewald.g++.4
../../../src/lmp_serial -in in.spin.iron_dipole_pppm
cp log.lammps log.${DATE}.spin.iron_dipole_pppm.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.iron_dipole_pppm
cp log.lammps log.${DATE}.spin.iron_dipole_pppm.g++.4
rm log.lammps log.cite dump*.lammpstrj
cd ..

# bcc iron
cd iron/
../../../src/lmp_serial -in in.spin.iron
cp log.lammps log.${DATE}.spin.iron.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.iron
cp log.lammps log.${DATE}.spin.iron.g++.4
../../../src/lmp_serial -in in.spin.iron_cubic
cp log.lammps log.${DATE}.spin.iron_cubic.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.iron_cubic
cp log.lammps log.${DATE}.spin.iron_cubic.g++.4
rm log.lammps log.cite dump*.lammpstrj
cd ..

# fcc nickel
cd nickel/
../../../src/lmp_serial -in in.spin.nickel
cp log.lammps log.${DATE}.spin.nickel.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.nickel
cp log.lammps log.${DATE}.spin.nickel.g++.4
../../../src/lmp_serial -in in.spin.nickel_cubic
cp log.lammps log.${DATE}.spin.nickel_cubic.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.nickel_cubic
cp log.lammps log.${DATE}.spin.nickel_cubic.g++.4
rm log.lammps log.cite dump*.lammpstrj
cd ..

# read restart
cd read_restart/
../../../src/lmp_serial -in in.spin.write_restart
cp log.lammps log.${DATE}.spin.write_restart.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.write_restart
cp log.lammps log.${DATE}.spin.write_restart.g++.4
../../../src/lmp_serial -in in.spin.restart
cp log.lammps log.${DATE}.spin.restart.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.restart
cp log.lammps log.${DATE}.spin.restart.g++.4
../../../src/lmp_serial -in in.spin.read_data
cp log.lammps log.${DATE}.spin.read_data.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.read_data
cp log.lammps log.${DATE}.spin.read_data.g++.4
rm log.lammps log.cite dump*.lammpstrj
cd ..

# setforce
cd setforce_spin/
../../../src/lmp_serial -in in.spin.setforce
cp log.lammps log.${DATE}.spin.setforce.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.setforce
cp log.lammps log.${DATE}.spin.setforce.g++.4
rm log.lammps log.cite dump*.lammpstrj
cd ..

# spin minimizers
cd spinmin/
../../../src/lmp_serial -in in.spin.bfo_min
cp log.lammps log.${DATE}.spin.bfo_min.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.bfo_min
cp log.lammps log.${DATE}.spin.bfo_min.g++.4
../../../src/lmp_serial -in in.spin.bfo_min_cg
cp log.lammps log.${DATE}.spin.bfo_min_cg.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.bfo_min_cg
cp log.lammps log.${DATE}.spin.bfo_min_cg.g++.4
../../../src/lmp_serial -in in.spin.bfo_min_lbfgs
cp log.lammps log.${DATE}.spin.bfo_min_lbfgs.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.bfo_min_lbfgs
cp log.lammps log.${DATE}.spin.bfo_min_lbfgs.g++.4
../../../src/lmp_serial -in in.spin.iron_min
cp log.lammps log.${DATE}.spin.iron_min.g++.1
mpirun -np 4 ../../../src/lmp_mpi -in in.spin.iron_min
cp log.lammps log.${DATE}.spin.iron_min.g++.4
rm log.lammps log.cite dump*.lammpstrj
cd ..
