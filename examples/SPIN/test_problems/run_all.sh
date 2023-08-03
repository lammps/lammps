#!/bin/bash

# test 1: damping and exchange 
cd validation_damped_exchange/
./run-test-exchange.sh
rm dump.data res_lammps.dat res_llg.dat
cd ..

# test 2: damping and Zeeman
cd validation_damped_precession/
./run-test-prec.sh
rm res_lammps.dat res_llg.dat
cd ..

# test 3: langevin, damping and Zeeman, low damping (0.01)
cd validation_langevin_precession_d0.01/
./run-test-prec.sh
rm average_spin test-prec-spin.in res_lammps.dat res_langevin.dat
cd ..

# test 4: langevin, damping and Zeeman, large damping (1.0)
cd validation_langevin_precession_d1.0/
./run-test-prec.sh
rm average_spin test-prec-spin.in res_lammps.dat res_langevin.dat
cd ..

# test 5: NVE run, test Etot preservation
cd validation_nve/
./run-test-nve.sh
rm nve_spin_lattice.pdf res_lammps.dat
cd ..
