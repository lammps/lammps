#!/bin/sh
# run test problems.
EXE=${1-../../src/lmp_mpi}
DATE=$(date +%1d%b%y)

# dimer

${EXE} -in in.water_dimer.amoeba -log log.${DATE}.water_dimer.amoeba.g++.1
mpirun -np 4 ${EXE} -in in.water_dimer.amoeba -log log.${DATE}.water_dimer.amoeba.g++.4

${EXE} -in in.water_dimer.hippo -log log.${DATE}.water_dimer.hippo.g++.1
mpirun -np 4 ${EXE} -in in.water_dimer.hippo -log log.${DATE}.water_dimer.hippo.g++.4

# hexamer

${EXE} -in in.water_hexamer.amoeba -log log.${DATE}.water_hexamer.amoeba.g++.1
mpirun -np 4 ${EXE} -in in.water_hexamer.amoeba -log log.${DATE}.water_hexamer.amoeba.g++.4

${EXE} -in in.water_hexamer.hippo -log log.${DATE}.water_hexamer.hippo.g++.1
mpirun -np 4 ${EXE} -in in.water_hexamer.hippo -log log.${DATE}.water_hexamer.hippo.g++.4

# water box

${EXE} -in in.water_box.amoeba -log log.${DATE}.water_box.amoeba.g++.1
mpirun -np 4 ${EXE} -in in.water_box.amoeba -log log.${DATE}.water_box.amoeba.g++.4

${EXE} -in in.water_box.hippo -log log.${DATE}.water_box.hippo.g++.1
mpirun -np 4 ${EXE} -in in.water_box.hippo -log log.${DATE}.water_box.hippo.g++.4

# ubiquitin

${EXE} -in in.ubiquitin -log log.${DATE}.ubiquitin.g++.1
mpirun -np 4 ${EXE} -in in.ubiquitin -log log.${DATE}.ubiquitin.g++.4

