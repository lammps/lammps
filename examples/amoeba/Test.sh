# run test problems

# dimer

../../src/lmp_mpi < in.water_dimer.amoeba
mv log.lammps log.water_dimer.amoeba.1.test
mv dump.water_dimer dump.water_dimer.amoeba.1.test

mpirun -np 4 ../../src/lmp_mpi < in.water_dimer.amoeba
mv log.lammps log.water_dimer.amoeba.4.test
mv dump.water_dimer dump.water_dimer.amoeba.4.test

../../src/lmp_mpi < in.water_dimer.hippo
mv log.lammps log.water_dimer.hippo.1.test
mv dump.water_dimer dump.water_dimer.hippo.1.test

mpirun -np 4 ../../src/lmp_mpi < in.water_dimer.hippo
mv log.lammps log.water_dimer.hippo.4.test
mv dump.water_dimer dump.water_dimer.hippo.4.test

# hexamer

../../src/lmp_mpi < in.water_hexamer.amoeba
mv log.lammps log.water_hexamer.amoeba.1.test
mv dump.water_hexamer dump.water_hexamer.amoeba.1.test

mpirun -np 4 ../../src/lmp_mpi < in.water_hexamer.amoeba
mv log.lammps log.water_hexamer.amoeba.4.test
mv dump.water_hexamer dump.water_hexamer.amoeba.4.test

../../src/lmp_mpi < in.water_hexamer.hippo
mv log.lammps log.water_hexamer.hippo.1.test
mv dump.water_hexamer dump.water_hexamer.hippo.1.test

mpirun -np 4 ../../src/lmp_mpi < in.water_hexamer.hippo
mv log.lammps log.water_hexamer.hippo.4.test
mv dump.water_hexamer dump.water_hexamer.hippo.4.test

# water box

../../src/lmp_mpi < in.water_box.amoeba
mv log.lammps log.water_box.amoeba.1.test
mv dump.water_box dump.water_box.amoeba.1.test

mpirun -np 32 ../../src/lmp_mpi < in.water_box.amoeba
mv log.lammps log.water_box.amoeba.32.test
mv dump.water_box dump.water_box.amoeba.32.test

../../src/lmp_mpi < in.water_box.hippo
mv log.lammps log.water_box.hippo.1.test
mv dump.water_box dump.water_box.hippo.1.test

mpirun -np 32 ../../src/lmp_mpi < in.water_box.hippo
mv log.lammps log.water_box.hippo.32.test
mv dump.water_box dump.water_box.hippo.32.test

# ubiquitin

../../src/lmp_mpi < in.ubiquitin
mv log.lammps log.ubi.1.test
mv dump.ubi dump.ubi.1.test

mpirun -np 32 ../../src/lmp_mpi < in.ubiquitin
mv log.lammps log.ubi.32.test
mv dump.ubi dump.ubi.32.test

