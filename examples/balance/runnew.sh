# ------------
# existing scripts

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance
mv log.lammps log.new.balance.g++.4

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance.bond.fast
mv log.lammps log.new.balance.bond.fast.g++.4

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance.bond.slow
mv log.lammps log.new.balance.bond.slow.g++.4

# ------------
# new scripts

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance.clock.dynamic
mv log.lammps log.new.balance.clock.dynamic.g++.4

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance.clock.static
mv log.lammps log.new.balance.clock.static.g++.4

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance.group.dynamic
mv log.lammps log.new.balance.group.dynamic.g++.4

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance.group.static
mv log.lammps log.new.balance.group.static.g++.4

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance.kspace
mv log.lammps log.new.balance.kspace.g++.4

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance.neigh.dynamic
mv log.lammps log.new.balance.neigh.dynamic.g++.4

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance.neigh.rcb
mv log.lammps log.new.balance.neigh.rcb.g++.4

mpirun -np 4 /home/sjplimp/lammps/src/lmp_g++ < in.balance.neigh.static
mv log.lammps log.new.balance.neigh.static.g++.4

mpirun -np 2 /home/sjplimp/lammps/src/lmp_g++ < in.balance.var.dynamic
mv log.lammps log.new.balance.var.dynamic.g++.2

