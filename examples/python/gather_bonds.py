from lammps import lammps
import numpy as np
import ctypes as ct
lmp = lammps()

lmp.command('atom_style full')
lmp.command('read_data data.spce')

natoms = lmp.get_natoms()
nbonds = lmp.get_nbonds()
print("atom count: ", natoms)
print("bond count: ", nbonds)

bonds = np.array(lmp.gather_bonds(), dtype=ct.c_int).reshape(nbonds, 3)

print("bond list (type,id1,id2): ")
print(bonds)
