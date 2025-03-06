from mpi4py import MPI
from lammps import lammps

L = lammps()
L.file('in.melt')


if MPI.COMM_WORLD.rank == 0:
    pe = L.get_thermo("pe")
    print("Potential Energy:", pe)
