
import sys
from lammps import lammps

lmp=lammps()
lmp=lammps(cmdargs=['-nocite','-sf','opt','-log','none'])

has_mpi=False
try:
    from mpi4py import MPI
    has_mpi=True
except:
    sys.exit(0)

mycomm=MPI.Comm.Split(MPI.COMM_WORLD, 0, 1)
lmp=lammps(comm=mycomm)

