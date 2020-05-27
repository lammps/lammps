"""Made by Charlie Sievers Ph.D. Candidate, UC Davis, Donadio Lab 2019"""
# from mpi4py import MPI
from lammps import lammps
import numpy as np

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()

""" LAMMPS  VARIABLES """

# data files
infile = "silicon_input_file.lmp"
ff_file = "ff-silicon.lmp"

# full output useful for testing
lmp = lammps()

# reduced output useful reducing IO for production runs
# lmp = lammps(cmdargs=["-screen", "none", "-log", "none"])

# lammps commands
lmp.command("atom_style full")
lmp.command("units metal")
lmp.command("processors * * *")
lmp.command("neighbor 1 bin")
lmp.command("boundary p p p")

# read data and force field file
lmp.command("read_data {}".format(infile))
lmp.file("{}".format(ff_file))

lmp.command("dynamical_matrix all eskm 0.000001 file dynmat.dat")

dynmat = np.loadtxt("dynmat.dat")
dynlen = int(3*np.sqrt(len(dynmat)/3))
dynmat = dynmat.reshape((dynlen, dynlen))

eigvals, eigvecs = np.linalg.eig(dynmat)

# frequencies in THz
omegas = np.sqrt(np.abs(eigvals))
print(omegas)
