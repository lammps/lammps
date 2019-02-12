from ase import Atoms, Atom
from ase.calculators.lammpslib import LAMMPSlib
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

GaAs = Atoms([Atom('Ga', (0.0, 0.0, 0.0)),
              Atom('As', (1.413425, 1.413425, 1.413425))],
             cell=[(0.0, 2.82685, 2.82685), (2.82685, 0.0, 2.82685), (2.82685, 2.82685, 0.0)],
             pbc=True,)

cmds = ["pair_style bop", "pair_coeff * * ../../../../../potentials/GaAs.bop.table Ga As",
        "comm_modify cutoff 12"]

mends = ["info system",
         "dynamical_matrix all eskm 0.000001 file dynmat.dat binary no",
         "neigh_modify delay 0"]

N = 5
GaAs = GaAs.repeat([N, N, N])

lammps = LAMMPSlib(lmpcmds=cmds, atom_types={'Ga': 1, 'As': 2}, amendments=mends, log_file='lammps.log')

GaAs.set_calculator(lammps)
GaAs.get_potential_energy()

if rank == 0:
    dynmat = np.loadtxt("dynmat.dat")
    dynmat = dynmat.reshape(([int(3*(len(dynmat)/3)**0.5), int(3*(len(dynmat)/3)**0.5)]))
    eigv = np.linalg.eigvals(dynmat)
    eigv.sort()
    eigv = np.sqrt(np.abs(eigv))/(2*np.pi)
    plt.hist(eigv, 80)
    plt.xlabel('Frequency (THz)')
    plt.show()

