from ase import Atoms, Atom
from ase.calculators.lammpslib import LAMMPSlib
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

GaN = Atoms([Atom('Ga', (1.59, 0.917986928012, 0.0)),
             Atom('Ga', (1.59, -0.917986928012, 2.583)),
             Atom('N', (1.59, 0.917986928012, 1.98891)),
             Atom('N', (1.59, -0.917986928012, 4.57191))],
            cell=[(1.59, -2.75396078403, 0.0), (1.59, 2.75396078403, 0.0), (0.0, 0.0, 5.166)],
            pbc=True)

cmds = ["pair_style tersoff", "pair_coeff * * ../../../../../potentials/GaN.tersoff Ga N"]

mends = ["info system",
         "dynamical_matrix all eskm 0.000001 file dynmat.dat binary no",
         "neigh_modify delay 0"]

N = 6
GaN = GaN.repeat([N, N, N])

lammps = LAMMPSlib(lmpcmds=cmds, atom_types={'Ga': 1, 'N': 2}, amendments=mends, log_file='lammps.log')

GaN.set_calculator(lammps)
GaN.get_potential_energy()

if rank == 0:
    dynmat = np.loadtxt("dynmat.dat")
    dynmat = dynmat.reshape(([int(3*(len(dynmat)/3)**0.5), int(3*(len(dynmat)/3)**0.5)]))
    eigv = np.linalg.eigvals(dynmat)
    eigv.sort()
    eigv = np.sqrt(np.abs(eigv))/(2*np.pi)
    plt.hist(eigv, 80)
    plt.xlabel('Frequency (THz)')
    plt.show()
