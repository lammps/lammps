from ase import Atoms, Atom
from ase.calculators.lammpslib import LAMMPSlib
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

quartz = Atoms(
    [Atom('Si', (1.1545226, -1.99969180169, 0.0)),
     Atom('Si', (1.1545226, 1.99969180169, 3.6036)),
     Atom('Si', (2.6069548, 2.15247249027e-16, 1.8018)),
     Atom('O', (1.6724232, -0.624132037742, 0.64378314)),
     Atom('O', (1.6724232, 0.624132037742, 2.9598186618)),
     Atom('O', (2.1623026, -2.49695388906, 4.2473849418)),
     Atom('O', (3.5392742, 1.13629495821, 1.1580150582)),
     Atom('O', (3.5392742, -1.13629495821, 2.4455813382)),
     Atom('O', (2.1623026, 2.49695388906, 4.76161686))],
    cell=[(2.458, -4.257380885, 0.0), (2.458, 4.257380885, 0.0), (0.0, 0.0, 5.4054)],
    pbc=True,
    )

# number of repeats
N = 3
quartz = quartz.repeat([N, N, N])

header = ['units metal',
          'atom_style charge',
          'atom_modify map array sort 0 0']

cmds = ["pair_style      buck/coul/long 10.0 8.0",
        "pair_coeff	1 1      0      1         0",
        "pair_coeff	1 2   18003.7572 0.20520 133.5381",
        "pair_coeff	2 2    1388.7730 0.36232 175.0000",
        "kspace_style ewald 1.0e-12",
        "set type 1 charge 2.4",
        "set type 2 charge -1.2"]

mends = ["dynamical_matrix all eskm 0.000001 file dynmat.dat binary no",
         "neigh_modify delay 0"]


lammps = LAMMPSlib(lmpcmds=cmds, lammps_header=header, amendments=mends, log_file='lammps.log')

quartz.set_calculator(lammps)
quartz.get_potential_energy()

if rank == 0:
    dynmat = np.loadtxt("dynmat.dat")
    dynmat = dynmat.reshape(([int(3*(len(dynmat)/3)**0.5), int(3*(len(dynmat)/3)**0.5)]))
    eigv = np.linalg.eigvals(dynmat)
    eigv.sort()
    plt.hist(33*np.sqrt(np.abs(eigv))/(2*np.pi), 80)
    plt.xlabel('Frequency (cm-1)')
    plt.show()

