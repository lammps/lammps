#!/usr/bin/env python
# mbar.py - calculate free energy from compute mbar results using pymbar

import sys
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import pymbar

parser = ArgumentParser(description='Compute free energy using MBAR from alchemical simulation data.')
parser.add_argument("units", help="unit system can be lj, real or si")
parser.add_argument("temperature", type=float, help="The temperature of the system")
parser.add_argument("datafile", help="File with u_kln data in .npy format")
parser.add_argument("-g", "--grid", type=float, nargs='+', metavar='VALUE',
                    help="values of the perturbed parameter at each state, used as x axis "
                    "of the plot: either 'lo hi' for equally spaced values or one value "
                    "per state (default: state index)")

args = parser.parse_args()

# same unit systems and constants as fep.py
r_value = {'lj': 1.0, 'real': 0.0019872036, 'si': 8.31446}
e_units = {'lj': 'energy units', 'real': 'kcal/mol', 'si': 'J/mol'}

if args.units in r_value:
    kT = r_value[args.units] * args.temperature
else:
    sys.exit("The provided units keyword is not valid")

u_kln = np.load(args.datafile)  # shape (nstates, nstates, nsamples)
(nstates, _, _) = u_kln.shape

# x axis of the plot: the grid of the perturbed parameter, if provided
if args.grid is None:
    grid = np.arange(1, nstates + 1)
    xlabel = 'state'
elif len(args.grid) == 2:
    grid = np.linspace(args.grid[0], args.grid[1], nstates)
    xlabel = r'$\lambda$'
elif len(args.grid) == nstates:
    grid = np.array(args.grid)
    xlabel = r'$\lambda$'
else:
    sys.exit(f"The grid must have 2 values (lo hi) or one value per state ({nstates})")

## Subsample data to extract uncorrelated equilibrium timeseries
N_k = np.zeros([nstates], np.int32) # number of uncorrelated samples
for k in range(nstates):
    [nequil, g, Neff_max] = pymbar.timeseries.detect_equilibration(u_kln[k,k,:])
    # discard the initial non-equilibrium part, then subsample the remainder
    u_kln_equil = u_kln[k,:,nequil:]
    indices = pymbar.timeseries.subsample_correlated_data(u_kln[k,k,nequil:], g=g)
    N_k[k] = len(indices)
    u_kln[k,:,0:N_k[k]] = u_kln_equil[:,indices]

# Compute free energy differences
mbar = pymbar.MBAR(u_kln, N_k)

# If this fails try setting compute_uncertainty to false
# See this issue: https://github.com/choderalab/pymbar/issues/419
results = mbar.compute_free_energy_differences(compute_uncertainty=True)

deltaf = results['Delta_f'][0,nstates-1]
udeltaf = results['dDelta_f'][0,nstates-1]

print("Free energy change")
print(deltaf, "+/-", udeltaf, 'kT')
deltaf *= kT
udeltaf *= kT
print(deltaf, "+/-", udeltaf, e_units[args.units])

deltafs = np.array([ results['Delta_f'][0,k] - results['Delta_f'][0,k-1] for k in range(1,nstates) ])

fig, ax = plt.subplots()

ax.plot(grid[:nstates-1], deltafs * kT, marker='o')
ax.set(xlabel=xlabel, ylabel=rf'$\Delta G$ [{e_units[args.units]}]')
fig.savefig('deltaG_vs_lambda.png')
print('Plot saved to deltaG_vs_lambda.png')
