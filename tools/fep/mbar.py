#!/usr/bin/env python

from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import pymbar

RGAS = 0.00831446261815324

parser = ArgumentParser(description='Compute free energy using MBAR from alchemical simulation data.')
parser.add_argument("temperature", type=float, help="The temperature of the system [K]")
parser.add_argument("datafile", help="File with u_kln data in .npy format")

args = parser.parse_args()
temperature = float(args.temperature)

kT = RGAS * temperature   # kJ/mol

u_kln = np.load(args.datafile)  # shape (nstates, nstates, nsamples)
(nstates, _, _) = u_kln.shape

## Subsample data to extract uncorrelated equilibrium timeseries
N_k = np.zeros([nstates], np.int32) # number of uncorrelated samples
for k in range(nstates):
    [nequil, g, Neff_max] = pymbar.timeseries.detect_equilibration(u_kln[k,k,:])
    indices = pymbar.timeseries.subsample_correlated_data(u_kln[k,k,:], g=g)
    N_k[k] = len(indices)
    u_kln[k,:,0:N_k[k]] = u_kln[k,:,indices].T

# Compute free energy differences
mbar = pymbar.MBAR(u_kln, N_k)

# If this fails try setting compute_uncertainty to false
# See this issue: https://github.com/choderalab/pymbar/issues/419
results = mbar.compute_free_energy_differences(compute_uncertainty=True)

deltaf = results['Delta_f'][0,nstates-1]
udeltaf = results['dDelta_f'][0,nstates-1]

print("Free energy change")
print(deltaf, "+/-", udeltaf, ' [kT]')
deltaf *= kT
udeltaf *= kT
print(deltaf, "+/-", udeltaf, ' kJ/mol')

deltafs = np.array([ results['Delta_f'][0,k] - results['Delta_f'][0,k-1] for k in range(1,nstates) ])

lambdas = np.linspace(0.0, 1.0, nstates)

fig, ax = plt.subplots()

ax.plot(lambdas[:nstates-1], deltafs * kT, marker='o')
ax.set(xlabel='λ', ylabel='ΔG [kJ/mol]')
fig.savefig('deltaG_vs_lambda.png')
print('Plot saved to deltaG_vs_lambda.png')
