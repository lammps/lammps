#!/usr/bin/env python3

#Program fitting the exchange interaction
#Model curve: Bethe-Slater function
import numpy as np, pylab, tkinter
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from decimal import *
import sys, string, os


argv = sys.argv
if len(argv) != 3:
  print("Syntax: ./plot_precession.py res_lammps.dat res_langevin.dat")
  sys.exit()

lammps_file = sys.argv[1]
langevin_file = sys.argv[2]

T_lmp,S_lmp,E_lmp = np.loadtxt(lammps_file, skiprows=0, usecols=(0,2,3),unpack=True)
T_lan,S_lan,E_lan = np.loadtxt(langevin_file, skiprows=0, usecols=(0,1,2),unpack=True)

plt.figure()
plt.subplot(211)
plt.ylabel('<Sz>')
plt.plot(T_lmp, S_lmp, 'b-', label='LAMMPS')
plt.plot(T_lan, S_lan, 'r--', label='Langevin')

plt.subplot(212)
plt.ylabel('E (in eV)')
plt.plot(T_lmp, E_lmp, 'b-', label='LAMMPS')
plt.plot(T_lan, E_lan, 'r--', label='Langevin')

plt.xlabel('T (in K)')
pylab.xlim([0,20.0])
plt.legend()
plt.show()
