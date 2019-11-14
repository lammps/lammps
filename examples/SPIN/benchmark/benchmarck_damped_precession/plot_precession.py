#!/usr/bin/env python3

import numpy as np, pylab, tkinter
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from decimal import *
import sys, string, os


argv = sys.argv
if len(argv) != 3:
  print("Syntax: ./plot_precession.py res_lammps.dat res_llg.dat")
  sys.exit()

lammps_file = sys.argv[1]
llg_file = sys.argv[2]

t_lmp,Sx_lmp,Sy_lmp,Sz_lmp = np.loadtxt(lammps_file, skiprows=0, usecols=(1,2,3,4),unpack=True)
t_llg,Sx_llg,Sy_llg,Sz_llg = np.loadtxt(llg_file, skiprows=0, usecols=(0,1,2,3),unpack=True)

plt.figure()
plt.subplot(311)
plt.ylabel('Sx')
plt.plot(t_lmp, Sx_lmp, 'b-', label='LAMMPS')
plt.plot(t_llg, Sx_llg, 'r--', label='LLG')

plt.subplot(312)
plt.ylabel('Sy')
plt.plot(t_lmp, Sy_lmp, 'b-', label='LAMMPS')
plt.plot(t_llg, Sy_llg, 'r--', label='LLG')

plt.subplot(313)
plt.ylabel('Sz')
plt.plot(t_lmp, Sz_lmp, 'b-', label='LAMMPS')
plt.plot(t_llg, Sz_llg, 'r--', label='LLG')

plt.xlabel('time (in ps)')
plt.legend()
plt.show()
