#!/usr/bin/env python3

import numpy as np, pylab, tkinter
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from decimal import *
import sys, string, os


argv = sys.argv
if len(argv) != 3:
  print("Syntax: ./plot_precession.py res_nvt_spin.dat res_nvt_lattice.dat")
  sys.exit()

dirname   = os.path.join(os.getcwd(), "Feb_07")
nvtspin_file = sys.argv[1]
nvtlatt_file = sys.argv[2]

ts,tmags,temps = np.loadtxt(nvtspin_file,skiprows=0, usecols=(1,2,3),unpack=True)
tl,tmagl,templ = np.loadtxt(nvtlatt_file,skiprows=0, usecols=(1,2,3),unpack=True)

fig = plt.figure(figsize=(8,8))
ax1  = plt.subplot(2,1,1)
ax2  = plt.subplot(2,1,2)

ax1.plot(ts, tmags, 'r-', label='Spin temp. (thermostat)')
ax1.plot(ts, temps, 'g-', label='Lattice temp.')
ax1.set_yscale("log")
ax1.set_ylabel("T (K)")
ax1.legend(loc=3)

ax2.plot(tl, tmagl, 'r-', label='Spin temp.')
ax2.plot(tl, templ, 'g-', label='Lattice temp. (thermostat)')
ax2.set_yscale("log")
ax2.set_ylabel("T (K)")
ax2.legend(loc=3)

plt.xlabel('Time (in ps)')
plt.legend()
plt.show()

fig.savefig(os.path.join(os.getcwd(), "nve_spin_lattice.pdf"), bbox_inches="tight")
plt.close(fig)
