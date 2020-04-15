#!/usr/bin/env python3

import numpy as np, pylab, tkinter
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from decimal import *
import sys, string, os


argv = sys.argv
if len(argv) != 3:
  print("Syntax: ./plot_precession.py res_lammps.dat")
  sys.exit()

dirname   = os.path.join(os.getcwd(), "Feb_07")
lammps_file = sys.argv[1]
llg_file = sys.argv[2]

t,tmag,temp,e_mag,e_kin,e_pot,e_tot = np.loadtxt(lammps_file,skiprows=0, usecols=(1,2,3,4,5,6,7),unpack=True)

fig = plt.figure(figsize=(8,8))
ax1  = plt.subplot(3,1,1)
ax2  = plt.subplot(3,1,2)
ax3  = plt.subplot(3,1,3)

ax1.plot(t, e_tot, 'k-', label='Total energy')
ax1.plot(t, e_pot, 'r-', label='Potential energy')
ax1.set_ylabel("E (eV)")
ax1.legend(loc=3)

ax2.plot(t, e_kin, 'b-', label='Kinetic energy')
ax2.plot(t, e_mag, 'g-', label='Magnetic energy')
ax2.set_ylabel("E (eV)")
ax2.legend(loc=3)

ax3.plot(t, temp, 'b--', label='Latt. temperature')
ax3.plot(t, tmag, 'r--', label='Spin temperature')
ax3.set_ylabel("T (K)")
ax3.legend(loc=3)

plt.xlabel('Time (in ps)')
plt.legend()
plt.show()

fig.savefig(os.path.join(os.getcwd(), "nve_spin_lattice.pdf"), bbox_inches="tight")
plt.close(fig)
