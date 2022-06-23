"""
compute_snap_dgrad.py
Purpose: Demonstrate extraction of descriptor gradient (dB/dR) array from compute snap.
         Show that dBi/dRj components summed over neighbors i yields same output as regular compute snap with dgradflag=0.
         This shows that the dBi/dRj components extracted with dgradflag=1 are correct.
Serial syntax:
    python compute_snap_dgrad.py
Parallel syntax:
    mpirun -np 2 python compute_snap_dgrad.py
"""

from __future__ import print_function
import sys
import ctypes
import numpy as np
from lammps import lammps, LMP_TYPE_ARRAY, LMP_STYLE_GLOBAL

# get MPI settings from LAMMPS

lmp = lammps()
me = lmp.extract_setting("world_rank")
nprocs = lmp.extract_setting("world_size")

cmds = ["-screen", "none", "-log", "none"]
lmp = lammps(cmdargs=cmds)

def run_lammps(dgradflag):
    lmp.command("clear")
    lmp.command("units metal")
    lmp.command("boundary	p p p")
    lmp.command("atom_modify	map hash")
    lmp.command(f"lattice         bcc {latparam}")
    lmp.command(f"region		box block 0 {nx} 0 {ny} 0 {nz}")
    lmp.command(f"create_box	{ntypes} box")
    lmp.command(f"create_atoms	{ntypes} box")
    lmp.command("mass 		* 180.88")
    lmp.command("displace_atoms 	all random 0.01 0.01 0.01 123456")
    # Pair style
    snap_options=f'{rcutfac} {rfac0} {twojmax} {radelem1} {radelem2} {wj1} {wj2} rmin0 {rmin0} quadraticflag {quadratic} bzeroflag {bzero} switchflag {switch} bikflag {bikflag} dgradflag {dgradflag}'
    lmp.command(f"pair_style 	zero {rcutfac}")
    lmp.command(f"pair_coeff 	* *")
    lmp.command(f"pair_style 	zbl {zblcutinner} {zblcutouter}")
    lmp.command(f"pair_coeff 	* * {zblz} {zblz}")
    # set up compute snap generating global array
    lmp.command(f"compute 	snap all snap {snap_options}")
    # Run
    lmp.command(f"thermo 		100")
    lmp.command(f"run {nsteps}")

# Declare simulation/structure variables
nsteps=0
nrep=2
latparam=2.0
ntypes=2
nx=nrep
ny=nrep
nz=nrep

# Declare compute snap variables
twojmax=8
m = (twojmax/2)+1
rcutfac=1.0
rfac0=0.99363
rmin0=0
radelem1=2.3
radelem2=2.0
wj1=1.0
wj2=0.96
quadratic=0
bzero=0
switch=0
bikflag=1
dgradflag=1

# Declare reference potential variables
zblcutinner=4.0
zblcutouter=4.8
zblz=73

# Number of descriptors
if (twojmax % 2 == 0):
    nd = int(m*(m+1)*(2*m+1)/6)
else:
    nd = int(m*(m+1)*(m+2)/3)
if me == 0:
    print(f"Number of descriptors based on twojmax : {nd}")

# Run lammps with dgradflag on
if me == 0:
    print("Running with dgradflag on")
run_lammps(1)

# Get global snap array
lmp_snap = lmp.numpy.extract_compute("snap",0, 2)

# Extract dBj/dRi (includes dBi/dRi)
natoms = lmp.get_natoms()
fref1 = lmp_snap[0:natoms,0:3].flatten()
eref1 = lmp_snap[-1,0] #lmp_snap[-6,0]
dbdr_length = np.shape(lmp_snap)[0]-(natoms) - 1 # Length of neighborlist pruned dbdr array
dBdR = lmp_snap[natoms:(natoms+dbdr_length),3:(nd+3)]
force_indices = lmp_snap[natoms:(natoms+dbdr_length),0:3].astype(np.int32)

# Sum over neighbors j for each atom i, like dgradflag=0 does.
array1 = np.zeros((3*natoms,nd))
a = 0
for k in range(0,nd):
    for l in range(0,dbdr_length):
        j = force_indices[l,0]
        i = force_indices[l,1]
        array1[3*(i)+a,k] += dBdR[l,k]
        a = a+1
        if (a>2):
            a=0

# Run lammps with dgradflag off
print("Running with dgradflag off")
run_lammps(0)

# Get global snap array
lmp_snap = lmp.numpy.extract_compute("snap",0, 2)
natoms = lmp.get_natoms()
fref2 = lmp_snap[natoms:(natoms+3*natoms),-1]
eref2 = lmp_snap[0,-1]
array2 = lmp_snap[natoms:natoms+(3*natoms), nd:-1]

# Sum the arrays obtained from dgradflag on and off.
summ = array1 + array2
#np.savetxt("sum.dat", summ)

print(f"Maximum difference in descriptor sums: {np.max(summ)}")
print(f"Maximum difference in reference forces: {np.max(fref1-fref2)}")
print(f"Difference in reference energy: {np.max(eref1-eref2)}")
