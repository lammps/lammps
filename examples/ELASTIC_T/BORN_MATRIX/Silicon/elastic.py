#!/usr/bin/env python -i
# preceding line should have path for Python on your machine

# elastic.py
# Purpose: demonstrate elastic constant calculation for
# two different crystal supercells, one with non-standard orientation
#
# Syntax:  elastic.py
#          uses in.elastic as LAMMPS input script

from __future__ import print_function
from elastic_utils import *

np.set_printoptions(precision = 3, suppress=True)

# get MPI settings from LAMMPS

lmp = lammps()
me = lmp.extract_setting("world_rank")
nprocs = lmp.extract_setting("world_size")

# cubic diamond lattice constants 

alat = 5.457

# define the cubic diamond orthorhombic supercell
# with 8 atoms

basisstring = ""
origin = np.zeros(3)
bond = np.ones(3)*0.25
b = origin
basisstring += "basis %g %g %g " % (b[0],b[1],b[2])
b = bond
basisstring += "basis %g %g %g " % (b[0],b[1],b[2])

for i in range(3):
    b = 2*bond
    b[i] = 0
    basisstring += "basis %g %g %g " % (b[0],b[1],b[2])
    b += bond
    basisstring += "basis %g %g %g " % (b[0],b[1],b[2])

hmat = np.eye(3)

varlist = {
    "logsuffix":"ortho",
    "a":alat,
    "a1x":hmat[0,0],
    "a2x":hmat[0,1],
    "a2y":hmat[1,1],
    "a3x":hmat[0,2],
    "a3y":hmat[1,2],
    "a3z":hmat[2,2],
    "l":alat,
    "basis":basisstring,
    "nlat":3,
}

cmdargs = gen_varargs(varlist)
cij_ortho = calculate_cij(cmdargs)

# define the cubic diamond triclinic primitive cell
# with 2 atoms

basisstring = ""
origin = np.zeros(3)
bond = np.ones(3)*0.25
b = origin
basisstring += "basis %g %g %g " % (b[0],b[1],b[2])
b = bond
basisstring += "basis %g %g %g " % (b[0],b[1],b[2])

hmat1 = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]]).T/np.sqrt(2)

# rotate primitive cell to LAMMPS orientation
# (upper triangular)

qmat, rmat = np.linalg.qr(hmat1)
ss = np.diagflat(np.sign(np.diag(rmat)))
rot = ss @ qmat.T
hmat2 = ss @ rmat

varlist = {
    "logsuffix":"tri",
    "a":alat,
    "a1x":hmat2[0,0],
    "a2x":hmat2[0,1],
    "a2y":hmat2[1,1],
    "a3x":hmat2[0,2],
    "a3y":hmat2[1,2],
    "a3z":hmat2[2,2],
    "l":alat/2**0.5,
    "basis":basisstring,
    "nlat":5,
}

cmdargs = gen_varargs(varlist)
cij_tri = calculate_cij(cmdargs)

if me == 0:
    print("\nPython output:")
    print("C_ortho = \n",cij_ortho)
    print()
    print("C_tri = \n",cij_tri)
    print()

    cij_tri_rot = rotate_cij(cij_tri, rot.T)

    print("C_tri(rotated back) = \n",cij_tri_rot)
    print()

    print("C_ortho-C_tri = \n", cij_ortho-cij_tri_rot)
    print()
