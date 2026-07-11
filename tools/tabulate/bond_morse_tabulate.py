#!/usr/bin/env python
# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   LAMMPS Development team: developers@lammps.org
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

from tabulate import BondTabulate

################################################################################
import math

def morse_energy(r):
    depth = 1.0
    width = 2.0
    rzero = 1.2
    ralpha = math.exp(-width*(r-rzero))
    f = depth * (-1.0 + (1.0 - ralpha) * (1.0 - ralpha))
    return f

def morse_force(r):
    depth = 1.0
    width = 2.0
    rzero = 1.2
    ralpha = math.exp(-width*(r-rzero))
    f = -2.0 * depth * width * (1.0 -ralpha) * ralpha
    return f

################################################################################

if __name__ == "__main__":
    btable = BondTabulate(morse_energy, morse_force, units='lj')
    btable.run('MORSE')
