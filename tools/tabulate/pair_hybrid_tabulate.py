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

from tabulate import PairTabulate

################################################################################
import math

# hybrid potential using Morse for repulsive and LJ for attractive.
# using tanh() as smooth switching function between the two.

def hybrid_energy(r):
    depth = 1.0
    width = 5.0
    epsilon = 1.0
    sigma = 1.0
    rzero = 1.2

    f1 = 4.0*epsilon*(math.pow(sigma/r,12.0) - math.pow(sigma/r,6.0))

    ralpha = math.exp(-width*(r-rzero))
    f2 = depth * (-1.0 + (1.0 - ralpha) * (1.0 - ralpha))

    switch = 0.5*math.tanh(10.0*(r - rzero)) + 0.5
    f = switch*f1 + (1.0 - switch)*f2
    return f

################################################################################

if __name__ == "__main__":
    ptable = PairTabulate(hybrid_energy, units='lj', comment='Morse repulsion + LJ attraction')
    ptable.run('MORSE_LJ')
