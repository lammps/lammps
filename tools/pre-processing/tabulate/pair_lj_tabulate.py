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
epsilon = 1.0
sigma = 1.0

def lj_energy(r):
    f = 4.0*epsilon*(math.pow(sigma/r,12.0) - math.pow(sigma/r,6.0))
    return f

def lj_force(r):
    epsilon = 1.0
    sigma = 1.0
    f = -4.0*epsilon*(-12.0*math.pow(sigma/r,12.0)/r + 6.0*math.pow(sigma/r,6.0)/r)
    return f
################################################################################

if __name__ == "__main__":
    ptable = PairTabulate(lj_energy, lj_force)
    ptable.run('LJ_11')
    epsilon = 1.0
    sigma = 1.5
    ptable.run('LJ_12')
