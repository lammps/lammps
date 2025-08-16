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

from tabulate import WallTabulate

################################################################################
import math

k = 100.0
rzero = 4.0

def harmonic_force(r):
    dr = r - rzero
    f = -2.0 * k * dr
    return f

def harmonic_energy(r):
    dr = r - rzero
    f = k * dr * dr
    return f

################################################################################

if __name__ == "__main__":
    wtable = WallTabulate(harmonic_energy, harmonic_force, units='real')
    wtable.run('HARMONIC')
