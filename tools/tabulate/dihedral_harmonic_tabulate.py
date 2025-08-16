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

from tabulate import DihedralTabulate

################################################################################
import math

def harmonic_energy(phi):
    k = 100.0
    # the force constant in LAMMPS is in energy per radians^2 so convert from degrees to radians
    deg2rad = math.pi / 180.0
    f = k * (1 - math.cos(2.0 * deg2rad * phi))
    return f

################################################################################

if __name__ == "__main__":
    dtable = DihedralTabulate(harmonic_energy, units='metal')
    dtable.run('HARM')
