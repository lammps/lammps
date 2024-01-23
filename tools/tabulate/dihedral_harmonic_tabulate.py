#!/usr/bin/env python

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
