#!/usr/bin/env python

from tabulate import AngleTabulate

################################################################################
import math

def harmonic_energy(theta):
    k = 100.0
    thetazero = 120.0
    # the force constant in LAMMPS is in energy per radians^2 so convert from degrees to radians
    deg2rad = math.pi / 180.0
    t = (theta - thetazero) * deg2rad
    f = k * t * t
    return f

################################################################################

if __name__ == "__main__":
    atable = AngleTabulate(harmonic_energy, units='real')
    atable.run('HARM')
