#!/usr/bin/env python

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
