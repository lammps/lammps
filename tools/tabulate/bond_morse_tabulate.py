#!/usr/bin/env python

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
