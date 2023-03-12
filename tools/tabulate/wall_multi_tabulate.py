#!/usr/bin/env python

from tabulate import WallTabulate
import os, sys
################################################################################
import math

k = 100.0
rzero = 4.0
epsilon = 0.02
sigma = 2.0
depth = 20.0
width = 2.0
r0 = 1.2

def harmonic_force(r):
    dr = r - rzero
    f = -2.0 * k * dr
    return f

def harmonic_energy(r):
    dr = r - rzero
    f = k * dr * dr
    return f

def lj126_force(r):
    f = -4.0*epsilon*(-12.0*math.pow(sigma/r,12.0)/r + 6.0*math.pow(sigma/r,6.0)/r)
    return f

def lj126_energy(r):
    f = 4.0*epsilon*(math.pow(sigma/r,12.0) - math.pow(sigma/r,6.0))
    return f

def morse_energy(r):
    ralpha = math.exp(-width*(r-r0))
    f = depth * (-1.0 + (1.0 - ralpha) * (1.0 - ralpha))
    return f

def morse_force(r):
    ralpha = math.exp(-width*(r-r0))
    f = -2.0 * depth * width * (1.0 -ralpha) * ralpha
    return f


################################################################################

if __name__ == "__main__":
    fname = 'walltab.dat'
    if os.path.exists(fname):
        os.remove(fname)
    sys.argv.append('--filename')
    sys.argv.append(fname)
    sys.argv.append('--num-points')
    sys.argv.append('400')
    sys.argv.append('--inner')
    sys.argv.append('0.01')
    sys.argv.append('--outer')
    sys.argv.append('4.0')
    wtable = WallTabulate(harmonic_energy, harmonic_force, units='real')
    wtable.run('HARMONIC')
    wtable = WallTabulate(lj126_energy, lj126_force, units='real')
    wtable.run('LJ126')
    wtable = WallTabulate(morse_energy, morse_force, units='real')
    wtable.run('MORSE')
