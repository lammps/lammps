#!/usr/bin/env python

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
