#!/usr/bin/env python3

import numpy as np, pylab, tkinter
import matplotlib.pyplot as plt
import mpmath as mp

mub=5.78901e-5          # Bohr magneton (eV/T)
kb=8.617333262145e-5    # Boltzman constant (eV/K)
g=2.0                   # Lande factor (adim)

Hz=10.0                # mag. field (T)

#Definition of the Langevin function
def func(t):
    return mp.coth(g*mub*Hz/(kb*t))-1.0/(g*mub*Hz/(kb*t))

npoints=200
ti=0.01
tf=20.0
for i in range (0,npoints): 
    temp=ti+i*(tf-ti)/npoints
    print('%lf %lf %lf' % (temp,func(temp),-g*mub*Hz*func(temp)))
