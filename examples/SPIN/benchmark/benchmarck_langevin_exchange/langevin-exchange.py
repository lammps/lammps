#!/usr/bin/env python3

import numpy as np, pylab, tkinter
import matplotlib.pyplot as plt
import mpmath as mp

mub=5.78901e-5          # Bohr magneton (eV/T)
kb=8.617333262145e-5    # Boltzman constant (eV/K)
J0=0.05                 # per-neighbor exchange interaction (eV)
z=6                     # number of NN (bcc)
g=2.0                   # Lande factor (adim)
Hz=10.0                # mag. field (T)

#Definition of the Langevin function
def func(sm,t):
    return mp.coth(z*J0*sm/(kb*t))-1.0/(z*J0*sm/(kb*t))

npoints=200
tolerance=1e-5
ti=0.01
tf=2000.0
sz=1.0
szg=0.5
for i in range (0,npoints):
    temp=ti+i*(tf-ti)/npoints
    count=0
    sz=1.0
    szg=0.5
    while (abs(sz-szg)) >= tolerance:
        sz=szg
        szg=func(sz,temp)
        count+=1
    emag=-z*J0*sz*sz
    print('%d %lf %lf %lf %lf' % (temp,szg,sz,emag,count))
