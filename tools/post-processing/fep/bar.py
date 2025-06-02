#!/usr/bin/env python
# bar.py - Bennet's acceptance ratio method for free energy calculation

import sys
import math

if len(sys.argv) < 4:
    print("Bennet acceptance ratio method")
    print("usage: bar.py temperature datafile01 datafile10 [delf_lo delf_hi]")
    print("  datafile01 contains (U_1 - U_0)_0 in 2nd column")
    print("  datafile10 contains (U_0 - U_1)_1 in 2nd column")
    print("  (first column is index, time step, etc. and is ignored)")
    print("  delf_lo and delf_hi are optional guesses bracketing the solution")
    sys.exit()

if len(sys.argv) == 6:
    delf_lo = float(sys.argv[4])
    delf_hi = float(sys.argv[5])
else:
    delf_lo = -50.0
    delf_hi =  50.0
    
def fermi(x):
    if x > 100:
        return 0.0
    else:
        return 1.0 / (1.0 + math.exp(x))

def avefermi(eng, delf):
    ave = 0.0
    n = 0
    for du in eng:
        ave += fermi((du + delf) / rt)
        n += 1
    return ave / n

def bareq(delf):
    ave0 = avefermi(eng01, -delf)
    ave1 = avefermi(eng10,  delf)
    return ave1 - ave0

def bisect(func, xlo, xhi, xtol = 1.0e-4, maxit = 20):
    if xlo > xhi:
        aux = xhi
        xhi = xlo
        xlo = aux
    if func(xlo) * func(xhi) > 0.0:
        print("error: root not bracketed by interval")
        sys.exit(2)
    for i in range(maxit):
        sys.stdout.write('.')
        sys.stdout.flush()
        xmid = (xlo + xhi) / 2.0
        if func(xlo) * func(xmid) < 0.0:
            xhi = xmid
        else:
            xlo = xmid
        if xhi - xlo < xtol:
            return xmid
    return xmid

print("Bennet acceptance ratio method")
print(sys.argv[1], " K")
rt = 0.008314 / 4.184 * float(sys.argv[1])

eng01 = []                                # read datafiles
with open(sys.argv[2], 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        eng01.append(float(line.split()[1]))

eng10 = []
with open(sys.argv[3], 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        eng10.append(float(line.split()[1]))

sys.stdout.write("solving")
delf = bisect(bareq, delf_lo, delf_hi)
print(".")

ave0 = avefermi(eng01, -delf)
ave1 = avefermi(eng10,  delf)

print("<...>0 = ", ave0)
print("<...>1 = ", ave1)
print("deltaA = ", delf)
