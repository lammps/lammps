#!/usr/bin/env python
# fdti.py - integrate compute fep results using the trapezoidal rule

import sys
import math

if len(sys.argv) < 3:
    print("Finite Difference Thermodynamic Integration (Mezei 1987)")
    print("Trapezoidal integration of compute_fep results at equally-spaced points")
    print("usage: fdti.py temperature hderiv < out.fep")
    sys.exit()

rt = 0.008314 / 4.184 * float(sys.argv[1]) # in kcal/mol
hderiv = float(sys.argv[2])

line = sys.stdin.readline()
while line.startswith("#"):
    line = sys.stdin.readline()

v = 1.0
tok = line.split()
if len(tok) == 4:
    v = float(tok[3])
lo = -rt * math.log(float(tok[2]) / v)
    
i = 1
sum = 0.0
for line in sys.stdin:
    tok = line.split()
    if len(tok) == 4:
        v = float(tok[3])
    hi = - rt * math.log(float(tok[2]) / v)
    sum += (hi + lo) / (2 * hderiv)
    lo = hi
    i += 1

print(sum/(i - 1))    # int_0^1: divide by i - 1 == multiply by delta
