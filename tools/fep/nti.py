#!/usr/bin/env python
# nti.py - integrate compute fep results using the trapezoidal rule

import sys
import math

if len(sys.argv) < 3:
    print("Thermodynamic Integration with Numerical Derivative")
    print("Trapezoidal integration of compute_fep results at equally-spaced points")
    print("usage: nti.py temperature hderiv < out.fep")
    sys.exit()

hderiv = float(sys.argv[2])

line = sys.stdin.readline()
while line.startswith("#"):
    line = sys.stdin.readline()

tok = line.split()
lo = float(tok[1])

i = 0
sum = 0.0
for line in sys.stdin:
    tok = line.split()
    hi = float(tok[1])
    sum += (hi + lo) / (2 * hderiv)
    lo = hi
    i += 1

print(sum/(i - 1))    # int_0^1: divide by i - 1 == multiply by delta
