#!/usr/bin/env python
# fep.py - calculate free energy from compute fep results

import sys, math

if len(sys.argv) < 2:
    print "Free Energy Perturbation"
    print "usage: fep.py temperature < fep.lmp"
    sys.exit()

rt = 0.008314 / 4.184 * float(sys.argv[1])

v = 1.0
sum = 0.0
for line in sys.stdin:
    if line.startswith("#"):
        continue
    tok = line.split()
    if len(tok) == 4:
        v = float(tok[3])
    sum += math.log(float(tok[2]) / v)

print -rt * sum
