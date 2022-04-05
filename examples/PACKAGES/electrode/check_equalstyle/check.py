#!/bin/python3

import numpy as np

const = np.loadtxt("const.csv", skiprows=1)
equal = np.loadtxt("equal.csv", skiprows=1)
ramp = np.loadtxt("ramp.csv", skiprows=1)

for matrix in [const, equal, ramp]:
    assert matrix.shape == (6, 3), matrix.shape

# equal follows const
diff0 = abs(np.sum(const - equal))
assert diff0 < 1e-8, diff0

# equal starts same as ramp but ends differently
diff1 = abs(np.sum(equal[0] - ramp[0]))
assert diff1 < 1e-8, diff1
diff2 = np.sum(abs(equal[-1] - ramp[-1]))
assert diff2 > 1e-8, diff2
