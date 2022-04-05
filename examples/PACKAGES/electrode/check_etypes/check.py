#!/bin/python3

import numpy as np

ref = np.loadtxt("ref.csv", skiprows=1)
etypes = np.loadtxt("etypes.csv", skiprows=1)

for matrix in [ref, etypes]:
    assert matrix.shape == (6, 3), matrix.shape

diff0 = abs(np.sum(ref - etypes))
assert diff0 < 1e-8, diff0
