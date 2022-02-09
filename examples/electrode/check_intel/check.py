#!/bin/python3

import numpy as np

nointel = np.loadtxt("ref_nointel.csv", skiprows=1)
intel = np.loadtxt("ref_intel.csv", skiprows=1)

for matrix in [nointel, intel]:
    assert matrix.shape == (6, 3), matrix.shape

diff0 = np.sum(np.abs(nointel - intel))
assert diff0 < 5e-4, diff0
