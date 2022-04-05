#!/bin/python3

import numpy as np

onestep = np.loadtxt("onestep.csv", skiprows=1)
twostep = np.loadtxt("twostep.csv", skiprows=1)

for matrix in [onestep, twostep]:
    assert matrix.shape == (288, 288), matrix.shape

diff = abs(np.sum(onestep - twostep))
assert diff < 1e-11, diff
