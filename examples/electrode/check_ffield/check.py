#!/bin/python3

import numpy as np

ref = np.loadtxt("ref.csv", skiprows=1)
ffield = np.loadtxt("ffield.csv", skiprows=1)
ffield_flip = np.loadtxt("ffield_flip.csv", skiprows=1)

for matrix in [ref, ffield, ffield_flip]:
    assert matrix.shape == (6, 3), matrix.shape

# ref is close to ffield
diff = np.sum(np.abs(ref - ffield))
assert diff < 3e-4, diff

# ffield and ffield_flip are identical
diff = np.sum(np.abs(ffield_flip - ffield))
assert diff < 1e-9, diff
