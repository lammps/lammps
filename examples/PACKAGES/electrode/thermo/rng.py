#!/usr/bin/env python3

import numpy as np

rng = np.random.default_rng(7)
vals = rng.standard_normal(10 ** 6 + 1)
np.savetxt("random.inp", vals)
