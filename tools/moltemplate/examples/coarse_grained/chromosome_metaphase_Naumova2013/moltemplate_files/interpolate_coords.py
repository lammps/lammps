#!/usr/bin/env python


err_msg = """
Usage:

   interpolate_coords.py Ndesired [scale] < coords_orig.raw > coords.raw

Example:

   interpolate_coords.py 30118 3.0 < coords_orig.raw > coords.raw

   # (Note: 30117 ~= 128000/4.25, but using 30118 makes interpolation cleaner.
   #        See the supplemental section of Naumova et al Science 2013, p 18.)

"""


import sys
from math import floor

# Parse the argument list:
if len(sys.argv) <= 1:
    sys.stderr.write("Error:\n\nTypical Usage:\n\n"+err_msg+"\n")
    exit(1)

n_new = int(sys.argv[1])

if len(sys.argv) > 2:
    scale = float(sys.argv[2])
else:
    scale = 1.0

coords_orig = []

lines = sys.stdin.readlines()

for line in lines:
    tokens = line.split()
    if (len(tokens) > 0):
        coords_orig.append(list(map(float, tokens)))
        g_dim = len(tokens)

n_orig = len(coords_orig)

if n_orig < 2:
    sys.stderr.write("Error:\n\nInput file contains less than two lines of coordinates\n")
    exit(1)

if n_new < 2:
    sys.stderr.write("Error:\n\nOutput file will contain less than two lines of coordinates\n")
    exit(1)

coords_new = [[0.0 for d in range(0, g_dim)] for i in range(0, n_new)]

for i_new in range(0, n_new):
    I_orig = (i_new) * (float(n_orig-1) / float(n_new-1))
    i_orig = int(floor(I_orig))
    i_remainder = I_orig - i_orig

    if (i_new < n_new-1):
        for d in range(0, g_dim):
            coords_new[i_new][d] = scale*(coords_orig[i_orig][d]
                                          +
                                          i_remainder*(coords_orig[i_orig+1][d]-
                                                       coords_orig[i_orig][d]))
    else:
        for d in range(0, g_dim):
            coords_new[i_new][d] = scale*coords_orig[n_orig-1][d]

    # print the coordates
    for d in range(0, g_dim-1):
        sys.stdout.write(str(coords_new[i_new][d]) + ' ')
    sys.stdout.write(str(coords_new[i_new][g_dim-1]) + "\n")
