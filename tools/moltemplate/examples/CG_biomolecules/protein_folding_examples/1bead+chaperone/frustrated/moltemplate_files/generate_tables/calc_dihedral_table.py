#!/usr/bin/env python

# Calculate a table of dihedral angle interactions used in the alpha-helix
# and beta-sheet regions of the frustrated protein model described in
#  provided in figure 8 of the supplemental materials section of:
# AI Jewett, A Baumketner and J-E Shea, PNAS, 101 (36), 13192-13197, (2004)
# Note that the "A" and "B" parameters were incorrectly reported to be 
# 5.4*epsilon and 6.0*epsilon.  The values used were 5.6 and 6.0 epsilon.
# The phiA and phiB values were 57.29577951308232 degrees (1 rad) 
# and 180 degrees, respectively.  Both expA and expB were 6.0.
#
# To generate the table used for the alpha-helix (1 degree resolution) use this:
#  ./calc_dihedral_table.py 6.0 57.29577951308232 6  5.6 180 6  0.0 359  360
# To generate the table used for the beta-sheets (1 degree resolution) use this:
#  ./calc_dihedral_table.py 5.6 57.29577951308232 6  6.0 180 6  0.0 359  360
#
# (If you're curious as to why I set the location of the minima at phi_alpha
#  to 1.0 radians (57.2957795 degrees), there was no particularly good reason.
#  I think the correct value turns out to be something closer to 50 degrees.)


from math import *
import sys


# The previous version included the repulsive core term
def U(phi, A, phiA, expA, B, phiB, expB, use_radians=False):
    conv_units = pi/180.0
    if use_radians:
        conv_units = 1.0
    termA = pow(cos(0.5*(phi-phiA)*conv_units), expA)
    termB = pow(cos(0.5*(phi-phiB)*conv_units), expB)
    return -A*termA - B*termB

# The previous version included the repulsive core term
def F(phi, A, phiA, expA, B, phiB, expB, use_radians=False):
    conv_units = pi/180.0
    if use_radians:
        conv_units = 1.0
    termA = (0.5*sin(0.5*(phi-phiA)*conv_units) * 
             expA * pow(cos(0.5*(phi-phiA)*conv_units), expA-1.0))
    termB = (0.5*sin(0.5*(phi-phiB)*conv_units) * 
             expB * pow(cos(0.5*(phi-phiB)*conv_units), expB-1.0))
    return -conv_units*(A*termA + B*termB)

if len(sys.argv) != 10:
    sys.stderr.write("Error: expected 9 arguments:\n"
                     "\n"
                     "Usage: "+sys.argv[0]+" A phiA expA B phiB expB phiMin phiMax N\n\n")
    sys.exit(-1)

A       = float(sys.argv[1])
phiA    = float(sys.argv[2])
expA    = float(sys.argv[3])
B       = float(sys.argv[4])
phiB    = float(sys.argv[5])
expB    = float(sys.argv[6])
phi_min = float(sys.argv[7])
phi_max = float(sys.argv[8])
N       =   int(sys.argv[9])

for i in range(0,N):
    phi = phi_min + i*(phi_max - phi_min)/(N-1)
    U_phi = U(phi, A, phiA, expA, B, phiB, expB, use_radians=False)
    F_phi = F(phi, A, phiA, expA, B, phiB, expB, use_radians=False)
    print(str(i+1)+' '+str(phi)+' '+str(U_phi)+' '+str(F_phi))

