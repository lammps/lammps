#!/usr/bin/env python

# Calculate a table of pairwise energies and forces between atoms in the
# protein and a chaperone provided in the supplemental materials section of:
# AI Jewett, A Baumketner and J-E Shea, PNAS, 101 (36), 13192-13197, (2004)
# This is stored in a tabulated force field with a singularity at a distance R.
#
# To calculate the table for interaction between
# ...the chaperone and a hydrophobic bead (2004 PNAS paper), use this table:
#   ./calc_chaperone_table.py 1.0 1.0 6.0 0.475 0.0 5.9 1181
# ...the chaperone and a hydrophilic bead (2004 PNAS paper), use this table:
#   ./calc_chaperone_table.py 1.0 1.0 6.0 0.0   0.0 5.9 1181
# ...the chaperone and a hydrophobic bead (2006 JMB paper), use this table:
#   ./calc_chaperone_table.py 1.0 1.0 3.0 0.60  3.1 8.0 981 True
# ...the chaperone and a hydrophilic bead (2006 JMB paper), use this table:
#   ./calc_chaperone_table.py 1.0 1.0 3.0 0.0   3.1 8.0 981 True

from math import *
import sys

def U(r, eps, sigma, R, h):
    #print('r='+str(r)+' eps='+str(eps)+' s='+str(sigma)+' R='+str(R)+' h='+str(h))
    # Formula is undefined at r=0, but you can take the limit:
    if r <= 0:
        return 4.0*pi*R*R*4.0*eps*(pow((sigma/R), 12.0)
                               - h*pow((sigma/R), 6.0))
    xp = sigma/(r+R)
    xm = sigma/(r-R)
    term10 = pow(xm, 10.0) - pow(xp, 10.0)
    term4  = pow(xm, 4.0)  - pow(xp, 4.0)
    return 4.0*pi*eps*(R/r) * (0.2*term10 - 0.5*h*term4)

def F(r, eps, sigma, R, h):
    # Formula is undefined at r=0, but you can take the limit:
    if r <= 0:
        return 0.0
    product_term_a = U(r, eps, sigma, R, h) / r
    ixp = (r+R)/sigma
    ixm = (r-R)/sigma
    dix_dr = 1.0/sigma
    term10 = (10.0/sigma)*(pow(ixm, -11.0) - pow(ixp, -11.0))
    term4  =  (4.0/sigma)*(pow(ixm, -5.0)  - pow(ixp, -5.0))
    product_term_b = 4.0*eps*pi*(R/r) * (0.2*term10 - 0.5*h*term4)
    return product_term_a + product_term_b


class InputError(Exception):
    """ A generic exception object containing a string for error reporting.

    """
    def __init__(self, err_msg):
        self.err_msg = err_msg
    def __str__(self):
        return self.err_msg
    def __repr__(self):
        return str(self)

if len(sys.argv) < 8:
    sys.stderr.write("Error: expected 7 arguments:\n"
                     "\n"
                     "Usage: "+sys.argv[0]+" epsilon sigma R h rmin rmax N\n\n")
    sys.exit(-1)

epsilon = float(sys.argv[1])
sigma   = float(sys.argv[2])
R       = float(sys.argv[3])
h       = float(sys.argv[4])
rmin    = float(sys.argv[5])
rmax    = float(sys.argv[6])
N       = int(sys.argv[7])

subtract_Urcut = False
if len(sys.argv) == 9:
    subtract_Urcut = True
rcut    = rmax

for i in range(0,N):
    r = rmin + i*(rmax-rmin)/(N-1)
    U_r = U(r, epsilon, sigma, R, h)
    F_r = F(r, epsilon, sigma, R, h)
    if subtract_Urcut:
        U_r -= U(rcut, epsilon, sigma, R, h)
        if (r >= rcut) or (i==N-1):
            U_r = 0.0
            F_r = 0.0
    print(str(i+1)+' '+str(r)+' '+str(U_r)+' '+str(F_r))

