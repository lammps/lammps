# DATE: 2008-01-29 CONTRIBUTOR: Unknown CITATION: Bere and Serra, Phil Mag 86, 2159 (2006)
# GaN potential: A. Bere, and A. Serra, Phil. Mag. 86, 2159(2006)
# note that the parameters for this literature potential are pairwise
# so that there are some flexibility in the way the 
# parameters can be entered. As one way, we assume that 
# lambda_ijk is equal to lambda_ik and eps_ijk is 
# equal to sqrt(lambda_ij*eps_ij*lambda_ik*eps_ik)/lambda_ik, 
# and all other parameters in the ijk line are for ik.

# These entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless
#
#         eps      sigma a lambda gamma  cos(theta)     A    B    p   q  tol
#
Ga Ga Ga 1.2000000 2.100 1.6 32.5 1.2 -0.333333333333 7.917 0.72 4.0 0.0 0.0
N  N  N  1.2000000 1.300 1.8 32.5 1.2 -0.333333333333 7.917 0.72 4.0 0.0 0.0
Ga Ga N  1.6136914 0.0   0.0 32.5 0.0 -0.333333333333 0.0   0.0  0.0 0.0 0.0
Ga N  N  2.1700000 1.695 1.8 32.5 1.2 -0.333333333333 7.917 0.72 4.0 0.0 0.0
N  Ga Ga 2.1700000 1.695 1.8 32.5 1.2 -0.333333333333 7.917 0.72 4.0 0.0 0.0
N  Ga N  1.6136914 0.0   0.0 32.5 0.0 -0.333333333333 0.0   0.0  0.0 0.0 0.0
N  N  Ga 1.6136914 0.0   0.0 32.5 0.0 -0.333333333333 0.0   0.0  0.0 0.0 0.0
Ga N  Ga 1.6136914 0.0   0.0 32.5 0.0 -0.333333333333 0.0   0.0  0.0 0.0 0.0
