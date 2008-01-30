# Stillinger-Weber parameters for CdTe
# Parameters are from:
#
#	Wang, Stroud, and Markworth, PhysRevB, 40,5, 1989
#
# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3, 
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q, tol

# Note that in LAMMPS, two-body parameters for element i and j are
# specified by the values on the i-j-j line. Two-body
# values on the i-j-i or i-i-j lines are ignored by LAMMPS, and so 
# are set to zero here.

Cd Cd Cd 	1.03 2.51  1.80  25.0  1.20  -0.333333333333
         	5.1726 0.8807  4.0  0.0 0.0
Cd Cd Te 	0.0  0.0  0.0  25.0  1.20  -0.333333333333
         	0.0 0.0 0.0 0.0 0.0
Cd Te Cd 	0.0  0.0  0.0  25.0  1.20  -0.333333333333
         	0.0 0.0 0.0 0.0  0.0
Cd Te Te 	1.03 2.51  1.80  25.0  1.20  -0.333333333333
         	7.0496 0.6022  4.0  0.0 0.0 
Te Cd Cd 	1.03 2.51  1.80  25.0  1.20  -0.333333333333
         	7.0496 0.6022  4.0  0.0 0.0
Te Cd Te 	0.0  0.0  0.0  25.0  1.20  -0.333333333333
         	0.0 0.0 0.0 0.0 0.0
Te Te Cd 	0.0  0.0  0.0  25.0  1.20  -0.333333333333
         	0.0 0.0 0.0 0.0 0.0
Te Te Te 	1.03 2.51  1.80  25.0  1.20  -0.333333333333
         	8.1415 0.6671  4.0  0.0 0.0 
