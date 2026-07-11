# DATE: 2007-06-11 UNITS: metal CONTRIBUTOR: Unknown CITATION: Wang, Stroud and Markworth, Phys Rev B, 40, 3129 (1989).

# CdTe Stillinger-Weber potential: Z. Q. Wang, D. Stroud,
# and A. J. Markworth, Phys. Rev. B, 40, 3129(1989).

# The Stillinger-Weber parameters given in the literature are pair
# specific. While most of the parameters are indeed pairwise parameters
# according to their definition, the parameters epsilon and lambda
# should be viewed as three-body dependent. Here we assume that the
# the three-body epsilon and lambda is a geometric mean of the pairwise
# epsilon and lambda.

# In lammps, the parameters for the ij pair are entered in
# the ijj three-body line. There is no unique way to convert pair
# parameters to three body parameters so the example here represents
# only one way. The three-body parameters epsilon_ijk can be calculated
# from the literature pair parameters using epsilon_ijk =
# sqrt(lambda_ij*epsilon_ij*lambda_ik*epsilon_ik)/lambda_ik, and the
# results are directly entered in this table. Obviously, this
# conversion does not change the two-body parameters epsilon_ijj. 

# The twobody ik pair parameters are entered on the i*k lines, where *
# can be any species. This is consistent with the LAMMPS requirement
# that twobody ik parameters be defined on the ikk line. Entries on all
# the other i*k lines are ignored by LAMMPS

# These entries are in LAMMPS "metal" units: epsilon = eV;
# sigma = Angstroms; other quantities are unitless

#     epsilon sigma a lambda gamma   cos(theta)     A      B     p   q  tol
Cd Cd Cd 1.03 2.51 1.80 25.0 1.20 -0.333333333333 5.1726 0.8807 4.0 0.0 0.0
Te Te Te 1.03 2.51 1.80 25.0 1.20 -0.333333333333 8.1415 0.6671 4.0 0.0 0.0
Cd Cd Te 1.03 0.0  0.0  25.0 0.0  -0.333333333333 0.0    0.0    0.0 0.0 0.0
Cd Te Te 1.03 2.51 1.80 25.0 1.20 -0.333333333333 7.0496 0.6022 4.0 0.0 0.0
Te Cd Cd 1.03 2.51 1.80 25.0 1.20 -0.333333333333 7.0496 0.6022 4.0 0.0 0.0
Te Cd Te 1.03 0.0  0.0  25.0 0.0  -0.333333333333 0.0    0.0    0.0 0.0 0.0
Te Te Cd 1.03 0.0  0.0  25.0 0.0  -0.333333333333 0.0    0.0    0.0 0.0 0.0
Cd Te Cd 1.03 0.0  0.0  25.0 0.0  -0.333333333333 0.0    0.0    0.0 0.0 0.0
