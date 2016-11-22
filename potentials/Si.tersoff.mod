# DATE: 2013-07-26 CONTRIBUTOR: Unknown CITATION: Kumagai, Izumi, Hara and Sakai, Comp Mat Sci, 39, 457 (2007)

# Modified Tersoff potential (named MOD potential) parameters for various elements and mixtures
# multiple entries can be added to this file, LAMMPS reads the ones it needs
# these entries are in LAMMPS "metal" units:
#   A,B = eV; lambda1,lambda2 = 1/Angstroms; R,D = Angstroms
#   other quantities are unitless

# MOD. This is the Si parameterization from the paper:
# [1] T. Kumagai, S. Izumi, S. Hara, S. Sakai, Comp. Mat. Sci., 39, 457 (2007) 

# format of a single entry (one or more lines):
#   element 1, element 2, element 3, 
#   beta, alpha, h, eta, 
#   beta_ters, lambda2, B, R, D, lambda1, A, n
#   c1,  c2,  c3,  c4,  c5  

# Notes
# 1) beta_ters must be equal to 1.0
# 2) R = (R1+R2)/2, where R1 and R2 determined at [1] (R1 = 2.7A, R2 = 3.3A)
# 3) D = (R2-R1)/2, where R1 and R2 determined at [1] (R1 = 2.7A, R2 = 3.3A)
# 4) n = 1.0/(2*delta), where delta determined at [1] (eta*delta = 0.53298909)

#mod
Si  Si  Si  1.0  2.3890327  -.365  1.0 
            1.0    1.345797  121.00047  3.0   0.3   3.2300135  3281.5905  .93810551 
            0.20173476  730418.72  1000000.0  1.0  26.0 
