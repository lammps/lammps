# DATE: 2026-01-27 CONTRIBUTOR: Roman Groger
# CITATION: R. Groger, J. Fikar, Acta Mater. XX, XXXX (2026)
#
# Modified Tersoff potential parameters for GaN describing short-range interactions in ionic-covalent
# model of GaN developed within the framework of the Streitz-Mintmire formalism. This file must be
# combined with the electrostatic energy to form a hybrid potential. The electrostatic parameters
# of Ga and N are stored in the file GaN.streitz.
#
# The parameters of homobonds are taken from Nord et al., J. Phys. Cond. Mat 15, 5649 (2003).
# Fitting of Ga-N parameters and the development of the ionic-covalent model for GaN is described in
#   Groger & Fikar, Acta Mater. XX, XXXX (2026)
#
# Multiple entries can be added to this file, LAMMPS reads the ones it needs.
# These entries are in LAMMPS "metal" units:
#   A,B = eV; lambda1,lambda2 = 1/Angstroms; R,D = Angstroms
#   other quantities are unitless
#
# Format of a single entry (one or more lines):
#   element 1, element 2, element 3, beta, alpha, h, eta, beta_ters, lambda2, B, R, D, lambda1, A, n, c1,  c2,  c3,  c4,  c5

Ga Ga Ga 1.0 1.846 -0.3013 1.0 1.0 1.4497 410.132 2.87 0.15 1.60916 535.199 1.0 0.007874 0.05149559604622222 0.5625 0.0 0.0
N N N 1.0 0.0 -0.045238 1.0 1.0 2.38426 423.769 2.2 0.2 3.55779 1044.77 1.0 0.76612 0.5998480604393894 0.040690958400000005 0.0 0.0
Ga Ga N 1.0 0.0 -0.42139396832132586 1.0 1.0 0.0 0.0 2.9 0.2 0.0 0.0 1.0 0.0009358372029931801 0.07556968022478684 0.42407692241770484 0.0 0.0
Ga N N 1.0 0.0 -0.42139396832132586 1.0 1.0 2.7220181847973004 4626.922179188084 2.9 0.2 3.0086205179630774 7426.478709035391 1.0 0.0009358372029931801 0.07556968022478684 0.42407692241770484 0.0 0.0
N Ga Ga 1.0 0.0 -0.42139396832132586 1.0 1.0 2.7220181847973004 4626.922179188084 2.9 0.2 3.0086205179630774 7426.478709035391 1.0 0.0009358372029931801 0.07556968022478684 0.42407692241770484 0.0 0.0
N Ga N 1.0 0.0 -0.045238 1.0 1.0 0.0 0.0 2.2 0.2 0.0 0.0 1.0 0.76612 0.5998480604393894 0.040690958400000005 0.0 0.0
N N Ga 1.0 0.0 -0.42139396832132586 1.0 1.0 0.0 0.0 2.9 0.2 0.0 0.0 1.0 0.0009358372029931801 0.07556968022478684 0.42407692241770484 0.0 0.0
Ga N Ga 1.0 1.846 -0.3013 1.0 1.0 0.0 0.0 2.87 0.15 0.0 0.0 1.0 0.007874 0.05149559604622222 0.5625 0.0 0.0
