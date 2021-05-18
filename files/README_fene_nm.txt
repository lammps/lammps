Bond fene/nm

Example

bond_style      fene/nm
bond_coeff      1 2.25344 1.5 1.0 1.12246 2 6

Corresponding Files: bond_fene_nm.cpp bond_fene_nm.h

Description:
Style fene/nm is based on the standard fene which consists of two components: 1. an attractive logarithmic term and 2. a repulsive Lennard-Joines (12-6) term. 

Unlike the standard fene style, fene/nm substitutes the standard (12-6) LJ term with a modified generalized (n-m) LJ potential (style nm/split).

The following coefficients must be defined for each bond type via the bond_coeff command as in the example above, or in the data file or restart files read by the read_data or read_restart commands:
-Spring Constant K
-Maximum bond extension R_0
-Energy
-Monomer Diameter σ 

 
