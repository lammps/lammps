# Pair nm/split:

Template: Variables defined below. 
pair_style	nm/split r_c
pair_coeff	* * E_0 r_0 n m

Example: Core-softening of LJ potential:
pair_style      nm/split 1.12246
pair_coeff      * * 1.0 1.1246 12 6
pair_coeff      * * 1.0 1.1246 11 6
pair_coeff      * * 1.0 1.1246 10 6
				.
				.
				.
pair_coeff      * * 1.0 1.1246  2 6

Corresponding Files: pair_nm_split.cpp  pair_nm_split.h

Description:
Style nm-split computes monomer-monomer interactions based on the generalized Lennard-Jones (LJ) potential by (Clarke). Here n is the parameter that controls the strength of the repulsive component of the interaction, m controls the strength of the attractive component. Note, nm/split is a modification to the already existing pair style nm/cut.
Unlike nm/cut, style nm/split employs the standard LJ (n=12, m=6) potential above the 2^1/6σ (i.e. when forces are attractive). Applications and references to style nm-split can be found in (Hoy) where nm-split was used in conjunction with the fix bond/swap command to equilibrate long-chain Kremer-Grest bead-spring polymer melts. 

For all of the nm pair styles, the following coefficients must be defined for each pair of atoms types via the pair_coeff command as in the examples above, or in the data file or restart files read by the read_data or read_restart commands.
-Depth of energy minimum  E_0
-Position of energy minimum and r_0 ( r_0 = 2^1/6 σ is the standard LJ value)
-Cutoff distance r_c
-n (repulsive)
-m (attractive)


Mixing, shift, table, tail correction, restart, rRESPA info
These pair styles do not support mixing. Thus, coefficients for all ij pairs must be specified explicitly.
All of the nm pair styles supports the pair_modify shift option for the energy of the pair interaction.
All of the nm pair styles support the pair_modify tail option for adding a long-range tail correction to the energy and pressure for the N-M portion of the pair interaction.
All of the nm pair styles write their information to binary restart files, so pair_style and pair_coeff commands do not need to be specified in an input script that reads a restart file.
All of the nm pair styles can only be used via the pair keyword of the run_style respa command. They do not support the inner, middle, outer keywords.

References:
(Clarke) Clarke and Smith, J Chem Phys, 84, 2290 (1986).
(Hoy) Hoy and Kröger, Phys. Rev. Lett. 124, 147801 (2020). 

