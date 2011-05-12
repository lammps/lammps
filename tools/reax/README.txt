=== ReaxFF tools ===
===============================

The programs in this folder can be used to analyze the
output of simulations using the ReaxFF potentials;

mol_fra.c: reads the output of fix reax/bonds
	   and identifies fragments
   Compile it using a C compiler
   To test, run it with Cutoff.dic and bonds.reax
   Contact: Aidan Thompson

bondConnectCheck.f90: reads the output of fix reax/bonds.
   Does not do fragment analysis.
   Compile it using FORTRAN compiler
   To test, run it with bonds.reax
   Contact: Paul Liangliang Huang <lhuang4@ncsu.edu>

reaxc_bond.pl: reads the bonding information in the
                .trj file produced by pair_style reax/c and
                outputs molecule counts for each frame.  

