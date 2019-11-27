This package contains a LAMMPS implementation of coarse-grained
models of DNA, which can be used to model sequence-specific
DNA strands.

Please cite [1] and the relevant oxDNA articles in any publication 
that uses this package.

See the doc pages and [2,3,4] for the individual bond and pair styles. 
The packages contains also a new Langevin-type rigid-body integrator,
which has also its own doc page and is explained in [5].

[1] O. Henrich, Y. A. Gutierrez-Fosado, T. Curk, T. E. Ouldridge,
"Coarse-grained simulation of DNA using LAMMPS", 
Eur. Phys. J. E 41, 57 (2018).

[2] T. Ouldridge, A. Louis, J. Doye, "Structural, mechanical, 
and thermodynamic properties of a coarse-grained DNA model",
J. Chem. Phys. 134, 085101 (2011).

[3] T.E. Ouldridge, Coarse-grained modelling of DNA and DNA 
self-assembly, DPhil. University of Oxford (2011).

[4] B.E. Snodin, F. Randisi, M. Mosayebi, et al., Introducing
Improved Structural Properties and Salt Dependence into a Coarse-Grained
Model of DNA, J. Chem. Phys. 142, 234901 (2015).

[5] R. Davidchack, T. Ouldridge, M. Tretyakov, "New Langevin and 
gradient thermostats for rigid body dynamics", J. Chem. Phys. 142, 
144114 (2015).

Example input and data files can be found in
/examples/USER/cgdna/examples/oxDNA/ and /oxDNA2/. Python setup 
tools which create single straight or helical DNA strands as
well as DNA duplexes or arrays of duplexes can be found in
/examples/USER/cgdna/util/. A technical report with more information
on the models, the structure of the input and data file, the setup tool
and the performance of the LAMMPS-implementation of oxDNA can be found
in /doc/src/PDF/USER-CGDNA.pdf.

IMPORTANT NOTE: This package can only be used if LAMMPS is compiled
with the MOLECULE and ASPHERE packages.  These should be included in
the LAMMPS build by typing "make yes-asphere yes-molecule" prior to
the usual compilation (see the "Including/excluding packages" section
of the LAMMPS manual).

The creator of this package is:

Dr Oliver Henrich
University of Strathclyde, Glasgow, UK
oliver d o t henrich a t strath d o t ac d o t uk


--------------------------------------------------------------------------

** Bond styles provided by this package:

bond_oxdna_fene.cpp:  backbone connectivity, a modified FENE potential

bond_oxdna2_fene.cpp: corresponding bond style in oxDNA2 (see [3])

** Pair styles provided by this package:

pair_oxdna_excv.cpp:  excluded volume interaction between the nucleotides

pair_oxdna_stk.cpp:  stacking interaction between consecutive nucleotides
                     on the same strand

pair_oxdna_hbond.cpp:  hydrogen-bonding interaction between complementary
                       nucleotides on different strands, e.g. A-T and C-G

pair_oxdna_xstk.cpp:  cross-stacking interaction between nucleotides

pair_oxdna_coaxstk.cpp:  coaxial stacking interaction between nucleotides


pair_oxdna2_excv.cpp, pair_oxdna2_coaxstk.cpp:
                     corresponding pair styles in oxDNA2 (see [3])

pair_oxdna2_dh.cpp:  Debye-Hueckel electrostatic interaction between backbone
                     sites

** Fixes provided by this package:

fix_nve_dotc_langevin.cpp:  fix for Langevin-type rigid body integrator "C"
                            in above Ref. [3] 

fix_nve_dot.cpp:  NVE-type rigid body integrator without noise
