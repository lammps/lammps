This example was intended to demonstrate the flexibility of LAMMPS and
moltemplate.

This is a relatively complex example containing two different types of
coarse-grained (united-atom) molecules.  This simulation uses the 3-body
(non-pairwise-additive) coarse-grained "mW" water model:
Molinero, V. and Moore, E.B., J. Phys. Chem. B 2009, 113, 4008-4016
Simulations using the "mW" water model can be several orders of magnitude
faster than simulations using simple all-atom models such as SPCE or TIP3P.

The united-atom TraPPE force field was used for the cyclododecane molecules.

Any force-field available in LAMMPS can be used with moltemplate. New force-fields are added by end users regularly. For a current list, see:
http://lammps.sandia.gov/doc/Section_commands.html#pair-style-potentials

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step 2)
README_run.sh


-------- REQUIREMENTS: ---------
  This example requires the "MANYBODY" package.
  If lammps complains of a missing pair style enter "make yes-MANYBODY"
  into the shell before compiling lammps.  For details see:
  http://lammps.sandia.gov/doc/Section_start.html#start_3
