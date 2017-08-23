This example is a simple simulation of a long alkane chain,
in a vacuum at room temperature using the OPLSAA force field.

1) Create the "system.data", "system.in.init", and "system.in.settings"
files which LAMMPS will read by running:

moltemplate.sh system.lt


2) Run LAMMPS in this order:

lmp_mpi -i run.in.min   # minimize the energy (to avoid atom overlap) before...
lmp_mpi -i run.in.nvt   # running the simulation at constant temperature

(The name of the LAMMPS executable, eg "lmp_mpi", may vary.)

---- Details ----

The "Alkane50" molecule, as well as the "CH2", and "CH3" monomers it contains
use the OPLSAA force-field.  This means that when we define these molecules,
we only specify the atom names, bond list, and coordinates.
We do not have to list the atom charges, angles, dihedrals, or impropers.
The rules for creating atomic charge and angle topology are contained in
the "loplsaa.lt" file created by step 3) above.  The "ch2group.lt",
"ch3group.lt", and "alkane50.lt" files all refer to "loplsaa.lt",
(as well as the "OPLSAA" force-field object which it defines).  Excerpt:

import "loplsaa.lt"
CH2 inherits OPLSAA { ...
CH3 inherits OPLSAA { ...
Alkane50 inherits OPLSAA { ...

Alternatively, you can manually define a list of angles, dihedrals, and
improper interactions in these files, instead of asking the force-field
to generate them for you.  You can also specify some of the angles and
dihedrals explicitly, and let the force-field handle the rest.
(Many of the examples which come with moltemplate do this.)
