This example contains a mixture of water(SPCE) and methane.
The methane molecules use OPLSAA force-field, but the water molecules do not.

---- Details ----

The methane molecules in this example use the OPLSAA force-field.
This means that the database of force-field parameters in "oplsaa.lt"
will be used to generate angles, dihedrals, and impropers.
The "moltemplate_files/methane.lt" file
contains these lines which refer to OPLSAA:

import "oplsaa.lt"
Methane inherits OPLSAA { ...

However the "SPCE" (water) molecules does NOT use a database to look up the
force-field parameters for this tiny molecule.
Instead, the "moltemplate_files/spce.lt" file declares all of the angle
interactions, atom properties and force-field parameters for water explicitly.
(Consequently, it makes no mention of "oplsaa.lt" or "OPLSAA".)

-------- Instructions: ---------

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step 2)
README_run.sh
