This example demonstrates how to build a simulation containing a box of methane.
(Not a very interesting example.)

---- Details ----

The methane molecules in this example use the OPLSAA force-field.
This means that the database of force-field parameters in "oplsaa.lt"
will be used to generate angles, dihedrals, and impropers.
The "moltemplate_files/methane.lt" file
contains these lines which refer to OPLSAA:

import "oplsaa.lt"
Methane inherits OPLSAA { ...

-------- Instructions: ---------

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step 2)
README_run.sh
