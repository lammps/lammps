This example is a simple simulation of 288 hexadecane molecules in a box at
room temperature and atmospheric pressure.  Please read the WARNING.TXT file.

-------- REQUIREMENTS: ---------
This example requires building LAMMPS with the "USER-MISC" package.
(because it uses dihedral_style fourier)
To do this, type "make yes-user-misc" before compiling LAMMPS.
 http://lammps.sandia.gov/doc/Section_start.html#start_3

More detailed instructions on how to build LAMMPS input files and 
run a short simulation are provided in other README files:

step 1) to setup the LAMMPS input files, run this file:
README_setup.sh

step 2) to run LAMMPS, follow the instructions in this file:
README_run.sh

------------ NOTE: There are two versions of this example. ----------------

Both examples use the same force-field parameters.

1) In this directory, all of the force-field parameters are listed explicitly
in the "alkanes.lt" file (located in the "moltemplate_files" directory).
This allows the user to manually control all of the force-field details.

2) However, there is an alternate version of this example in the 
"../AMBER_force_field_examples" directory.
In that version, the force-fields are loaded from a much larger file named
"gaff.lt" which contains all of the parameters in the AMBER GAFF force-field
database.  The "gaff.lt" is similar to the "alkanes.lt" file except that
it is larger (because it contains information for nearly all small organic
molecules).  It is located in a different directory (in the "common" directory).

Relying on "gaff.lt" frees the user from the drudgery of manually specifying
all of these force-field details for every molecule.  (However, the user must
be careful to choose @atom-type names which match AMBER GAFF conventions,
such as "c3" and "h1", in this example.)
