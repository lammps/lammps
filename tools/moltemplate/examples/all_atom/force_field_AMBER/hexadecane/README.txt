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

      (Currently there is a bug which makes this step slow.
       I'll fix it later -Andrew 2013-10-15.)

step 2) to run LAMMPS, follow the instructions in this file:
README_run.sh

------------ NOTE: There are two versions of this example. ----------------

Both examples use the same force-field parameters.

1)
In this version, the force-field parameters are loaded from the "gaff.lt" file
(located in the "force_fields" subdirectory of the moltemplate distribution).
This frees the user from the drudgery of manually specifying all of these
force-field details for every molecule.  (However, the user must be careful
to choose @atom-type names which match AMBER GAFF conventions,
such as the "c3" and "h1" atoms, in this example.)

2)
Alternately, there is another "hexadecane" example in the "all_atom_examples"
directory.  In that example, force-field parameters are loaded from a file
named "alkanes.lt" (instead of "gaff.lt").  The "alkanes.lt" file contains
only the excerpts from "gaff.lt" which are relevant to the hydrocarbon
molcules used in that example.  ("gaff.lt" contains parameters for most
small organic molecules, not just hydrocarbons.)
In this way, by editing "alkanes.lt", the user can manually control all of the
force-field details in the simulation.  (Without feeling as though they are
relying on some kind of mysterious "black box" to do it for them.)

