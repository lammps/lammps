Unfortunately, moltemplate does not come with a file containing OPLSAA 
paramters which is ready to use.  You must build it yourself.
This directory has tools and instructions to explain how to do this.

-----------------------------

When you want to run a new simulation, you must download the full 
"oplsaa.prm" force-field file (from the TINKER web site) and manually 
delete all the atom types which you do not need.  (See below.)
Then you must use the "oplsaa_moltemplate.py" script to create 
"oplsaa.lt" file which moltemplate.sh needs.  Then you must run moltemplate.sh.
You must do this for each new simulation you plan to run which uses OPLSAA.

---- Details: ----

   Download the original "oplsaa.prm" file here:
http://dasher.wustl.edu/tinker/distribution/params/oplsaa.prm
   or here:
http://dasher.wustl.edu/ffe/distribution/params/oplsaa.prm

and save this file as "oplsaa_subset.prm".  Then you must EDIT THIS FILE
so that it only contains atom types you plan to have in your simulation
(see below).  Finally, you must run the opls_moltemplate.py script this way:


python oplsaa_moltemplate.py  oplsaa_subset.prm


This will create a file named "oplsaa.lt"
Look over the newly created "oplsaa.lt" file.
Then, if necessary, move this file to wherever you plan to run moltemplate.

There is a directory containing an example how to use this code located here:
examples/all_atom_examples/OPLSAA_force_field_examples/ethylene/moltemplate_files/oplsaa_lt_generator/README.TXT

----- DETAILS: Editing the "oplsaa_subset.prm" file -------

Again, before you run "oplsaa_moltemplate.py", you must edit the "oplsaa.prm"
file (or "oplsaa_subset.prm file) and eliminate atom types which do not 
correspond to any of the atoms in your simulation.  This means you must
look for lines near the beginning of this file which begin with the word "atom"
and refer to atom types which appear in the simulation you plan to run.  All
other lines (beginning with the word "atom") must be deleted or commented out.
(Leave the rest of the file alone.)

For example:
If you were working with ethylene and benzene you would delete every line
beginning with the word "atom", except for these four lines:

# for ethylene:
atom         88   47    CM    "Alkene H2-C="                 6    12.011    3
atom         89   46    HC    "Alkene H-C="                  1     1.008    1
# for benzene:
atom         90   48    CA    "Aromatic C"                   6    12.011    3
atom         91   49    HA    "Aromatic H-C"                 1     1.008    1

Then you are ready to run oplsaa_moltemplate.py on this file.

(You can try to delete more irrelevant information, but be careful.
 It is not necessary, and it is easy to make mistakes.)


----- Using the "oplsaa.lt" file -----

Once you have created the "oplsaa.lt" file, you can create files (like 
ethylene.lt) which define molecules that refer to these atom types.
Here is an excerpt from "ethyelene.lt":

Ethylene inherits OPLSAA {
  write('Data Atoms') { 
    list of atoms goes here ...
  }
  write('Data Bond List') { 
    list of bonds goes here ...
  }
}

And then run moltemplate.


-----------  CHARGE: -----------

By default, the OPLSAA force-field assigns atom charge according to atom type.
When you run moltemplate, it will create a file named "system.in.charges",
containing commands like:

set type 2 charge -0.42
set type 3 charge 0.21

(This assumes your main moltemplate file is named "system.lt".  If it was
named something else, eg "polymer.lt", then the file created by moltemplate
will be named "polymer.in.charges".)

Include these commands somewhere in your LAMMPS input script 
(or use the LAMMPS "include" command to load the commands in system.in.charges)

Note that the atom numbers (eg "2", "3") in this file will not match the 
OPLS atom numbers.  (Check the output_ttree/ttree_assignments.txt file,
created by moltemplate, to see a table of "@atom" type numbers translated 
from OPLSAA into LAMMPS.)

----------- CREDIT -----------

If you use these tools and you publish a paper using OPLSAA, please also cite 
the TINKER program.  (Because the original "oplsaa.prm" file which we depend on
is distributed with TINKER.)  I think these are the relevant citations:

1) Ponder, J. W., & Richards, F. M. (1987). "An efficient newton‐like method for molecular mechanics energy minimization of large molecules. Journal of Computational Chemistry", 8(7), 1016-1024.

2) Ponder, J. W, (2004) "TINKER: Software tools for molecular design", http://dasher.wustl.edu/tinker/

-------------------------------

Jason Lambert and Andrew Jewett
April, 2014

Please email bugs to jewett.aij@gmail.com and jlamber9@gmail.com

