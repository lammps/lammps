This directory contains instructions for creating a a moltemplate file
("oplsaa.lt") containing force-field definitions relevant to the "Alkane50"
molecule.  (However, these instructions should work for other molecules too.)

--- Instructions ---

First, check and see if there is an "oplsaa_subset.prm" file present.
If not, then download this file:

http://dasher.wustl.edu/tinker/distribution/params/oplsaa.prm
   This file is also available here:
http://dasher.wustl.edu/ffe/distribution/params/oplsaa.prm

and save this file as "oplsaa_subset.prm".  Then you must EDIT THIS FILE
so that it only contains atom types you plan to have in your simulation
(see below for details).  Then run the opls_moltemplate.py script this way:


python oplsaa_moltemplate.py  oplsaa_subset.prm


This will create a file named "oplsaa.lt"
Look over the newly created "oplsaa.lt" file.
Then move this file to wherever you plan to run moltemplate.  For example:

mv -f oplsaa.lt ..

----- DETAILS: Editing the "oplsaa_subset.prm" file -------

Again, before you run "oplsaa_moltemplate.py", you must edit the "oplsaa.prm"
file (or "oplsaa_subset.prm file) and eliminate atom types which do not 
correspond to any of the atoms in your simulation.  This means you must
look for lines near the beginning of this file which begin with the word "atom"
and refer to atom types which appear in the simulation you plan to run.  All
other lines (beginning with the word "atom") must be deleted or commented out.
(Leave the rest of the file unmodified!)

For example:
If you were working with a simple alkane chain, you would delete every line
beginning with the word "atom", except for these three lines:


atom         80   13    CT    "Alkane CH3-"                  6    12.011    4
atom         81   13    CT    "Alkane -CH2-"                 6    12.011    4
atom         85   46    HC    "Alkane H-C"                   1     1.008    1

Then you are ready to run oplsaa_moltemplate.py on this file.

(Note: The atom type numbers, like "89", "81", "85", "13", "46", etc... may vary
       depending on when you downloaded the "oplsaa.prm" file.  Be sure to check
       the descriptions of each atom type after you download it: "Alkane CH3-")


----- Using the "oplsaa.lt" file -----

Once you have created the "oplsaa.lt" file, you can create files (like 
ethylene.lt) which define molecules that refer to these atom types.
Here is an excerpt from "methane.lt":

import "oplsaa.lt"
Methane inherits OPLSAA {
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
the TINKER program.  (Because these examples use the "oplsaa.prm" file which
is distributed with TINKER.)  I think these are the relevant citations:

1) Ponder, J. W., & Richards, F. M. (1987). "An efficient newton‐like method for molecular mechanics energy minimization of large molecules. Journal of Computational Chemistry", 8(7), 1016-1024.

2) Ponder, J. W, (2004) "TINKER: Software tools for molecular design", http://dasher.wustl.edu/tinker/

-------------------------------

Andrew Jewett and Jason Lambert
May, 2014

Please email bugs to jewett.aij@gmail.com and jlamber9@gmail.com
