This directory contains instructions for creating a a moltemplate file
("oplsaa.lt") containing force-field definitions relevant to ethylene example.
Create "oplsaa.lt" this way:

./oplsaa_moltemplate.py  oplsaa_subset.txt

Then move the newly created opls.lt file to the directory
where you plan to run moltemplate.

There is a directory containing an example how to use this code located here:
examples/all_atom_examples/OPLSAA_force_field_examples/ethylene/

----------- GENERAL DOCUMENTATION -----------
The "oplsaa_subset.txt" file should contain force-field information
for a subset of the atoms in the oplsaa forcefield.  This directory
also contains the python code ("oplsaa_moltemplate.py") used to
generate that file, as well as the full oplsaa force field file: "oplsaa.txt".

NOTE: THE ORIGINAL "oplsaa.txt" FILE IS LOCATED IN THE "common/opls/" DIRECTORY.
      (I deleted it here to save space.)

python oplsaa_moltemplate.py  <force field file name>

First, back-up the oplsaa.txt file and then delete all the unnecessary atoms
from the file.

It is ideal for you to run this on a subset of the OPLS forcefield relevant 
to your problem.  (It is possible for you to generate a full OPLS force field
moltemplate file but that demands a lot of time and generates a file that is
on the order 100 gigabytes, which moltemplate can not read.)  
Save yourself time and energy, and make a copy of the oplsaa.txt.  Then edit 
this file and delete or comment out any line beginning with the word "atom"
which is not relevant to your problem.  (However keep the rest of the file.)

For example:
Let's assume this new file is called "oplsaa_subset.txt".  If you were 
working with ethylene you would delete every atom except for atom 85(C=) 
and 86(H-C=) and use that as the input file for oplsaa_moltemplate.py

atom         85   37    CM    "Alkene H2-C="                 6    12.011    3
atom         86   36    HC    "Alkene H-C="                  1     1.008    1

You do not need to make any other changes other than that.  Then run:

python oplsaa_moltemplate.py oplsaa_subset.txt

   ...this will create a file named "oplsaa.lt"

Look over the newly created "oplsaa.lt" file.  Then copy it to the "common" directory (where other force-field files like gaff.lt are located, ...or just copy it to wherever you plan to run moltemplate).  

Then you can create molecules that refer to it using:

Ethylene inherits OPLSAA {
  # Definition of Ethylene molecule goes here.
}

CHARGE:

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

-Jason Lambert
February, 2014

Please email bugs to jewett.aij@gmail.com

