This example is a simple simulation of a long alkane chain,
in a vacuum at room temperature using the OPLSAA force field.  


NOTE: This particular file contains instructions for how to build molecules 
      using the OPLSAA force-field.  However, moltemplate is not limited to
      OPLSAA.  Moltemplate allows users to access any of the force-field 
      styles available in LAMMPS (including custom, user-defined force-fields 
      implemented by editing the LAMMPS source code).


-------- INSTRUCTIONS FOR USING OPLSAA WITH YOUR OWN MOLECULES: --------

1) Download the "oplsaa.prm" file containing the OPLSAA force field
parameters.  I do not have permission to distribute this file, 
but you can download the latest version from one of these URLS:
http://dasher.wustl.edu/tinker/distribution/params/oplsaa.prm
http://dasher.wustl.edu/ffe/distribution/params/oplsaa.prm

2) Create the "oplsaa_subset.prm" file by making a copy of the "oplsaa.prm"
file, renaming it to "oplsaa_subset.prm", and deleting the atoms you don't need.
For example, if you are building a simple alkane chain, you would delete every
line beginning with the word "atom", except for these three lines:

atom         80   13    CT    "Alkane CH3-"                  6    12.011    4
atom         81   13    CT    "Alkane -CH2-"                 6    12.011    4
atom         85   46    HC    "Alkane H-C"                   1     1.008    1

(Leave the rest of the file unmodified.)  


3) Create the "oplsaa.lt" file using this command:

oplsaa_moltemplate.py oplsaa_subset.prm

   (Credit to Jasen Lambert for contributing this useful script.)


4) Create the "system.data", "system.in.init", and "system.in.settings"
files which LAMMPS will read by running:

moltemplate.sh system.lt


5) 
To run LAMMPS, you must make sure LAMMPS was built with the "USER-MISC" package.
(because oplsaa_moltemplate.py uses dihedral_style fourier)
To do this, type "make yes-user-misc" before compiling LAMMPS.
 http://lammps.sandia.gov/doc/Section_start.html#start_3


6) Run LAMMPS in this order:

lmp_g++ -i run.in.min
lmp_g++ -i run.in.npt
lmp_g++ -i run.in.nvt

(Replace "lmp_g++" with the name of the LAMMPS executable you are using.)

---- Details ----

The "Alkane50" molecule, as well as the "CH2", and "CH3" monomers it contains
use the OPLSAA force-field.  This means that when we define these molecules, 
we only specify the atom names, bond list, and coordinates.  
We do not have to list the atom charges, angles, dihedrals, or impropers.
The rules for creating atomic charge and angle topology are contained in
the "oplsaa.lt" file created by step 3) above.  The "ch2group.lt", 
"ch3group.lt", and "alkane50.lt" files all refer to "oplsaa.lt", 
(as well as the "OPLSAA" force-field object which it defines).  Excerpt:

import "oplsaa.lt"
CH2 inherits OPLSAA { ... 
CH3 inherits OPLSAA { ... 
Alkane50 inherits OPLSAA { ... 

Alternatively, you can manually define a list of angles, dihedrals, and
improper interactions in these files, instead of asking the force-field
to generate them for you.  You can also specify some of the angles and
dihedrals explicitly, and let the force-field handle the rest.
(Many of the molecule examples which come with moltemplate do this.)
