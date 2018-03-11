This directory contains many .LT files molecules
which have been parameterized with the MARTINI force field
and converted into moltemplate format.

------- CITE -----------------------------
NOTE: We extracted the parameters in the MARTINI force field from the files
distributed with the "EMC" tool.  If you use these .lt files, please also cite:
P. J. in ‘t Veld and G. C. Rutledge, Macromolecules 2003, 36, 7358.
---------------------------------------------
THESE FILES HAVE NOT BEEN CAREFULLY TESTED.  PLEASE REPORT BROKEN EXAMPLES TO:
jewett.aij -at- gmail.com, or report issues at github.com/jewettaij/moltemplate
---------------------------------------------
PLEASE SHARE ANY NEW EXAMPLES YOU CREATE WITH THE COMMUNITY, either by emailing:
jewett.aij -at- gmail.com, or a pull request at github.com/jewettaij/moltemplate
---------------------------------------------

How to use these files:

Look at the "DOPC_bilayer_preformed" example.
In particular, look at the "moltemplate_files/system.lt" file.

This example contains only one kind of lipid, but you can create a mixture
by replacing the "lipids = new DPPC" command with a command like:
lipids = new random([DPPC, DSPC], [195,195], 1234) 

----- comments -----

Several of the examples in the "MARTINI_examples" directory are limited to
only one kind of lipid.  In these examples, the force field parameters were
hard coded inside the definition of the lipid molecule
(specifically, the DPPC.lt file).
This makes the examples slightly easier to read and understand because
everything is contained in the same file, but not useful for general use.

This directory, on the other hand, contains more general .LT files useful
for building multi-component membranes with multiple types of lipids and
sterols.  Because most of these molecules share many of the the same atom
types and force field parameters, all of this information has been saved
in a separate file, "martini.lt", which is located in the
"moltemplate/force_fields" subdirectory (distributed with moltemplate).

The conversion to MOLTEMPLATE (.LT) format was done by
David Stelter and Saeed Momeni Bashusqeh, converting the
EMC files distributed by Pieter J. in 't Veld.

--- Generalizing to other Lipids ---

More Lipids and Sterols can be downloaded at:
http://cgmartini.nl/index.php/force-field-parameters/lipids
(http://cgmartini.nl)
in gromacs .itp and .GRO formats and converted into moltemplate format.

For inspiration how this should be done, download the DPPC molecule files from
http://cgmartini.nl/index.php/force-field-parameters/lipids2/351-lipid.html?dir=PC&lipid=DOPC
and compare these files with the DOPC.lt file in this directory.

We copied the coordinates from the DOPC.gro file into the DOPC.lt file,
and (attempted to) make sure that the atom types matched the atom types
in the "martini.lt" file (which should also match the names in the .ITP files).

--- Generalizing to other kinds of molecules (eg. amino acids ---

The "martini.lt" file currently only contains the definitions of atoms
used to simulate lipids and sterols.
To simulate other molecules such as proteins using the MARTINI force field,
you will need to create a more comprehensive "martini.lt" file which includes
these new particle types.  One way to do that is to download the .PRM files
(EMC format) for the molecule types you are interested in, and include them
in the list of .PRM files when you run the "emcprm2lt.py" conversion script.

For more details how this can be done, go to the "force_field" subdirectory
and look for the "martini.lt" file.  Additional .PRM files are located in
the "martini" subdirectory in that folder.

