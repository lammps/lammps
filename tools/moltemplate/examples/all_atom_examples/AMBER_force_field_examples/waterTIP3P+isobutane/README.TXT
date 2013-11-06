This is an example of how to use "canned" force-fields in like GAFF in LAMMPS.
In this example, we specify only the atom names, bond connectivity,
(and coordinates and charge), and use moltemplate to
load the GAFF parameters from an external file (gaff.lt)
(...instead of specifying them explicitly in the molecule definition).

The simulation consists of a mixture of isobutane and water.
Over time (less than 1 ns), the two molecules phase-separate.

The GAFF parameters are applied only to the isobutane molecule.
(The water molecule paramters are defined explicitly in common/tip3p_2004.lt)
For this to work, make sure you have defined the MOLTEMPLATE_PATH
environment variable and set it to "common".  See manual for more details.)

# WARNING: THIS IS A PRELIMINARY EXAMPLE WHICH USES AMBER'S GAFF FORCE FIELD.
#          THIS FEATURE IS CURRENTLY BEING TESTED (AS OF 2013-8-08).
#          THE ABILITY TO DETECT AND ASSIGN GAFF FORCE FIELD PARAMETERS
#          MOLECULES ACCORDING TO ATOM TYPE IS AN EXPERIMENTAL FEATURE
#          AND CURRENTLY PROBABLY HAS BUGS (IN THE DIHEDRALS AND IMPROPERS).
#
#          (In addition, I am embarassed to admit I do not understand
#           atom nomenclature, and I am not sure if I am using
#           the correct GAFF atom names in the isobutane molecule.)
#          PLEASE REPORT BUGS AND/OR SEND CORRECTIONS.  -A 2013-8-08
#
# -------- REQUIREMENTS: ---------
# 1) This example requires building LAMMPS with the "USER-MISC" package.
#    (because it makes use of "gaff.lt" which uses dihedral_style fourier)
#    To do this, type "make yes-user-misc" before compiling LAMMPS.
#   http://lammps.sandia.gov/doc/Section_start.html#start_3
# 2) You must define your MOLTEMPLATE_PATH environment variable
#    and set it to the "common" subdirectory of your moltemplate distribution.
#    (See the "Installation" section in the moltemplate manual.)

More detailed instructions on how to build LAMMPS input files and 
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step 2)
README_run.sh
