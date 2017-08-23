# -------- REQUIREMENTS: ---------
# 1) This example requires the "USER-MISC" package.  (Use "make yes-USER-MISC")
#   http://lammps.sandia.gov/doc/Section_start.html#start_3
# 2) It also may require additional features and bug fixes for LAMMPS.
#   be sure to download and copy the "additional_lammps_code" from
#   http://moltemplate.org     (upper-left corner menu)
# 3) Unpack it
# 4) copy the .cpp and .h files to the src folding of your lammps installation.
# 5) Compile LAMMPS.

This is an example of a very simple coarse-grained protein.

This example contains a 1-bead (C-alpha model) representation of the
"unfrustrated" 4-helix bundle model used in this paper:
G. Bellesia, AI Jewett, and J-E Shea, Protein Science, Vol19 141-154 (2010)

In this model, there are three atom-types (bead-types), H, L, and N
representing one amino-acid each.  The "H" beads represent the hydrophobic
amino acids, and are attracted to eachother with a strength of "1.0"
(in dimensionless units of "epsilon").  The "L" and "N" atoms are
hydrophilic and purely repulsive, and only differ in their secondary-structure
propensity (ie their dihedral parameters).

The dihedral-interaction is bi-stable with two deep local minima (corresponding
to helix-like and sheet-like secondary structure).  You can adjust the bias
in favor of one minima or another by modifying the angle-shift parameter in
the appropriate "dihedral_coeff" command in the other .lt file.

A definition for the 4-sheet beta-barell protein model is also included.
If you want to simulate that molecule instead, then edit the "system.lt"
file (in the "moltemplate_files" subdirectory), and replace this line:
prot = new 4HelixBundle
   with
prot = new 4SheetBundle

-------------
Instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step2)
README_run.sh
