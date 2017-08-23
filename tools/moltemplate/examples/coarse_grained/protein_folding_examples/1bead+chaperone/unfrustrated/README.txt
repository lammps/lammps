# This directory demonstrates how to run a short simulation of
# the "unfrustrated" coarse-grained protein model used in:
# AI Jewett, A Baumketner and J-E Shea, PNAS, 101 (36), 13192-13197, (2004)
# (http://www.pnas.org/content/101/36/13192)
#
# During this short simulation (run.in.nvt) the protein evolves
# from an unfolded initial conformation to the folded state.
#
# -------- REQUIREMENTS: ---------
# 1) This example requires the "USER-MISC" package.  (Use "make yes-USER-MISC")
#   http://lammps.sandia.gov/doc/Section_start.html#start_3
# 2) It also may require additional features and bug fixes for LAMMPS.
#   be sure to download and copy the "additional_lammps_code" from
#   http://moltemplate.org     (upper-left corner menu)
# 3) Unpack it
# 4) copy the .cpp and .h files to the src folding of your lammps installation.
# 5) Compile LAMMPS.

-------------
Instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step2)
README_run.sh
