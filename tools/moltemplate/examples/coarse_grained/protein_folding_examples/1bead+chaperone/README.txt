# This directory contains examples of how to run a short simulation of a
# coarse-grained protein-like polymer, folding in the presence and absence of
# a chaperone (modeled as an attractive or repulsie spherical shell).
#
# The protein models and the chaperone models are described and used here:
# AI Jewett, A Baumketner and J-E Shea, PNAS, 101 (36), 13192-13197, (2004)
# (http://www.pnas.org/content/101/36/13192)
# ...and also here:
# AI Jewett and J-E Shea, J. Mol. Biol, Vol 363(5), (2006)
#
# (In the "frustrated+minichaperone" directory, the protein is
#  placed outside the chaperone sphere, as opposed to inside.)
#
# -------- REQUIREMENTS: ---------
# 1) These examples require the "USER-MISC" package.  (Use "make yes-USER-MISC")
#   http://lammps.sandia.gov/doc/Section_start.html#start_3
# 2) They also may require additional features and bug fixes for LAMMPS.
#   be sure to download and copy the "additional_lammps_code" from
#   http://moltemplate.org     (upper-left corner menu)
# 3) Unpack it
# 4) copy the .cpp and .h files to the src folding of your lammps installation.
# 5) Compile LAMMPS.
#

-------------
Instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files in each directory.

step 1)
README_setup.sh

step2)
README_run.sh
