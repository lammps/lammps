# This directory demonstrates how to run a long simulation of
# the "frustrated" coarse-grained protein confined in a frustrated
# coarse-grained chaperonin (R=6, h=0.475) as described in:
# AI Jewett, A Baumketner and J-E Shea, PNAS, 101 (36), 13192-13197, (2004)
# (http://www.pnas.org/content/101/36/13192)
#
# Note: If you want to use a "hydrophilic" chaperone (with h=0.0
#       instead of h=0.475), then replace the word "CHAP_INTERIOR_H0.475"
#       (at the end of "system.lt") with "CHAP_INTERIOR_H0"
#
# Because this process takes a long time (even with the help of the chaperone)
# I save the data relatively infrequently.
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
