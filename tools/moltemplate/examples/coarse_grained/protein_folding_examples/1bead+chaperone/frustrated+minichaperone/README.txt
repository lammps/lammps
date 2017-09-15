# This directory demonstrates how to run a long simulation of
# the "frustrated" coarse-grained protein in the presence of one
# or more coarse-graine small ("mini") chaperones (R=3, h=0.6) as described in:
#
# AI Jewett and J-E Shea, J. Mol. Biol, Vol 363(5), (2006)
#   and earlier in:
# AI Jewett, A Baumketner and J-E Shea, PNAS, 101 (36), 13192-13197, (2004)
#
# Because this process takes a long time (even with the help of the chaperone)
# I save the data relatively infrequently.
#
# Note: In the 2006 paper, only one protein and one chaperone was simulated.
#       In this example, 8 proteins and 8 chaperones were simulated.
#
# Note: In this case, the chaperones appear to catalyze aggregation.
#       This is due to an artifact in the protein model.  That model
#       was not designed to study aggregation.  However the simulation
#       is suitable for making pretty pictures (to show off moltemplate).
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
