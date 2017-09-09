# This example shows how to build a multicomponent spherical vesicle.
# The lipid bilayer is composed of two different lipids (DPPC and DLPC),
# The vesicle contains 120 trans-membrane protein inclusions.
#
# ---------------- Prerequisites: ------------------
# You must run packmol to generate the coordinates beforehand.
# Afterwards, move and rename the final coordinate file to "../system.xyz"
# To do this, check the README.sh file in the ../packmol_files directory.
# (or follow these instructions below)
#
#  cd ../packmol_files
#  packmol < step1_proteins.inp
#  packmol < step2_innerlayer.inp
#  packmol < step3_outerlayer.inp
#  cp step3_outerlayer.xyz ../system.xyz
#
#  These steps could take a few hours.
#
# --- After you have done that, you can run moltemplate using this command: ---

moltemplate.sh system.lt -xyz ../system.xyz

