# The version of "2bead.lt" in this directory is defines two types of
# molecules (named "2bead/H" and "2bead/P").
#
# However, there is another version of this file which is easier to understand.
# I recommend reading that file first.
# It is located at "simplified_version_one_residue/2bead.lt".
# It defines only one type of molecule (named "2bead")
# It is much simpler.

# ------

# Use these commands to generate the LAMMPS input script and data file
# (and other auxilliary files):

# run moltemplate

moltemplate.sh system.lt

# This will generate various files with names ending in *.in* and *.data. 

