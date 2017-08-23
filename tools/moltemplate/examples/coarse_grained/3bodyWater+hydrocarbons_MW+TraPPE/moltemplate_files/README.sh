# Use this command to generate the LAMMPS input files:

moltemplate.sh -a "@atom:/WatMW/mW 1" system.lt

# The -a argument insures that the "mW" atom type is assigned to "1".
# (This is necessary for the pair_coeff command to work.
#  See system.lt for details.)

# Note: To get rid of the annoying "atom_style unspecified warnings,
# use the "-atomstyle" command line argument, as in:
# moltemplate.sh -atomstyle full  -a "@atom:/WatMW/mW 1" system.lt
