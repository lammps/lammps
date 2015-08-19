# This is a small version of a carbon-nanotube, water capillary system.
# It was inspired by this paper:
#
#    Laurent Joly, J. Chem. Phys. 135(21):214705 (2011)
#
# Note: To investigate the behavior from that paper, you would have to increase 
#       the spacing between the two graphene sheets to prevent the water from
#       making contact with the lower graphene wall. 
#
# Requirements: 1) Set your $MOLTEMPLATE_PATH variable
#               2) The "RIGID" LAMMPS package may be needed later
# To run this system at constant pressure, it might help to compile LAMMPS with
# the optional RIGID package, and use "fix rigid" on the carbon.  (Optional.)
#
# Also, if you have not yet done this set your MOLTEMPLATE_PATH environment 
# variable to access it.  (See installation instructions.)  
# Most likely some of the files in this example (like graphene.lt, tip3p2004.lt)
# are not in this directory, but are in the "common" directory.
# Set MOLTEMPLATE_PATH to point to the "common" directory.
#
# -----------------------------------------------------------
#
# To run moltemplate, use:

moltemplate.sh system.lt

# If you have VMD installed, you can run it this way to visualize the system:
#
# moltemplate.sh -vmd system.lt
#
# Either way, this will generate the following files (which LAMMPS reads):
#    system.data, system.in, system.in.init, system.in.settings
