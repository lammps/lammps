# run moltemplate this way

moltemplate.sh system.lt

# This will generate various files with names ending in *.in* and *.data
# which are needed by LAMMPS.

#  ------ Other versions: --------
#
# If you are using the "other_versions/charmm/1beadProtSci2010.lt" file,
# then you must run moltemplate this way:
#
# moltemplate.sh -overlay-dihdedrals system.lt
