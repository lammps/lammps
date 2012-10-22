
LAMMPS_BINARY="lmp_ubuntu"  # change this if your binary has a different name

# To verify that the files created by ltemplify.py and moltemplate.sh are
# valid LAMMPS files, you can start a short simulation with them this way:
$LAMMPS_BINARY -i run.in.nvt

# NOTE: BECAUSE ALL OF THE ORIGINAL FORCE FIELD PARAMETERS WERE INTENTIONALLY
#       REMOVED, THE SYSTEM WILL MOVE IN A VERY UNREALISTIC WAY.
