# --- Running LAMMPS ---
#  -- Prerequisites: --
# The "run.in.nvt" LAMMPS input script links to the input 
# scripts and data files you hopefully have created earlier 
# with moltemplate.sh:
#   system.in.init, system.in.settings, system.data
# If not, carry out the instructions in "README_run_moltemplate.sh".
#
#  -- Instructions: --
# If "lmp_linux" is the name of the command you use to invoke lammps,
# then you would run lammps this way:

lmp_linux -i run.in.nvt

# NOTE: BECAUSE ALL OF THE ORIGINAL FORCE FIELD PARAMETERS WERE INTENTIONALLY
#       REMOVED, THE SYSTEM WILL MOVE IN A VERY UNREALISTIC WAY.
