# --- Running LAMMPS ---
#  -- Prerequisites: --
# The "run.in.nvt" file is a LAMMPS input script containing
# references to the input scripts and data files
# you hopefully have created earlier with moltemplate.sh:
#   system.in.init, system.in.settings, system.data
# If not, carry out the instructions in "README_setup.sh".
#
#  -- Instructions: --
# If "lmp_linux" is the name of the command you use to invoke lammps,
# then you would run lammps on these files this way:


lmp_linux -i run.in.nvt  # Run a simulation at constant volume

#or 

lmp_linux -i run.in.npt  # Run a simulation at constant pressure
                         # (Note: Constant pressure conditions have not been 
                         #  well tested.  The "run.in.npt" script may fail.)



# If you have compiled the MPI version of lammps, you can run lammps in parallel
#mpirun -np 4 lmp_linux -i run.in.nvt
#or
#mpirun -np 4 lmp_linux -i run.in.npt
# (assuming you have 4 processors available)
