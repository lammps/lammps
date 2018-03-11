# --- Running LAMMPS ---
#  -- Prerequisites: --
# The "run.in.nvt" file is a LAMMPS input script containing
# references to the input scripts and data files
# you hopefully have created earlier with moltemplate.sh:
#   system.in.init, system.in.settings, system.data, and table_int.dat
# If not, carry out the instructions in "README_setup.sh".
#
#  -- Instructions: --
# If "lmp_mpi" is the name of the command you use to invoke lammps,
# then you would run lammps on these files this way:


lmp_mpi -i run.in.npt  # Run a simulation at constant pressure (tension)

#or

lmp_mpi -i run.in.nvt  # Run a simulation at constant volume

#(Note: The constant volume simulation lacks pressure equilibration. These are
#       completely separate simulations. The results of the constant pressure
#       simulation are ignored when beginning the simulation at constant volume.
#       This can be fixed.  Read "run.in.nvt" for equilibration instructions.)





# If you have compiled the MPI version of lammps, you can run lammps in parallel
#mpirun -np 4 lmp_mpi -i run.in.npt
#or
#mpirun -np 4 lmp_mpi -i run.in.nvt
# (assuming you have 4 processors available)
