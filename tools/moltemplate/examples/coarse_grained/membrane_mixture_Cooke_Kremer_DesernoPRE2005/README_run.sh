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


# If you have compiled the MPI version of lammps, you can run lammps in parallel
#mpiexec -np 4 lmp_mpi -i run.in.npt
#or
#mpiexec -np 4 lmp_mpi -i run.in.nvt
# (assuming you have 4 processors available)
