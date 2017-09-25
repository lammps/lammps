# --- Running LAMMPS ---
#  -- Prerequisites: --
# The "run.in.nvt" file is a LAMMPS input script containing
# references to the input scripts and data files
# you hopefully have created earlier with MOLTEMPLATE and PACKMOL:
#   system.in.init, system.in.settings, system.in.coords, system.data,
#   and table_int.dat
# If not, carry out the instructions in "README_setup.sh".
#
#  -- Instructions: --
# If "lmp_mpi" is the name of the command you use to invoke lammps,
# then you would run lammps on these files this way:

lmp_mpi -i run.in.min  # Minimize the system (important, and very slow)

lmp_mpi -i run.in.nvt  # Run a simulation at constant volume



# If you have compiled the MPI version of lammps, you can run lammps in parallel
#mpirun -np 4 lmp_mpi -i run.in.min
#or
#mpirun -np 4 lmp_mpi -i run.in.nvt
# (assuming you have 4 processors available)
