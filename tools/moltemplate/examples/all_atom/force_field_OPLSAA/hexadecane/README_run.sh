# --- Running LAMMPS ---
#
# The 2 files "run.in.npt", and "run.in.nvt" are LAMMPS
# input scripts which link to the input scripts and data files
# you hopefully have created earlier with moltemplate.sh:
#   system.in.init, system.in.settings, system.data
# If not, carry out the instructions in "README_setup.sh".
#
#  -- Instructions: --
# If "lmp_mpi" is the name of the command you use to invoke lammps,
# then you would run lammps on these files this way:


lmp_mpi -i run.in.npt  # minimization and simulation at constant pressure
lmp_mpi -i run.in.nvt  # simulation at constant volume


# If you have compiled the MPI version of lammps, you can run lammps in parallel
#mpirun -np 4 lmp_mpi -i run.in.npt
#mpirun -np 4 lmp_mpi -i run.in.nvt
# (assuming you have 4 processors available)
