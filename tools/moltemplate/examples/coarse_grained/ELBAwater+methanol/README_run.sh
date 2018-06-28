# --- Running LAMMPS ---
#
# The file "run.in.npt" is a LAMMPS
# input script which link to the input scripts and data files
# you hopefully have created earlier with moltemplate.sh:
#   system.in.init, system.in.settings, system.data
# If not, carry out the instructions in "README_setup.sh".
#
#  -- Instructions: --
# If "lmp_mpi" is the name of the command you use to invoke lammps,
# then you would run lammps on these files this way:

lmp_mpi -i run.in.npt  # simulation at constant pressure

# (Note: The "lmp_mpi" program is also frequently called "lmp_ubuntu", 
#        and has many other names depending upon how it was compiled.)

# If you have compiled the MPI version of lammps, you can run lammps in parallel
#mpirun -np 4 lmp_mpi -i run.in.npt
# (assuming you have 4 processors available)
