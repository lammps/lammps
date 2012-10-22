# This is just an example.
#
# Note: The 3 files "run.in.min", "run.in.npt", and "run.in.nvt" are LAMMPS 
#       input scripts which link to the input scripts and data files you 
#       created earlier with moltemplate.sh:
#       system.in.init, system.in.settings, system.data

LAMMPS_COMMAND="lmp_ubuntu"

# Here "$LAMMPS_BINARY" is the name of the command you use to invoke lammps
# (such as lmp_linux, lmp_g++, lmp_mac, lmp_cygwin, etc...) Change if necessary.

# Run lammps using the following 3 commands:

"$LAMMPS_COMMAND" -i run.in.min    # minimize         (OPTIONAL)
"$LAMMPS_COMMAND" -i run.in.npt    # equilibrate the pressure
#"$LAMMPS_COMMAND" -i run.in.nvt    # production run   (OPTIONAL)

# Alternately, if you have MPI installed, try something like this:

#NUMPROCS=4
#mpirun -np $NUMPROCS "$LAMMPS_COMMAND" -i run.in.min # minimize  (OPTIONAL)
#mpirun -np $NUMPROCS "$LAMMPS_COMMAND" -i run.in.npt # equilibrate the pressure
#mpirun -np $NUMPROCS "$LAMMPS_COMMAND" -i run.in.nvt # production run (OPTIONAL)

