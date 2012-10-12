# This is just an example.
#
# Note: The "run.in.nvt" file is a LAMMPS input script which attempts to read 
#       the input scripts and data files you created with moltemplate:
#       system.in.init, system.in.settings, system.data

LAMMPS_COMMAND="lmp_ubuntu"

# Here "$LAMMPS_BINARY" is the name of the command you use to invoke lammps
# (such as lmp_linux, lmp_g++, lmp_mac, lmp_cygwin, etc...) Change if necessary.

# Run lammps using the following 3 commands:

"$LAMMPS_COMMAND" -i run.in.nvt

# Alternately, if you have MPI installed, try something like this:

#NUMPROCS=4
#mpirun -np $NUMPROCS "$LAMMPS_COMMAND" -i run.in.nvt
