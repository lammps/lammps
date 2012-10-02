# You would probably run lammps this way:
#
# lmp_ubuntu -i run.in.nvt

# The "run.in.nvt" file is a LAMMPS input scripts which refers
# to the input scripts & data files you created earlier when you ran moltemplate
# system.in.init, system.in.settings, system.data



# -----------------------------------



LAMMPS_COMMAND="lmp_ubuntu"

# Here "$LAMMPS_BINARY" is the name of the command you use to invoke lammps
# (such as lmp_linux, lmp_g++, lmp_mac, lmp_cygwin, etc...) Change if necessary.

# Run lammps using the following 3 commands:

"$LAMMPS_COMMAND" -i run.in.min    # minimize
"$LAMMPS_COMMAND" -i run.in.npt    # pressure equilibration
"$LAMMPS_COMMAND" -i run.in.nvt    # production run

# Alternately, if you have MPI installed, try something like this:

#NUMPROCS=4
#mpirun -np $NUMPROCS "$LAMMPS_COMMAND" -i run.in.min
#mpirun -np $NUMPROCS "$LAMMPS_COMMAND" -i run.in.nvt
#mpirun -np $NUMPROCS "$LAMMPS_COMMAND" -i run.in.nvt
