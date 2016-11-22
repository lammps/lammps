# -------- REQUIREMENTS: ---------
# 1) This example requires the "MANYBODY" package.
#    As of 2012-9, it is included by default, but this may change in the future.
#    If lammps complains of a missing pair style enter "make yes-MANYBODY"
#    into the shell before compiling lammps.  For details see:
#   http://lammps.sandia.gov/doc/Section_start.html#start_3
This is a relatively complex example containing two different types of
molecules, and a hybrid of Lennard-Jones and 3-body SW "pair" styles.

The cyclododecane molecule uses the 
TraPPE force field for hydrocarbon chains.  
The parameters for the TraPPE force field are
in a file named "trappe1998.lt" which should be
located in the MOLTEMPLATE_PATH.
(See moltemplate installation instructions.)

The water solvent is implemented using the 3-body single-particle 
coarse-grained "mW" water model:
Molinero, V. and Moore, E.B., J. Phys. Chem. B 2009, 113, 4008-4016

More detailed instructions on how to build LAMMPS input files and 
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step 2)
README_run.sh
