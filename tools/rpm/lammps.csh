# set environment for LAMMPS executables to find potential files
if ( "$?LAMMPS_POTENTIALS" == 0 ) setenv LAMMPS_POTENTIALS @DATADIR@/potentials
if ( "$?MSI2LMP_LIBRARY" == 0 ) setenv MSI2LMP_LIBRARY @DATADIR@/frc_files
