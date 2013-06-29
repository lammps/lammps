# set environment for LAMMPS executables to find potential files
if ( "$?LAMMPS_POTENTIALS" == 0 ) setenv LAMMPS_POTENTIALS /usr/share/lammps/potentials
if ( "$?BIOSYM_LIBRARY" == 0 ) setenv BIOSYM_LIBRARY /usr/share/lammps/biosym_frc_files
