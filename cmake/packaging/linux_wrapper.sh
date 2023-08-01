#!/bin/sh
# wrapper for bundled executables

BASEDIR=$(dirname "$0")
EXENAME=$(basename "$0")

# append to LD_LIBRARY_PATH to prefer local (newer) libs
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BASEDIR}/lib

# set some environment variables for LAMMPS etc.
LAMMPS_POTENTIALS=${BASEDIR}/share/lammps/potentials
MSI2LMP_LIBRARY=${BASEDIR}/share/lammps/frc_files
export LD_LIBRARY_PATH LAMMPS_POTENTIALS MSI2LMP_LIBRARY

exec "${BASEDIR}/bin/${EXENAME}" "$@"
