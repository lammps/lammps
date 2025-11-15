#!/bin/sh
# wrapper for bundled executables

# reset locale to avoid problems with decimal numbers
export LC_ALL=C

BASEDIR="$(dirname "$0")"
EXENAME="$(basename "$0")"

# save old settings (for restoring them later)
OLDPATH="${PATH}"
OLDLDLIB="${LD_LIBRARY_PATH}"

# prepend path to find our custom executables
PATH="${BASEDIR}/bin:${PATH}"

# append to LD_LIBRARY_PATH to prefer local (newer) libs
# however the LAMMPS shared library must always come first
# so that it does not get overridden by other installed versions
LD_LIBRARY_PATH="${BASEDIR}/libexec/lammps:${LD_LIBRARY_PATH}:${BASEDIR}/lib"

# set some environment variables for LAMMPS etc.
LAMMPS_POTENTIALS="${BASEDIR}/share/lammps/potentials"
MSI2LMP_LIBRARY="${BASEDIR}/share/lammps/frc_files"

# export everything
export LD_LIBRARY_PATH LAMMPS_POTENTIALS MSI2LMP_LIBRARY PATH OLDPATH OLDLDLIB

exec "${BASEDIR}/bin/${EXENAME}" "$@"
