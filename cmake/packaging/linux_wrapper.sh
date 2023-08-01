#!/bin/sh
# wrapper for bundled executables

BASEDIR=$(dirname "$0")
EXENAME=$(basename "$0")

# append to LD_LIBRARY_PATH to prefer local (newer) libs
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BASEDIR}/lib
export LD_LIBRARY_PATH

exec "${BASEDIR}/bin/${EXENAME}" "$@"
