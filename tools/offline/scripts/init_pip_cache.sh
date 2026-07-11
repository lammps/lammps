#!/bin/bash

if [ -z "${LOGGING_DIR}" ]
then
    echo "Must set LOGGING_DIR environment variable"
    exit 1
fi

if [ -z "${PIP_CACHE_DIR}" ]
then
    echo "Must set PIP_CACHE_DIR environment variable"
    exit 1
fi

set -x

mkdir -p "$PIP_CACHE_DIR"

# download packages that might be needed
cd "$PIP_CACHE_DIR"
pip3 download pip setuptools wheel
pip3 download -r "$LAMMPS_DIR/doc/utils/requirements.txt"
