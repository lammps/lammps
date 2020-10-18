#!/bin/sh

cat > example.py <<EOF
from pylammps import *

lmp = lammps_open_no_mpi(0,None,None)
ver = lammps_version(lmp)

print("LAMMPS version ",ver)
lammps_close(lmp)
EOF

PYTHONPATH=$PWD:${PYTHON_PATH-PWD}

export PYTHON_PATH

python example.py
