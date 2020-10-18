#!/bin/sh

cat > example.py <<EOF
from pylammps import *

lmp = lammps_open_no_mpi(0,None,None)
ver = lammps_version(lmp)
lammps_command(lmp, "units real")
lammps_command(lmp, "lattice fcc 2.5")
lammps_command(lmp, "region box block -5 5 -5 5 -5 5")
lammps_command(lmp, "create_box 1 box")
lammps_command(lmp, "create_atoms 1 box")

print("LAMMPS version ",ver)
print("Number of created atoms: %g" % lammps_get_natoms(lmp))
lammps_close(lmp)
EOF

PYTHONPATH=$PWD:${PYTHON_PATH-PWD}

export PYTHON_PATH

python example.py
