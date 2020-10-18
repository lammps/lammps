#!/bin/sh

cat > example.lua <<EOF
require("lualammps")

lmp = lualammps.lammps_open_no_mpi(0,nil,nil)
ver = lualammps.lammps_version(lmp)

print("LAMMPS version ", ver)
lualammps.lammps_close(lmp )
EOF

lua example.lua
