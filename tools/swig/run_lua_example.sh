#!/bin/sh

cat > example.lua <<EOF
require("lualammps")

lmp = lualammps.lammps_open_no_mpi(0,nil,nil)
ver = lualammps.lammps_version(lmp)

lualammps.lammps_command(lmp, "units real")
lualammps.lammps_command(lmp, "lattice fcc 2.5")
lualammps.lammps_command(lmp, "region box block -5 5 -5 5 -5 5")
lualammps.lammps_command(lmp, "create_box 1 box")
lualammps.lammps_command(lmp, "create_atoms 1 box")

print("LAMMPS version ", ver)
print("Number of created atoms: ", lualammps.lammps_get_natoms(lmp))
print("Current size of timestep: ", lualammps.double_p_value(lualammps.void_p_to_double_p(lualammps.lammps_extract_global(lmp,"dt"))))
lualammps.lammps_close(lmp)
EOF

lua example.lua
