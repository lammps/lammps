#!/bin/sh

cat > example.rb <<EOF
require 'rubylammps'

lmp = Rubylammps.lammps_open_no_mpi(0,nil,nil)
ver = Rubylammps.lammps_version(lmp)

Rubylammps.lammps_command(lmp, "units real")
Rubylammps.lammps_command(lmp, "lattice fcc 2.5")
Rubylammps.lammps_command(lmp, "region box block -5 5 -5 5 -5 5")
Rubylammps.lammps_command(lmp, "create_box 1 box")
Rubylammps.lammps_command(lmp, "create_atoms 1 box")

print("LAMMPS version ", ver, "\n")
print("Number of created atoms: ", Rubylammps.lammps_get_natoms(lmp), "\n")
Rubylammps.lammps_close(lmp)
EOF

export RUBYLIB=$PWD:${RUBYLIB-${PWD}}
ruby example.rb
