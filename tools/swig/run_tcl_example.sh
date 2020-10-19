#!/bin/sh

cat > example.tcl <<EOF
load ./tcllammps.so

set lmp [lammps_open_no_mpi 0 NULL NULL ]
set ver [lammps_version \$lmp]

lammps_command \$lmp "units real"
lammps_command \$lmp "lattice fcc 2.5"
lammps_command \$lmp "region box block -5 5 -5 5 -5 5"
lammps_command \$lmp "create_box 1 box"
lammps_command \$lmp "create_atoms 1 box"

set dt [double_p_value [void_p_to_double_p [lammps_extract_global \$lmp dt]]]
puts "LAMMPS version \$ver"
puts [format "Number of created atoms: %g" [lammps_get_natoms \$lmp]]
puts "Current size of timestep: \$dt"
lammps_close \$lmp
EOF

tclsh example.tcl
