#!/bin/sh

cat > example.tcl <<EOF
load ./tcllammps.so

set lmp [lammps_open_no_mpi 0 NULL NULL ]
set ver [lammps_version \$lmp]

puts "LAMMPS version \$ver"
lammps_close \$lmp
EOF

tclsh example.tcl
