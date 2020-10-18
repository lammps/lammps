#!/bin/sh

cat > example.pl <<EOF
use pllammps;

\$lmp = pllammps::lammps_open_no_mpi(0,undef,undef);
\$ver = pllammps::lammps_version(\$lmp);

pllammps::lammps_command(\$lmp, "units real");
pllammps::lammps_command(\$lmp, "lattice fcc 2.5");
pllammps::lammps_command(\$lmp, "region box block -5 5 -5 5 -5 5");
pllammps::lammps_command(\$lmp, "create_box 1 box");
pllammps::lammps_command(\$lmp, "create_atoms 1 box");

print("LAMMPS version ",\$ver,"\n");
print("Number of created atoms: ", pllammps::lammps_get_natoms(\$lmp), "\n");
pllammps::lammps_close(\$lmp)
EOF

PERL5LIB=$PWD:${PERL5LIB-PWD}

export PERL5LIB

perl example.pl
