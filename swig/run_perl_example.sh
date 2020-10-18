#!/bin/sh

cat > example.pl <<EOF
use pllammps;

\$lmp = pllammps::lammps_open_no_mpi(0,undef,undef);
\$ver = pllammps::lammps_version(\$lmp);

print("LAMMPS version ",\$ver,"\n");
pllammps::lammps_close(\$lmp)
EOF

PERL5LIB=$PWD:${PERL5LIB-PWD}

export PERL5LIB

perl example.pl
