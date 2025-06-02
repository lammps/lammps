#!/bin/sh
rm -f *.dump *.restart? *.restart
../charmm2lammps.pl -border=2.0 -cmap=36 -ions -water all36_prot 1gb1
cp ../../../potentials/charmm36.cmap .
