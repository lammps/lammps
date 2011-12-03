#!/bin/sh
#
# this script modifies the OpenMP directives in LAMMPS
# in a way, so that they are also accepted by the PGI
# compilers. This modification incurs a performance
# penalty, though, so it is not part of the regular code.
#
# to apply the hack. go to src/, install all packages
# as needed and then run: sh ../tools/hack_openmp_for_pgi.sh
# this needs to be repeated after each "make package-update"
#
# Axel Kohlmeyer <akohlmey@gmail.com>

for f in *.h *.cpp
do \
   sed -e '/#pragma omp/s/^\(.*default\)(none)\(.*\)$/\1(shared)\2/' \
       -e '/#pragma omp/s/shared([a-z0-9,_]\+)//' \
       -i.bak $f
done
