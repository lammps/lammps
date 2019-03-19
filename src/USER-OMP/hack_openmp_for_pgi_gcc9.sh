#!/bin/sh

# convert default(none) directives for OpenMP pragmas to default(shared) and remove shared() directive
# this allows compiling OpenMP pragmas in LAMMPS with compilers that don't support default(none) properly
# or require backward incompatible OpenMP 4 and OpenMP 5 semantics

for f in *.h *.cpp
do \
   sed -e '/#pragma omp/s/^\(.*default\)(none)\(.*\)$/\1(shared)\2/' \
       -e '/#pragma omp/s/shared([a-z0-9,_]\+)//' \
       -i.bak $f
done
