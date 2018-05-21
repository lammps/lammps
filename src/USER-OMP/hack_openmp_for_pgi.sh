#!/bin/sh

for f in *.h *.cpp
do \
   sed -e '/#pragma omp/s/^\(.*default\)(none)\(.*\)$/\1(shared)\2/' \
       -e '/#pragma omp/s/shared([a-z0-9,_]\+)//' \
       -i.bak $f
done
