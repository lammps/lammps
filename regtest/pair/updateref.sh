#!/bin/sh

for s in log.{simple,charged,overlay}* ; do \
  t=`echo $s | sed -e 's/log./ref./'`
  grep -v ^Loop $s | egrep -v '(LAMMPS|OpenMP|MPI|serial)' | grep -v 'Memory usage' | sed -e 's/-0\.0000000000/0.0000000000 /g' | gzip -9 > refoutput/${t}.gz
done
