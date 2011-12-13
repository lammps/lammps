#!/bin/sh

for s in *.data-*-* ; do \
  t=`echo $s | sed -e 's/.data-/.ref-/'`
  grep -v LAMMPS $s | sed -e 's/-0\.0000000000/0.0000000000 /g' > refoutput/${t}
done
