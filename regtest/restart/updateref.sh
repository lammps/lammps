#!/bin/sh

for s in *.data-*-* ; do \
  t=`echo $s | sed -e 's/.data-/.ref-/'`
  c=`echo $s | sed -e 's/.data-/.coeff-/'`
  grep -v LAMMPS $s | sed -e 's/-0\.0000000000/0.0000000000 /g' > refoutput/${t}
  test -f $c && grep -v LAMMPS $c | sed -e 's/-0\.0000000000/0.0000000000 /g' >> refoutput/${t}
done
