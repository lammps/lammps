#!/bin/sh

for r in *.restart-*-* ; do \
  t=`echo $r | sed -e 's/.restart-/.ref-/'`
  c=`echo $r | sed -e 's/.restart-/.coeff-/'`
  s=`echo $r | sed -e 's/.restart-/.data-/'`
  test -f $s || echo "no data file $s found for restart $r" && continue
  grep -v LAMMPS $s | sed -e 's/-0\.0000000000/0.0000000000 /g' > refoutput/${t}
  test -f $c && grep -v LAMMPS $c | sed -e 's/-0\.0000000000/0.0000000000 /g' >> refoutput/${t}
done
