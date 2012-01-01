#!/bin/sh

for s in *.restart-*-* ; do \
  t=`echo $s | sed -e 's/.restart-/.ref-/'`
  c=`echo $s | sed -e 's/.restart-/.coeff-/'`
  s=`echo $s | sed -e 's/.restart-/.data-/'`
  test -f $s || echo "no data file $s found for restart $r"
  test -f $s || continue
  grep -v LAMMPS $s | sed -e 's/-0\.0000000000/0.0000000000 /g' > ${t}
  test -f $c && grep -v LAMMPS $c | sed -e 's/-0\.0000000000/0.0000000000 /g' >> ${t}
  if [ -f refoutput/${t} ]
  then
     diff refoutput/${t} $t && rm $t
  else
     echo no reference output for $s
  fi
done
