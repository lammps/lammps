#!/bin/sh

for s in *.data-*-* ; do \
  t=`echo $s | sed -e 's/.data-/.ref-/'`
  c=`echo $s | sed -e 's/.data-/.coeff-/'`
  grep -v LAMMPS $s | sed -e 's/-0\.0000000000/0.0000000000 /g' > ${t}
  test -f $c && grep -v LAMMPS $c | sed -e 's/-0\.0000000000/0.0000000000 /g' >> ${t}
  if [ -f refoutput/${t} ]
  then
     diff refoutput/${t} $t && rm $t
  else
     echo no reference output for $s
  fi
done
