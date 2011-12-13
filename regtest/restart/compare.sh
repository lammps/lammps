#!/bin/sh

for s in *.data-*-* ; do \
  t=`echo $s | sed -e 's/.data-/.ref-/'`
  grep -v LAMMPS $s | sed -e 's/-0\.0000000000/0.0000000000 /g' > ${t}
  if [ -f refoutput/${t} ]
  then
     diff -u refoutput/${t} $t && rm $t
  else
     echo no reference output for $s
  fi
done
