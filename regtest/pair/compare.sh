#!/bin/sh

for s in log.* ; do \
  t=`echo $s | sed -e 's/log./ref./'`
  grep -v ^Loop $s | egrep -v '(LAMMPS|OpenMP|MPI|serial)' | grep -v 'Memory usage' | sed -e 's/-0\.0000000000/0.0000000000 /g' > $t
  if [ -f refoutput/${t}.gz ]
  then
     zdiff -u refoutput/${t}.gz $t && rm $t
  else
     echo no reference output for $s
  fi
done

for s in log.overlay* ; do \
  r=`echo $s | sed -e 's/overlay/charged/'`
  t=`echo $r | sed -e 's/log.charged/ref.ovc-charged/'`
  u=`echo $s | sed -e 's/log.overlay/ref.ovc-overlay/'`

  grep -v ^Loop $r | egrep -v '(LAMMPS|OpenMP|MPI|serial)' | grep -v 'Memory usage' | sed -e 's/-0\.0000000000/0.0000000000 /g' > $t
  grep -v ^Loop $s | egrep -v '(LAMMPS|OpenMP|MPI|serial)' | grep -v 'Memory usage' | sed -e 's/-0\.0000000000/0.0000000000 /g' > $u

  diff -u $t $u && rm $t $u
done

