#!/bin/sh

for s in log.* ; do \
  t=`echo $s | sed -e 's/log./ref./'`
  grep -v ^Loop $s | egrep -v '(OpenMP|MPI|serial)' | grep -v 'Memory usage' | sed -e 's/-0\.0000000000/0.0000000000 /g' > $t
  if [ -f refoutput/${t}.gz ]
  then
     zdiff -u refoutput/${t}.gz $t && rm $t
  else
     echo no reference output for $s
  fi
done
