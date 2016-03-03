#!/bin/sh 

dirs=`/bin/ls -d bench_*`

for d in $dirs
do \
  pushd $d

  s=log.0xOpenMP_1xMPI
  t=`echo $s | sed -e 's/log./ref./'`
  if [ ! -f $s ]
  then
    popd
    continue
  fi

  egrep -v '(LAMMPS|OpenMP|MPI|serial|Memory|Loop|Perform|FFT)' $s   \
    | sed -e 's/-0\.0000000000/0.0000000000 /g' -e '/^Section.*/,$d' \
    > $t

  if [ -f refoutput/${t}.gz ]
  then
    zdiff -u refoutput/${t}.gz $t && rm $t
  else
    echo no reference output for $s
  fi

  for l in log.?xOpenMP*; do \
    [ $l == $s ] && continue
    
    k=`echo $l | sed -e 's/log./ref./'`

    egrep -v '(LAMMPS|OpenMP|MPI|serial|Memory|Loop|Perform|FFT|/omp)' $l   \
      | sed -e 's/-0\.0000000000/0.0000000000 /g' -e '/^Section.*/,$d' \
      > $k

    if [ -f refoutput/${t}.gz ]
    then
      zdiff -u refoutput/${t}.gz $k && rm $k
    else
      echo no reference output for $l
    fi

  done
  popd
done
