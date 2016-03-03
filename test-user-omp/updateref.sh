#!/bin/sh

if [ $# == 0 ]
then
    dirs=`/bin/ls -d bench_*`
else
    dirs="$@"
fi

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

  egrep -v '(LAMMPS|OpenMP|MPI|serial|Memory|Loop|Perform|FFT)' $s \
    | sed -e 's/-0\.0000000000/0.0000000000 /g'   \
          -e '/^Section.*/,$d' -e '/^.*Pair time.*/,$d' \
    | gzip -9 > refoutput/${t}.gz
  popd
done
