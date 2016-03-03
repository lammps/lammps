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
  for s in log.*-serial
  do \
    t=`echo $s | sed -e 's/log./ref./' -e 's,-serial,,'`

    if [ ! -f $s ]
    then
      popd
      continue
    fi

    egrep -v '(LAMMPS|OpenMP|MPI|serial|Memory|Loop|Perform|FFT)' $s \
      | sed -e 's/-0\.0000000000/0.0000000000 /g'   \
            -e '/^Section.*/,$d' -e '/^.*Pair time.*/,$d' \
      | gzip -9 > ${t}.gz
  done
  popd
done
