#!/bin/sh

dirs=`/bin/ls -d bench_*`

for d in $dirs
do \
  pushd $d
  for f in [lr][oe][gf].*xOpenMP_*xMPI ; do \
    rm -f $f
  done
  popd
done
