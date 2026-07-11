#!/bin/bash -e
NT=$1
PROG="./KokkosCore_PerformanceTest_Mempool"
COMMON_ARGS="--kokkos-threads=$NT --fill_stride=1 --fill_level=70 --chunk_span=5 --repeat_inner=100"

postproc() {
cat log | head -n 1 | rev | cut -d ' ' -f 1 | rev >> xvals
cat log | tail -n 1 | rev | cut -d ' ' -f 1 | rev >> yvals
}

for yset in 1 2 3
do
  rm -f xvals yvals
  for x in 1 2 4 8 16 32
  do
    echo "yset $yset x factor $x"
    $PROG $COMMON_ARGS --alloc_size=`expr $x \* 1000000` --super_size=`expr $x \* 100000` > log
    postproc
  done
  rm -f yvals$yset
  mv yvals yvals$yset
done

rm -f datapoints
paste -d',' xvals yvals1 yvals2 yvals3 > datapoints
