#!/bin/bash -e
NT=$1
PROG="./KokkosCore_PerformanceTest_Mempool"
COMMON_ARGS="--kokkos-threads=$NT --fill_stride=1 --alloc_size=10027008 --super_size=65536 --repeat_inner=100 --chunk_span=4 --repeat_outer=10"

postproc() {
cat log | grep "fill ops per second" | rev | cut -d ' ' -f 2 | rev >> yvals_fill
cat log | grep "cycle ops per second" | rev | cut -d ' ' -f 2 | rev >> yvals_cycle
}

rm -f xvals yvals_fill yvals_cycle
for x in 75 95
do
  echo "test fill level $x"
  echo $x >> xvals
  $PROG $COMMON_ARGS --fill_level=$x 2>&1 | tee log
  postproc
done

rm -f datapoints
paste xvals yvals_fill yvals_cycle > datapoints.txt
