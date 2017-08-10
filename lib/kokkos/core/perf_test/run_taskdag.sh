#!/bin/bash -e
NT=$1
PROG="./KokkosCore_PerformanceTest_TaskDAG"
COMMON_ARGS="--kokkos-threads=$NT --alloc_size=10027008 --super_size=65536 --repeat_outer=10"

postproc() {
cat log | grep "tasks per second" | rev | cut -d ' ' -f 2 | rev >> yvals
}

rm -f xvals yvals
for x in 21 23
do
  echo "test input $x"
  echo $x >> xvals
  $PROG $COMMON_ARGS --input=$x 2>&1 | tee log
  postproc
done

rm -f datapoints.txt
paste xvals yvals > datapoints.txt

