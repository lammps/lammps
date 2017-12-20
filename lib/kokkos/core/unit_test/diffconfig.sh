#!/bin/bash
# CMake and Make tests run in separate directories
#   The mapping of ARCH to #define is very complicated
#       so diff is used instead of grepping
if test "`basename $PWD`" = "cmaketest"; then 
  outfile=$1
  resfile=../results/$1
else
  outfile=config/tmpstore/$1
  resfile=config/results/$1
fi

diff=`diff $outfile $resfile 2>&1 | grep -e define -e "such file"`
if test -z "$diff"; then 
  echo Passed
else
  echo Failed: $diff
fi
