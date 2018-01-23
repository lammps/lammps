#!/bin/bash
if test "`basename $PWD`" = "cmaketest"; then 
  outfile=$1
else
  outfile=config/tmpstore/$1
fi

grep_arch=`grep KOKKOS_ARCH    $outfile | grep $2 2>&1`
grep_devs=`grep KOKKOS_DEVICES $outfile | grep $3 2>&1`
if test -n "$grep_arch"; then 
  if test -n "$grep_devs"; then 
    echo Passed
  else
    echo Failed
  fi
else
  echo Failed
fi
