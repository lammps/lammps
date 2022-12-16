#!/bin/bash

LMP_BIN="$1"
NP="${2:-1}"
echo "MPI over $NP procs:"
for feat in conp etypes tf
do
  echo "Using base input file in.$feat:"
  echo "mat_inv, log excerpts:"
  logfile="log.algo_test.$NP.$feat"
  mpirun -np $NP $LMP_BIN -i in.$feat -l $logfile > /dev/null 2>&1
  grep -A2 'Per MPI rank' $logfile
  grep -B1 'Loop time' $logfile
  rm $logfile
  for cgtype in mat_cg cg
  do
    for tol in 1e-4 1e-5 1e-6
    do
      echo "$cgtype, tol = $tol, log excerpts:"
      logfile="log.algo_test.$NP.$feat.$cgtype.$tol"
      sed '/electrode/ s/$/ algo '"$cgtype"' '"$tol"'/' in.$feat > in.temp
      mpirun -np $NP $LMP_BIN -i in.temp -l $logfile > /dev/null 2>&1
      grep -A2 'Per MPI rank' $logfile
      grep -B1 'Loop time' $logfile
      rm $logfile
    done
  done
done
