#!/bin/sh

NPROCS=1
if [ $# -gt 0 ]; then
  NPROCS=$1
fi

bash ./clean.sh

python ./double-re-short.py $NPROCS $HOME/compile/lammps-icms/src/lmp_omp in.gREM > total_output.$NPROCS

exit 0
