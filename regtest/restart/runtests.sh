#!/bin/sh

if [ $# -lt 2 ]
then
   cat <<EOF

  usage: $0 <serial LAMMPS executable> <restart2data executable>

EOF
exit 0
fi

exe="$1"

export OMP_NUM_THREADS=1

for s in nosfx opt omp
do \
  for t in 01 02
  do \
    ${exe} -log log.restart2data-${s}-${t} -echo log -screen none \
      -in in.restart2data -var tst ${t} -var sfx ${s} -var r2d "$2"
  done
done
