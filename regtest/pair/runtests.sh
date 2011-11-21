#!/bin/sh

if [ $# -lt 1 ]
then
   cat <<EOF

  usage: $0 <command(s) to run LAMMPS>

EOF
exit 0
fi

exe="$@"

export OMP_NUM_THREADS=1

for tst in 01 02 ; do \
  for sfx in nosfx opt omp-1-no-neigh omp-1-neigh omp-2-no-neigh \
    omp-2-neigh omp-4-no-neigh omp-4-neigh ; do \
    ${exe} -log log.simple-${tst}-${sfx} -echo none -screen none \
           -in in.simple -var sfx ${sfx} -var tst ${tst}
  done
done

for tst in 01 02 ; do \
  for sfx in nosfx opt omp-1-no-neigh omp-1-neigh omp-2-no-neigh \
    omp-2-neigh omp-4-no-neigh omp-4-neigh ; do \
    ${exe} -log log.charged-${tst}-${sfx} -echo none -screen none \
           -in in.charged -var sfx ${sfx} -var tst ${tst}
  done
done

for tst in 01 02 ; do \
  for sfx in nosfx opt omp-1-no-neigh omp-1-neigh omp-2-no-neigh \
    omp-2-neigh omp-4-no-neigh omp-4-neigh ; do \
    ${exe} -log log.overlay-${tst}-${sfx} -echo none -screen none \
           -in in.overlay -var sfx ${sfx} -var tst ${tst}
  done
done

