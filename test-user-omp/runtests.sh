#!/bin/sh

if [ $# -lt 1 ]
then
  cat <<EOF

  usage: $0 <LAMMPS binary>

EOF
exit 0
fi

exe="$@"

dirs=`/bin/ls -d bench_*`

for d in $dirs
do \
  pushd $d
  in=`/bin/ls in.* | head -1`
  for s in 1 2 4 8
  do \
    mpirun -np ${s} -x OMP_NUM_THREADS=1 $exe -log log.0xOpenMP_${s}xMPI -in $in -echo none
  done
  t=1
  for s in 1 2 4 8
  do \
    mpirun -np ${s} -x OMP_NUM_THREADS=${t} $exe -log log.${t}xOpenMP_${s}xMPI -in $in -echo none -sf omp
  done
  t=2
  for s in 1 2 4 
  do \
    mpirun -np ${s} -x OMP_NUM_THREADS=${t} $exe -log log.${t}xOpenMP_${s}xMPI -in $in -echo none -sf omp
  done
  t=4
  for s in 1 2 
  do \
    mpirun -np ${s} -x OMP_NUM_THREADS=${t} $exe -log log.${t}xOpenMP_${s}xMPI -in $in -echo none -sf omp
  done
  t=8
  mpirun -np 1 -x OMP_NUM_THREADS=${t} $exe -log log.${t}xOpenMP_1xMPI -in $in -echo none -sf omp
  popd
done
