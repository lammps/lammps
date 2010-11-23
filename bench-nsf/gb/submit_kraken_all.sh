#!/bin/sh 
#PBS -A <project>
#PBS -N gb-all
#PBS -l size=192
#PBS -j oe
#PBS -l walltime=06:00:00
#PBS

set -v
set -x

cd $PBS_O_WORKDIR
n=`expr $PBS_NNODES / 12`

s=1
m=`expr 2 \* $s`
t=`expr 12 / $m`
size=`expr $n \* $m`
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken -in in.gb-cpu -var len 40.0 -log log.gb-cpu-${t}m${s}t
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken-omp -in in.gb-omp -var len 40.0 -log log.gb-omp-${t}m${s}t

s=2
m=`expr 2 \* $s`
t=`expr 12 / $m`
size=`expr $n \* $m`
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken -in in.gb-cpu -var len 40.0 -log log.gb-cpu-${t}m${s}t
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken-omp -in in.gb-omp -var len 40.0 -log log.gb-omp-${t}m${s}t

s=3
m=`expr 2 \* $s`
t=`expr 12 / $m`
size=`expr $n \* $m`
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken -in in.gb-cpu -var len 40.0 -log log.gb-cpu-${t}m${s}t
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken-omp -in in.gb-omp -var len 40.0 -log log.gb-omp-${t}m${s}t

s=6
m=`expr 2 \* $s`
t=`expr 12 / $m`
size=`expr $n \* $m`
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken -in in.gb-cpu -var len 40.0 -log log.gb-cpu-${t}m${s}t
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken-omp -in in.gb-omp -var len 40.0 -log log.gb-omp-${t}m${s}t

