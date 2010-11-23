#!/bin/sh 
#PBS -A TG-TRA100009
#PBS -N sds-all
#PBS -l size=96
#PBS -j oe
#PBS -l walltime=01:00:00

set -v
set -x

cd $PBS_O_WORKDIR
n=`expr $PBS_NNODES / 12`

s=1
m=`expr 2 \* $s`
t=`expr 12 / $m`
size=`expr $n \* $m`
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken -in in.sds-cpu -log log.sds-cpu-${s}m${t}t
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken-omp -in in.sds-omp -log log.sds-omp-${s}m${t}t

s=2
m=`expr 2 \* $s`
t=`expr 12 / $m`
size=`expr $n \* $m`
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken -in in.sds-cpu -log log.sds-cpu-${s}m${t}t
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken-omp -in in.sds-omp -log log.sds-omp-${s}m${t}t

s=3
m=`expr 2 \* $s`
t=`expr 12 / $m`
size=`expr $n \* $m`
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken -in in.sds-cpu -log log.sds-cpu-${s}m${t}t
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken-omp -in in.sds-omp -log log.sds-omp-${s}m${t}t

s=6
m=`expr 2 \* $s`
t=`expr 12 / $m`
size=`expr $n \* $m`
env OMP_NUM_THREADS=$t aprun -n $size -S $s -d $t ./lmp_kraken -in in.sds-cpu -log log.sds-cpu-${s}m${t}t

