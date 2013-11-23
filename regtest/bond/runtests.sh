#!/bin/sh

export OMP_NUM_THREADS=2

runcheck () {
exe=$1
sty=$2
coe=$3
arg=$4
log=$5

for n in 1 2 4 8
do \
    mpirun -np $n ${exe} -i in.bond_angle ${sty} ${coe} "${arg}" -log ${log}-${n} -echo none -screen none
    grep -A2 ^Step ${log}-${n} > out
    diff -u ref/${log} out && rm out ${log}-${n} && echo OK ${log} ${n} || echo ERROR ${log} ${n}
    mpirun -np $n ${exe} -i in.bond_angle ${sty} ${coe} "${arg}" -log ${log}-omp-${n} -echo none -screen none -suffix omp
    grep -A2 ^Step ${log}-omp-${n} > omp
    diff -u ref/${log} omp && rm omp ${log}-omp-${n} && echo OK ${log} omp ${n} || echo ERROR ${log} omp ${n}
done
}

rm -f *.data *.restart log.*

runcheck lmp_g++ '-v bstyle harmonic'    '-v bcoeff' '1 100.0 1.01'            log.harmonic-none
runcheck lmp_g++ '-v bstyle morse'       '-v bcoeff' '1 100.0 2.0 1.01'        log.morse-none
runcheck lmp_g++ '-v bstyle nonlinear'   '-v bcoeff' '1 100.0 1.01 2.0'        log.nonlinear-none

# first term of class2 is harmonic
runcheck lmp_g++ '-v bstyle class2'      '-v bcoeff' '1 1.01 100.0 0.0 0.0'    log.class2a-none
runcheck lmp_g++ '-v bstyle class2'      '-v bcoeff' '1 1.01 100.0 80.0 80.0'  log.class2b-none

runcheck lmp_g++ '-v bstyle fene'        '-v bcoeff' '1 50.0 1.01 0.5 0.9'     log.fene-none

# hack to handle styles with a '/'
mkdir -p fene
runcheck lmp_g++ '-v bstyle fene/expand' '-v bcoeff' '1 50.0 1.01 0.5 0.9 0.0' log.fene-expand-a-none
runcheck lmp_g++ '-v bstyle fene/expand' '-v bcoeff' '1 50.0 1.01 0.5 0.9 0.1' log.fene-expand-b-none
rm -r fene

mkdir -p harmonic/shift
#runcheck lmp_g++ '-v bstyle harmonic/shift' '-v bcoeff' '1 100.0 1.01 1.01'    log.harmonic-shift-a-none
#runcheck lmp_g++ '-v bstyle harmonic/shift' '-v bcoeff' '1 100.0 1.01 2.0'     log.harmonic-shift-b-none
#runcheck lmp_g++ '-v bstyle harmonic/shift/cut' '-v bcoeff' '1 100.0 1.01 1.01' log.harmonic-shift-cut-a-none
#runcheck lmp_g++ '-v bstyle harmonic/shift/cut' '-v bcoeff' '1 100.0 1.01 2.0'  log.harmonic-shift-cut-b-none

# XXX: need special case inputs for bond style quartic and table.



rm -f *.data *.restart
