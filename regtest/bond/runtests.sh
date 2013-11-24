#!/bin/sh

export OMP_NUM_THREADS=2

if [ $# -lt 1 ]
then
  cat <<EOF

  usage: $0 <LAMMPS MPI binary>

EOF
exit 0
fi

exe="$1"


runcheck () {
lmp=$1
sty=$2
coe=$3
arg=$4
log=$5

for n in 1 2 4 8
do \
    mpirun -x OMP_NUM_THREADS -np $n ${lmp} -i in.bond_angle ${sty} ${coe} "${arg}" -log ${log}-${n} -echo none -screen none
    grep -A2 ^Step ${log}-${n} | sed 's,-0.000000, 0.000000,g' > out
    diff -u ref/${log} out && rm out ${log}-${n} && echo OK ${log} ${n} || echo ERROR ${log} ${n}
    mpirun -x OMP_NUM_THREADS -np $n ${lmp} -i in.bond_angle ${sty} ${coe} "${arg}" -log ${log}-omp-${n} -echo none -screen none -suffix omp
    grep -A2 ^Step ${log}-omp-${n} | sed 's,-0.000000, 0.000000,g' > omp
    diff -u ref/${log} omp && rm omp ${log}-omp-${n} && echo OK ${log} omp ${n} || echo ERROR ${log} omp ${n}
done
}

rm -f *.data *.restart log.*

# hack to handle styles with a '/'
mkdir -p fene
mkdir -p harmonic/shift
mkdir -p none_cosine/shift
mkdir -p none_fourier

runcheck ${exe} '-v astyle harmonic'   '-v acoeff' '1 100.0 91.0'           log.none-harmonic
# first term of charmm angle is same as harmonic
runcheck ${exe} '-v astyle charmm'     '-v acoeff' '1 100.0 91.0  0.0 1.7'  log.none-charmm-a
runcheck ${exe} '-v astyle charmm'     '-v acoeff' '1 100.0 91.0 40.0 1.7'  log.none-charmm-b
# first term of quartic angle is same as harmonic
runcheck ${exe} '-v astyle quartic'    '-v acoeff' '1 91.0 100.0 0.0 0.0'   log.none-quartic-a
runcheck ${exe} '-v astyle quartic'    '-v acoeff' '1 91.0 100.0 -10.0 5.0' log.none-quartic-b
runcheck ${exe} '-v astyle cosine'     '-v acoeff' '1 10.0'                  log.none-cosine
runcheck ${exe} '-v astyle cosine/squared'   '-v acoeff' '1 200.0 101.0'     log.none-cosine-squared
runcheck ${exe} '-v astyle cosine/delta'     '-v acoeff' '1 200.0 101.0'     log.none-cosine-delta
runcheck ${exe} '-v astyle cosine/shift'     '-v acoeff' '1 200.0 101.0'     log.none-cosine-shift
runcheck ${exe} '-v astyle cosine/shift/exp' '-v acoeff' '1 100.0 101.0 2.0' log.none-cosine-shift-exp
runcheck ${exe} '-v astyle cosine/periodic'  '-v acoeff' '1 20.0 -1 4'       log.none-cosine-periodic
runcheck ${exe} '-v astyle fourier'          '-v acoeff' '1 50.0 2.0 1.0 .5' log.none-fourier-a
runcheck ${exe} '-v astyle fourier/simple'   '-v acoeff' '1 10.0 1.0 2.0'    log.none-fourier-simple-a
# these two are the same
runcheck ${exe} '-v astyle fourier'          '-v acoeff' '1 20.0 1.0 2.0 0.' log.none-fourier-b
runcheck ${exe} '-v astyle fourier/simple'   '-v acoeff' '1 20.0 2.0 1.0'    log.none-fourier-simple-b
# these two are the same
runcheck ${exe} '-v astyle fourier'          '-v acoeff' '1 20.0 1.0 0.0 2.' log.none-fourier-c
runcheck ${exe} '-v astyle fourier/simple'   '-v acoeff' '1 20.0 2.0 2.0'    log.none-fourier-simple-c
# same as cosine
runcheck ${exe} '-v astyle fourier'          '-v acoeff' '1 10.0 1.0 1.0 0.' log.none-fourier-d
runcheck ${exe} '-v astyle fourier/simple'   '-v acoeff' '1 10.0 1.0 1.0'    log.none-fourier-simple-d

# XXX: need special case input for angle styles class2 sdk table dipole

runcheck ${exe} '-v bstyle harmonic'    '-v bcoeff' '1 100.0 1.01'            log.harmonic-none
runcheck ${exe} '-v bstyle morse'       '-v bcoeff' '1 100.0 2.0 1.01'        log.morse-none
runcheck ${exe} '-v bstyle nonlinear'   '-v bcoeff' '1 100.0 1.01 2.0'        log.nonlinear-none
# first term of class2 is harmonic
runcheck ${exe} '-v bstyle class2'      '-v bcoeff' '1 1.01 100.0 0.0 0.0'    log.class2a-none
runcheck ${exe} '-v bstyle class2'      '-v bcoeff' '1 1.01 100.0 80.0 80.0'  log.class2b-none
runcheck ${exe} '-v bstyle fene'        '-v bcoeff' '1 50.0 1.01 0.5 0.9'     log.fene-none
runcheck ${exe} '-v bstyle fene/expand' '-v bcoeff' '1 50.0 1.01 0.5 0.9 0.0' log.fene-expand-a-none
runcheck ${exe} '-v bstyle fene/expand' '-v bcoeff' '1 50.0 1.01 0.5 0.9 0.1' log.fene-expand-b-none
runcheck ${exe} '-v bstyle harmonic/shift' '-v bcoeff' '1 100.0 1.01 2.5'     log.harmonic-shift-none
runcheck ${exe} '-v bstyle harmonic/shift/cut' '-v bcoeff' '1 100.0 1.01 2.5' log.harmonic-shift-cut-none

# XXX: need special case inputs for bond style quartic and table.

# cleanup
rm -r none_cosine
rm -r none_fourier
rm -r fene
rm -r harmonic

rm -f *.data *.restart
