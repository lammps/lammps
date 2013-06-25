#!/bin/sh
# automated build script to build windows installers from the lammps sources

MINGW_BUILD_DIR=${HOME}/mingw-cross
echo X-compiling LAMMPS for Windows in ${MINGW_BUILD_DIR}

for d in "${PWD}" "${PWD%/src}" "${PWD%/tools/mingw-cross}" "$1"
do \
    if test -d "${d}/.git"
    then
        LAMMPS_PATH="${d}"
        break
    fi
done

if test -z "${LAMMPS_PATH}"
then
    echo "'${PWD}' is not a suitable working directory"
    exit 1
fi

for d in 32bit 32bit-mpi 64bit 64bit-mpi
do \
  dir="${MINGW_BUILD_DIR}/${d}"
  test -d "${dir}" || mkdir -p "${dir}" || exit 2
done

pushd "${LAMMPS_PATH}"

datestr=$(date +%Y%m%d)
sed -e "s/@VERSION@/${datestr}/" tools/mingw/win32-serial.nsis \
    > ${MINGW_BUILD_DIR}/32bit/lammps.nsis
sed -e "s/@VERSION@/${datestr}/" tools/mingw/win64-serial.nsis \
    > ${MINGW_BUILD_DIR}/64bit/lammps.nsis

git archive -v --format=tar --prefix=lammps-current/ HEAD \
    README LICENSE doc/Manual.pdf doc/PDF src lib python  \
    examples/{README,dipole,peri,hugoniostat,colloid,crack,friction,msst,obstacle,body,sputter,pour,ELASTIC,neb,ellipse,flow,meam,min,indent,micelle,shear,srd,dreiding,eim,prd,rigid,COUPLE,peptide,melt,comb,tad,reax,USER/{awpmd,misc,phonon,cg-cmm}} \
    bench potentials tools/*.cpp tools/*.f \
    | tar -C ${MINGW_BUILD_DIR} -xvf -
popd

pushd ${MINGW_BUILD_DIR}/lammps-current

# start by building libraries
pushd lib

for d in atc awpmd colvars linalg meam poems voronoi
do \
    pushd $d
    make mingw32-cross mingw64-cross || exit 1
    popd
done

popd

# now to the main source dir
pushd src
pushd STUBS
make -f Makefile.mingw32-cross
make -f Makefile.mingw64-cross
popd

make yes-all no-kim no-gpu no-user-cuda no-reax no-user-atc no-user-awpmd
make mingw32-cross mingw64-cross || exit 2
cp lmp_mingw32-cross ${MINGW_BUILD_DIR}/32bit/lmp_serial.exe
cp lmp_mingw64-cross ${MINGW_BUILD_DIR}/64bit/lmp_serial.exe

popd

# now build some utilities
pushd ${MINGW_BUILD_DIR}/32bit

i686-w64-mingw32-g++ -o restart2data.exe -DLAMMPS_SMALLSMALL -O2 -march=i686 \
    -mtune=generic -mfpmath=387 -mpc64 ${LAMMPS_PATH}/tools/restart2data.cpp
cp restart2data.exe ../32bit-mpi/

i686-w64-mingw32-g++ -o binary2txt.exe -DLAMMPS_SMALLSMALL -O2 -march=i686 \
    -mtune=generic -mfpmath=387 -mpc64 ${LAMMPS_PATH}/tools/binary2txt.cpp
cp binary2txt.exe ../32bit-mpi/

i686-w64-mingw32-gfortran -o chain.exe -O2 -march=i686 -mtune=generic \
    -mfpmath=387 -mpc64 ${LAMMPS_PATH}/tools/chain.f
cp chain.exe ../32bit-mpi/chain.exe

cd ../64bit

x86_64-w64-mingw32-g++ -o restart2data.exe -DLAMMPS_SMALLBIG -O2 -march=core2 \
    -mtune=core2 -mpc64 -msse2 ${LAMMPS_PATH}/tools/restart2data.cpp
cp restart2data.exe ../64bit-mpi/

x86_64-w64-mingw32-g++ -o binary2txt.exe -DLAMMPS_SMALLBIG -O2 -march=core2 \
    -mtune=core2 -mpc64 -msse2 ${LAMMPS_PATH}/tools/binary2txt.cpp
cp binary2txt.exe ../64bit-mpi/

x86_64-w64-mingw32-gfortran -o ${MINGW_BUILD_DIR}/64bit/chain.exe -O2 \
    -march=core2 -mtune=core2 -mpc64 -msse2 ${LAMMPS_PATH}/tools/chain.f
cp chain.exe ../64bit-mpi/chain.exe

cd ../32bit
makensis lammps.nsis
cd ../64bit
makensis lammps.nsis

popd

