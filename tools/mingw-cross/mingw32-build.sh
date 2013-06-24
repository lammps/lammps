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
# XXX
#sed -e "/^Version/s/\(Version:[ 	]\+\)[0-9].*$/\1${datestr}/" tools/mingw/win32-serial.nsis > ${MINGW_BUILD_DIR}/32bit/lammps.nsis

git archive -v --format=tar --prefix=lammps-current/ HEAD \
    README LICENSE doc/Manual.pdf doc/PDF src lib python  \
    examples/{README,dipole,peri,hugoniostat,colloid,crack,friction,msst,obstacle,body,sputter,pour,ELASTIC,neb,ellipse,flow,meam,min,indent,micelle,shear,srd,dreiding,eim,prd,rigid,COUPLE,peptide,melt,comb,tad,reax,USER/{awpmd,misc,phonon,cg-cmm}} \
    bench potentials tools/*.cpp tools/*.f \
    | tar -C ${MINGW_BUILD_DIR} -xvf -
popd

pushd ${MINGW_BUILD_DIR}/lammps-current

# start by building libraries
pushd lib

#cd lib/atc
#make -f Makefile.g++ CC=g++ CCFLAGS="-fPIC -I../../src -I../../src/STUBS  ${RPM_OPT_FLAGS} %{bigintsize}" EXTRAMAKE=Makefile.lammps.linalg
#cd ../awpmd
#make -f Makefile.openmpi CC=g++ CCFLAGS="-fPIC -Isystems/interact/TCP/ -Isystems/interact -Iivutils/include ${RPM_OPT_FLAGS} %{bigintsize}" EXTRAMAKE=Makefile.lammps.linalg

for d in colvars linalg meam poems voronoi
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
i686-w64-mingw32-g++ -o ${MINGW_BUILD_DIR}/32bit/restart2data.exe -DLAMMPS_SMALLSMALL \
	-O2 -march=i686 -mtune=generic -mfpmath=387 -mpc64 tools/restart2data.cpp
cp ${MINGW_BUILD_DIR}/32bit/restart2data.exe ${MINGW_BUILD_DIR}/32bit-mpi/restart2data.exe
x86_64-w64-mingw32-g++ -o ${MINGW_BUILD_DIR}/64bit/restart2data.exe -DLAMMPS_SMALLBIG \
	-O2 -march=core2 -mtune=core2 -mpc64 -msse2 tools/restart2data.cpp
cp ${MINGW_BUILD_DIR}/64bit/restart2data.exe ${MINGW_BUILD_DIR}/64bit-mpi/restart2data.exe

i686-w64-mingw32-g++ -o ${MINGW_BUILD_DIR}/32bit/binary2txt.exe -DLAMMPS_SMALLSMALL \
	-O2 -march=i686 -mtune=generic -mfpmath=387 -mpc64 tools/binary2txt.cpp
cp ${MINGW_BUILD_DIR}/32bit/binary2txt.exe ${MINGW_BUILD_DIR}/32bit-mpi/binary2txt.exe
x86_64-w64-mingw32-g++ -o ${MINGW_BUILD_DIR}/64bit/binary2txt.exe -DLAMMPS_SMALLBIG \
	-O2 -march=core2 -mtune=core2 -mpc64 -msse2 tools/binary2txt.cpp
cp ${MINGW_BUILD_DIR}/64bit/binary2txt.exe ${MINGW_BUILD_DIR}/64bit-mpi/binary2txt.exe

i686-w64-mingw32-gfortran -o ${MINGW_BUILD_DIR}/32bit/chain.exe -DLAMMPS_SMALLSMALL \
	-O2 -march=i686 -mtune=generic -mfpmath=387 -mpc64 tools/chain.f
cp ${MINGW_BUILD_DIR}/32bit/chain.exe ${MINGW_BUILD_DIR}/32bit-mpi/chain.exe
x86_64-w64-mingw32-gfortran -o ${MINGW_BUILD_DIR}/64bit/chain.exe -DLAMMPS_SMALLBIG \
	-O2 -march=core2 -mtune=core2 -mpc64 -msse2 tools/chain.f
cp ${MINGW_BUILD_DIR}/64bit/chain.exe ${MINGW_BUILD_DIR}/64bit-mpi/chain.exe


