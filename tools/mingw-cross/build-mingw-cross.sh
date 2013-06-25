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

# clean up leftovers from an old build and rebuild directories
for d in 32bit 32bit-mpi 64bit 64bit-mpi lammps-current
do \
  dir="${MINGW_BUILD_DIR}/${d}"
  rm -rf ${dir}
  mkdir -p "${dir}" || exit 2
done

pushd "${LAMMPS_PATH}"

git archive -v --format=tar --prefix=lammps-current/ HEAD \
    README LICENSE doc/Manual.pdf doc/PDF src lib python  \
    examples/{README,dipole,peri,hugoniostat,colloid,crack,friction,msst,obstacle,body,sputter,pour,ELASTIC,neb,ellipse,flow,meam,min,indent,micelle,shear,srd,dreiding,eim,prd,rigid,COUPLE,peptide,melt,comb,tad,reax,USER/{awpmd,misc,phonon,cg-cmm}} \
    bench potentials tools/*.cpp tools/*.f tools/mingw-cross \
    | tar -C ${MINGW_BUILD_DIR} -xvf -
popd

pushd ${MINGW_BUILD_DIR}/lammps-current

# start by building libraries
pushd lib

for d in atc awpmd colvars linalg meam poems voronoi
do \
    pushd $d
    make mingw32-cross mingw64-cross mingw32-cross-mpi mingw64-cross-mpi || exit 3
    popd
done

popd

# now to the main source dir and build MPI stub libraries
pushd src
pushd STUBS
make -f Makefile.mingw32-cross
make -f Makefile.mingw64-cross
popd

# configure installed packages
make yes-all no-kim no-gpu no-user-cuda no-reax
make mingw32-cross mingw64-cross mingw32-cross-mpi mingw64-cross-mpi || exit 4
cp lmp_mingw32-cross ${MINGW_BUILD_DIR}/32bit/lmp_serial.exe
cp lmp_mingw64-cross ${MINGW_BUILD_DIR}/64bit/lmp_serial.exe
cp lmp_mingw32-cross-mpi ${MINGW_BUILD_DIR}/32bit-mpi/lmp_mpi.exe
cp lmp_mingw64-cross-mpi ${MINGW_BUILD_DIR}/64bit-mpi/lmp_mpi.exe

popd

# now build some utilities
pushd ${MINGW_BUILD_DIR}

i686-w64-mingw32-g++ -o 32bit/restart2data.exe -DLAMMPS_SMALLSMALL -O2 -march=i686 \
    -mtune=generic -mfpmath=387 -mpc64 lammps-current/tools/restart2data.cpp
x86_64-w64-mingw32-g++ -o 64bit/restart2data.exe -DLAMMPS_SMALLBIG -O2 -march=core2 \
    -mtune=core2 -mpc64 -msse2 lammps-current/tools/restart2data.cpp

i686-w64-mingw32-g++ -o 32bit/binary2txt.exe -DLAMMPS_SMALLSMALL -O2 -march=i686 \
    -mtune=generic -mfpmath=387 -mpc64 lammps-current/tools/binary2txt.cpp
x86_64-w64-mingw32-g++ -o 64bit/binary2txt.exe -DLAMMPS_SMALLBIG -O2 -march=core2 \
    -mtune=core2 -mpc64 -msse2 lammps-current/tools/binary2txt.cpp

i686-w64-mingw32-gfortran -o 32bit/chain.exe -O2 -march=i686 -mtune=generic \
    -mfpmath=387 -mpc64 lammps-current/tools/chain.f
x86_64-w64-mingw32-gfortran -o 64bit/chain.exe -O2 -march=core2 -mtune=core2 \
    -mpc64 -msse2 lammps-current/tools/chain.f

# assemble and customize installer scripts 
datestr=$(date +%Y%m%d)
cp lammps-current/tools/mingw-cross/win??-*.nsis .
cp lammps-current/tools/mingw-cross/EnvVarUpdate.nsh .
sed -i -e "s/@VERSION@/${datestr}/g" win??-*.nsis

# determine os vendor and release for installer tweaks.
vendor=$(grep  release /etc/issue | cut -d \  -f 1)
release=$(grep  release /etc/issue | cut -d \  -f 1)

# Fedora 19 ships with GCC-4.8.x which has different exception handling in libgcc
if [ "$vendor" = "Fedora" ] && [ $release -ge 19 ]
then
    sed -i -e "s/libgcc_s_sjlj-1.dll/libgcc_s_seh-1.dll/g" win64-*.nsis
fi

# build installers
makensis win32-serial.nsis
makensis win32-mpi.nsis
makensis win64-serial.nsis
makensis win64-mpi.nsis

popd

exit 0

# clean up build and temporary directories (not yet)
for d in 32bit 32bit-mpi 64bit 64bit-mpi lammps-current
do \
  dir="${MINGW_BUILD_DIR}/${d}"
  rm -rf ${dir}
  mkdir -p "${dir}" || exit 2
done
