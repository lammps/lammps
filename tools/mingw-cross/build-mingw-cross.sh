#!/bin/sh
# automated build script to build windows installers from the lammps sources

MINGW_BUILD_DIR=${HOME}/mingw-cross
NUMCPU=1
echo X-compiling LAMMPS for Windows in ${MINGW_BUILD_DIR} with ${NUMCPU} procs

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
for d in mingw32 mingw64 lammps-current
do \
  dir="${MINGW_BUILD_DIR}/${d}"
  rm -rf ${dir}
  mkdir -p "${dir}" || exit 2
done

pushd "${LAMMPS_PATH}"

git archive -v --format=tar --prefix=lammps-current/ HEAD \
    README LICENSE doc src lib python \
    examples/{README,ASPHERE,KAPPA,MC,VISCOSITY,dipole,peri,hugoniostat,colloid,crack,friction,msst,obstacle,body,sputter,pour,ELASTIC,neb,ellipse,flow,meam,min,indent,deposit,micelle,shear,srd,dreiding,eim,prd,rigid,COUPLE,peptide,melt,comb,tad,reax,balance,snap,USER/{atc,awpmd,misc,phonon,cg-cmm,sph,fep}} \
    bench potentials tools/*.cpp tools/*.f tools/mingw-cross tools/msi2lmp tools/createatoms tools/colvars \
    | tar -C ${MINGW_BUILD_DIR} -xvf -
popd

pushd ${MINGW_BUILD_DIR}/lammps-current

# need to build parallel wrapper libs first
pushd src/STUBS
make -f Makefile.mingw32-cross
make -f Makefile.mingw64-cross
popd

# start by building libraries
pushd lib

for d in atc awpmd colvars gpu linalg meam poems voronoi
do \
    pushd $d
    make -j${NUMCPU} mingw32-cross mingw64-cross mingw32-cross-mpi mingw64-cross-mpi || exit 3
    popd
done

popd

# now to the main source dir
pushd src

# configure installed packages
make yes-all no-kokkos no-kim no-user-cuda no-reax no-user-qmmm no-user-lb no-mpiio no-user-intel no-user-quip no-python no-user-h5md no-user-vtk
make -j${NUMCPU} mingw32-cross OMP=yes || exit 4
make -j${NUMCPU} mingw64-cross OMP=yes || exit 4
make yes-user-lb yes-mpiio
make -j${NUMCPU} mingw32-cross-mpi OMP=yes || exit 4
make -j${NUMCPU} mingw64-cross-mpi OMP=yes || exit 4
cp lmp_mingw32-cross ${MINGW_BUILD_DIR}/mingw32/lmp_serial.exe
cp lmp_mingw32-cross-mpi ${MINGW_BUILD_DIR}/mingw32/lmp_mpi.exe
cp lmp_mingw64-cross ${MINGW_BUILD_DIR}/mingw64/lmp_serial.exe
cp lmp_mingw64-cross-mpi ${MINGW_BUILD_DIR}/mingw64/lmp_mpi.exe

popd

# now build some utilities
pushd ${MINGW_BUILD_DIR}

TOOLDIR=lammps-current/tools
MINGW32FLAGS="-DLAMMPS_SMALLSMALL -O2 -march=i686  -mtune=generic -mfpmath=387 -mpc64 "
MINGW64FLAGS="-DLAMMPS_SMALLBIG   -O2 -march=core2 -mtune=core2   -mpc64 -msse2"

i686-w64-mingw32-g++   ${MINGW32FLAGS} -static -o mingw32/restart2data.exe ${TOOLDIR}/restart2data.cpp
x86_64-w64-mingw32-g++ ${MINGW64FLAGS} -static -o mingw64/restart2data.exe ${TOOLDIR}/restart2data.cpp

i686-w64-mingw32-g++   ${MINGW32FLAGS} -static -o mingw32/binary2txt.exe ${TOOLDIR}/binary2txt.cpp
x86_64-w64-mingw32-g++ ${MINGW64FLAGS} -static -o mingw64/binary2txt.exe ${TOOLDIR}/binary2txt.cpp

i686-w64-mingw32-gfortran   ${MINGW32FLAGS} -static -o mingw32/chain.exe ${TOOLDIR}/chain.f -lquadmath
x86_64-w64-mingw32-gfortran ${MINGW64FLAGS} -static -o mingw64/chain.exe ${TOOLDIR}/chain.f -lquadmath

i686-w64-mingw32-gfortran   ${MINGW32FLAGS} -static -o mingw32/createatoms.exe \
	${TOOLDIR}/createatoms/createAtoms.f -lquadmath
x86_64-w64-mingw32-gfortran ${MINGW64FLAGS} -static -o mingw64/createatoms.exe \
	${TOOLDIR}/createatoms/createAtoms.f -lquadmath

make -C ${TOOLDIR}/msi2lmp/src TARGET=${PWD}/mingw32/msi2lmp.exe \
	CC=i686-w64-mingw32-gcc CFLAGS="${MINGW32FLAGS}" clean
make -C ${TOOLDIR}/msi2lmp/src TARGET=${PWD}/mingw32/msi2lmp.exe \
	CC=i686-w64-mingw32-gcc CFLAGS="${MINGW32FLAGS}" LDFLAGS=-static
make -C ${TOOLDIR}/msi2lmp/src TARGET=${PWD}/mingw64/msi2lmp.exe \
	CC=x86_64-w64-mingw32-gcc CFLAGS="${MINGW64FLAGS}" clean
make -C ${TOOLDIR}/msi2lmp/src TARGET=${PWD}/mingw64/msi2lmp.exe \
	CC=x86_64-w64-mingw32-gcc CFLAGS="${MINGW64FLAGS}" LDFLAGS=-static

make -C ${TOOLDIR}/colvars EXT=.exe clean
make -C ${TOOLDIR}/colvars EXT=.exe CXX=i686-w64-mingw32-g++ CXXFLAGS="${MINGW32FLAGS}" LDFLAGS=-static
cp ${TOOLDIR}/colvars/*.exe ${PWD}/mingw32/
make -C ${TOOLDIR}/colvars EXT=.exe clean
make -C ${TOOLDIR}/colvars EXT=.exe CXX=x86_64-w64-mingw32-g++ CXXFLAGS="${MINGW64FLAGS}" LDFLAGS=-static
cp ${TOOLDIR}/colvars/*.exe ${PWD}/mingw64/

# assemble and customize installer scripts 
datestr=$(date +%Y%m%d)
cp ${TOOLDIR}/mingw-cross/lammps.nsis ${TOOLDIR}/mingw-cross/EnvVarUpdate.nsh .
cp ${TOOLDIR}/mingw-cross/Obj_mingw32/libOpenCL.dll mingw32
cp ${TOOLDIR}/mingw-cross/Obj_mingw64/libOpenCL.dll mingw64
cp ${TOOLDIR}/mingw-cross/Obj_mingw32/gzip.exe mingw32
cp ${TOOLDIR}/mingw-cross/Obj_mingw64/gzip.exe mingw64
cp ${TOOLDIR}/mingw-cross/Obj_mingw32/ffmpeg.exe mingw32
cp ${TOOLDIR}/mingw-cross/Obj_mingw64/ffmpeg.exe mingw64
cp lammps-current/lib/gpu/Obj_mingw32/ocl_get_devices mingw32/ocl_get_devices.exe
cp lammps-current/lib/gpu/Obj_mingw64/ocl_get_devices mingw64/ocl_get_devices.exe

# determine os vendor and release for installer tweaks.
if [ -f /etc/os-release ]
then \
  . /etc/os-release
  vendor="${NAME}"
  release="${VERSION_ID}"
else
  vendor=$(grep  release /etc/issue | cut -d \  -f 1)
  release=$(grep  release /etc/issue | cut -d \  -f 3)
fi
arch=$(uname -m)

# convert text files into CR/LF format.
unix2dos lammps-current/LICENSE lammps-current/README lammps-current/tools/msi2lmp/README
find lammps-current/{bench,examples,potentials} -type f -print | xargs unix2dos
find lammps-current/tools/msi2lmp/frc_files -type f -print | xargs unix2dos
# bulk rename README to README.txt
for f in $(find lammps-current/{tools,bench,examples,potentials} -name README -print)
do  mv -v $f $f.txt; done

# create up-to-date version of the manual
pushd lammps-current
make -C src pdf
popd

# build installers
LIBGCC=libgcc_s_sjlj-1.dll
makensis -DMINGW=/usr/i686-w64-mingw32/sys-root/mingw/bin/   \
    -DVERSION=${datestr} -DBIT=32 -DLIBGCC=${LIBGCC} lammps.nsis

# Fedora 19 ships with GCC-4.8.x which has different exception handling
# on 64-bit windows and thus uses a different name for libgcc
if [ "$vendor" = "Fedora" ] && [ $release -ge 19 ]
then
    LIBGCC=libgcc_s_seh-1.dll
fi
makensis -DMINGW=/usr/x86_64-w64-mingw32/sys-root/mingw/bin/ \
    -DVERSION=${datestr} -DBIT=64 -DLIBGCC=${LIBGCC} lammps.nsis

popd

exit 0

# clean up build and temporary directories (not yet)
for d in mingw32 mingw64 lammps-current
do \
  dir="${MINGW_BUILD_DIR}/${d}"
  rm -rf ${dir}
  mkdir -p "${dir}" || exit 2
done
