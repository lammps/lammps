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
    bench potentials tools/*.cpp tools/*.f tools/mingw-cross tools/msi2lmp \
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
    make mingw32-cross mingw64-cross mingw32-cross-mpi mingw64-cross-mpi || exit 3
    popd
done

popd

# now to the main source dir
pushd src

# configure installed packages
make yes-all no-kim no-user-cuda no-reax
make mingw32-cross mingw64-cross mingw32-cross-mpi mingw64-cross-mpi || exit 4
cp lmp_mingw32-cross ${MINGW_BUILD_DIR}/32bit/lmp_serial.exe
cp lmp_mingw64-cross ${MINGW_BUILD_DIR}/64bit/lmp_serial.exe
cp lmp_mingw32-cross-mpi ${MINGW_BUILD_DIR}/32bit-mpi/lmp_mpi.exe
cp lmp_mingw64-cross-mpi ${MINGW_BUILD_DIR}/64bit-mpi/lmp_mpi.exe

popd

# now build some utilities
pushd ${MINGW_BUILD_DIR}

TOOLDIR=lammps-current/tools
MINGW32FLAGS="-DLAMMPS_SMALLSMALL -O2 -march=i686  -mtune=generic -mfpmath=387 -mpc64"
MINGW64FLAGS="-DLAMMPS_SMALLBIG   -O2 -march=core2 -mtune=core2   -mpc64 -msse2"

i686-w64-mingw32-g++   ${MINGW32FLAGS} -o 32bit/restart2data.exe ${TOOLDIR}/restart2data.cpp
x86_64-w64-mingw32-g++ ${MINGW64FLAGS} -o 64bit/restart2data.exe ${TOOLDIR}/restart2data.cpp

i686-w64-mingw32-g++   ${MINGW32FLAGS} -o 32bit/binary2txt.exe ${TOOLDIR}/binary2txt.cpp
x86_64-w64-mingw32-g++ ${MINGW64FLAGS} -o 64bit/binary2txt.exe ${TOOLDIR}/binary2txt.cpp

i686-w64-mingw32-gfortran   ${MINGW32FLAGS} -o 32bit/chain.exe ${TOOLDIR}/chain.f
x86_64-w64-mingw32-gfortran ${MINGW64FLAGS} -o 64bit/chain.exe ${TOOLDIR}/chain.f

make -C ${TOOLDIR}/msi2lmp/src TARGET=${PWD}/32bit/msi2lmp.exe \
	CC=i686-w64-mingw32-g++ CFLAGS=${MINGW32FLAGS}
make -C ${TOOLDIR}/msi2lmp/src TARGET=${PWD}/64bit/msi2lmp.exe \
	CC=x86_64-w64-mingw32-g++ CFLAGS=${MINGW64FLAGS}

# assemble and customize installer scripts 
datestr=$(date +%Y%m%d)
cp ${TOOLDIR}/mingw-cross/lammps.nsis ${TOOLDIR}/mingw-cross/EnvVarUpdate.nsh .
cp ${TOOLDIR}/mingw-cross/Obj_mingw32/libOpenCL.dll 32bit
cp ${TOOLDIR}/mingw-cross/Obj_mingw64/libOpenCL.dll 64bit
cp lammps-current/lib/gpu/Obj_mingw32/ocl_get_devices 32bit/ocl_get_devices.exe
cp lammps-current/lib/gpu/Obj_mingw64/ocl_get_devices 64bit/ocl_get_devices.exe

# determine os vendor and release for installer tweaks.
vendor=$(grep  release /etc/issue | cut -d \  -f 1)
release=$(grep  release /etc/issue | cut -d \  -f 3)

# Fedora 19 ships with GCC-4.8.x which has different
# exception handling and thus uses a different name for libgcc
LIBGCC=libgcc_s_sjlj-1.dll
if [ "$vendor" = "Fedora" ] && [ $release -ge 19 ]
then
    LIBGCC=libgcc_s_seh-1.dll
fi

# convert text files into CR/LF format.
unix2dos lammps-current/LICENSE lammps-current/README
find lammps-current/{bench,examples,potentials} -type f -print | xargs unix2dos
find lammps-current/tools/msi2lmp/biosym_frc_files -type f -print | xargs unix2dos

# build installers
makensis -DMINGW=/usr/i686-w64-mingw32/sys-root/mingw/bin/   \
    -DVERSION=${datestr} -DBIT=32 -DLIBGCC=${LIBGCC} lammps.nsis
makensis -DMINGW=/usr/i686-w64-mingw32/sys-root/mingw/bin/   \
    -DVERSION=${datestr} -DBIT=32 -DLIBGCC=${LIBGCC} -DMPI=1 \
    lammps.nsis

makensis -DMINGW=/usr/x86_64-w64-mingw32/sys-root/mingw/bin/ \
    -DVERSION=${datestr} -DBIT=64 -DLIBGCC=${LIBGCC} lammps.nsis
makensis -DMINGW=/usr/x86_64-w64-mingw32/sys-root/mingw/bin/ \
    -DVERSION=${datestr} -DBIT=64 -DLIBGCC=${LIBGCC} -DMPI=1 \
    lammps.nsis

popd

exit 0

# clean up build and temporary directories (not yet)
for d in 32bit 32bit-mpi 64bit 64bit-mpi lammps-current
do \
  dir="${MINGW_BUILD_DIR}/${d}"
  rm -rf ${dir}
  mkdir -p "${dir}" || exit 2
done
