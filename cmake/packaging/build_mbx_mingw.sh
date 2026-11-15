#!/bin/bash
# Precompiled MBX for cross-compiling LAMMPS and the MBX plugin for Windows: compiling some MBX
# sources with optimization takes a long time, so cmake/Modules/Packages/MBX.cmake uses an
# archive of a compiled MBX where it lists one for the MBX version and the major version of
# the cross-compiler.  This helps to speed up testing of the installer builds quite a bit.
# This script makes that archive, compiled exactly as MBX.cmake compiles it (without MPI, with
# OpenMP), in the MinGW cross-compilation environment of the Windows installers (Fedora with the
# mingw64-* packages, e.g. tools/singularity/fedora41_musl_mingw.def).
#
#   build_mbx_mingw.sh <LAMMPS source folder> [<more CMake options>]
#
# It writes <MBX version>-mingw64-gcc<N>.tar.gz to the current folder and prints its SHA256.  Upload
# the archive to https://download.lammps.org/thirdparty/ and add the line that it prints to the list
# in MBX.cmake.
set -euo pipefail

if [ $# -lt 1 ] || [ ! -f "$1/cmake/Modules/Packages/MBX.cmake" ]; then
  sed -n '2,12p' "$0"
  exit 1
fi
src=$(cd "$1" && pwd)
shift
build=$(mktemp -d "${PWD}/mbx-mingw.XXXXXX")

mingw64-cmake -S "${src}/examples/PACKAGES/mbx/plugin" -B "${build}" -D LAMMPS_SOURCE_DIR="${src}/src" \
              -D CMAKE_BUILD_TYPE=Release -D BUILD_MPI=off -D BUILD_OMP=on -D MBX_MINGW_FROM_SOURCE=on "$@"
# the long compile is a few files, but the rest may as well run in parallel
MAKEFLAGS="-j$(nproc)" cmake --build "${build}" --target mbx_build

name=$(sed -n 's/^MBX_MINGW_ARCHIVE:INTERNAL=//p' "${build}/CMakeCache.txt")
stage="${build}/${name}"
mkdir "${stage}"
cp -a "${build}/mbx_build-prefix/include" "${build}/mbx_build-prefix/lib" "${stage}"
# the pkg-config file names the folders of this build; MBX.cmake does not use it
rm -rf "${stage}/lib/pkgconfig"
tar -C "${build}" --owner=0 --group=0 -czf "${name}.tar.gz" "${name}"
rm -rf "${build}"
echo "${name}.tar.gz"
echo "    set(MBX_MINGW_SHA256_${name} \"$(sha256sum "${name}.tar.gz" | cut -d' ' -f1)\")"
