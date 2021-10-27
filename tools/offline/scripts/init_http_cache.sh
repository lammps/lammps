#!/bin/bash

if [ -z "${HTTP_CACHE_DIR}" ]
then
    echo "Must set HTTP_CACHE_DIR environment variable"
    exit 1
fi

mkdir -p "$HTTP_CACHE_DIR"
mkdir -p "$HTTP_CACHE_DIR/potentials"
mkdir -p "$HTTP_CACHE_DIR/thirdparty"
cd $HTTP_CACHE_DIR

LAMMPS_DOWNLOADS_URL="https://download.lammps.org"
LAMMPS_POTENTIALS_URL="${LAMMPS_DOWNLOADS_URL}/potentials"
LAMMPS_THIRDPARTY_URL="${LAMMPS_DOWNLOADS_URL}/thirdparty"

###############################################################################
# potentials
POTENTIALS=(
    C_10_10.mesocnt.028de73ec828b7830d762702eda571c1
    TABTP_10_10.mesont.744a739da49ad5e78492c1fc9fd9f8c1
    C_10_10.mesocnt
    TABTP_10_10.mesont
)

echo "Dowloading potentials..."
for p in ${POTENTIALS[@]}
do
    if [ ! -f "$HTTP_CACHE_DIR/potentials/$p" ]
    then
        wget -O "$HTTP_CACHE_DIR/potentials/$p" "$LAMMPS_POTENTIALS_URL/$p"
    fi
done

###############################################################################
# thirdparty code
echo "Dowloading thirdparty tarballs..."

MPICH2_WIN64_DEVEL_URL="${LAMMPS_THIRDPARTY_URL}/mpich2-win64-devel.tar.gz"
MPICH2_WIN32_DEVEL_URL="${LAMMPS_THIRDPARTY_URL}/mpich2-win32-devel.tar.gz"
VORO_URL="${LAMMPS_THIRDPARTY_URL}/voro++-0.4.6.tar.gz"
OPENCL_LOADER_URL="${LAMMPS_THIRDPARTY_URL}/opencl-loader-2021.09.18.tar.gz"
SCAFACOS_FIX_URL="${LAMMPS_THIRDPARTY_URL}/scafacos-1.0.1-fix.diff"
GTEST_URL="https://github.com/google/googletest/archive/release-1.11.0.tar.gz"
YAML_URL="https://pyyaml.org/download/libyaml/yaml-0.2.5.tar.gz"
MATHJAX_URL="https://github.com/mathjax/MathJax/archive/3.1.3.tar.gz"
EIGEN3_URL="${LAMMPS_THIRDPARTY_URL}/eigen-3.4.0.tar.gz"
CUB_URL="https://github.com/NVlabs/cub/archive/1.12.0.tar.gz"
KOKKOS_URL="https://github.com/kokkos/kokkos/archive/3.4.01.tar.gz"
KIM_URL="https://s3.openkim.org/kim-api/kim-api-2.2.1.txz"
MSCG_URL="https://github.com/uchicago-voth/MSCG-release/archive/1.7.3.1.tar.gz"
PLUMED_URL="https://github.com/plumed/plumed2/releases/download/v2.7.2/plumed-src-2.7.2.tgz"
PACELIB_URL="https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2021.10.25.tar.gz"
LATTE_URL="https://github.com/lanl/LATTE/archive/v1.2.2.tar.gz"
SCAFACOS_URL="https://github.com/scafacos/scafacos/releases/download/v1.0.1/scafacos-1.0.1.tar.gz"
MDI_URL="https://github.com/MolSSI-MDI/MDI_Library/archive/v1.2.9.tar.gz"
N2P2_URL="https://github.com/CompPhysVienna/n2p2/archive/v2.1.4.tar.gz"

GTEST_FILENAME="gtest-1.11.0.tar.gz"
MATHJAX_FILENAME="mathjax-3.1.3.tar.gz"
CUB_FILENAME="cub-1.12.0.tar.gz"
KOKKOS_FILENAME="kokkos-3.4.01.tar.gz"
MSCG_FILENAME="mscg-1.7.3.1.tar.gz"
LATTE_FILENAME="latte-1.2.2.tar.gz"
N2P2_FILENAME="n2p2-2.1.4.tar.gz"
PACELIB_FILENAME="pacelib-2021.10.25.tar.gz"

TARBALLS=(
    MPICH2_WIN64_DEVEL_URL
    MPICH2_WIN32_DEVEL_URL
    VORO_URL
    OPENCL_LOADER_URL
    SCAFACOS_FIX_URL
    GTEST_URL
    YAML_URL
    MATHJAX_URL
    EIGEN3_URL
    CUB_URL
    KOKKOS_URL
    KIM_URL
    MSCG_URL
    PLUMED_URL
    PACELIB_URL
    LATTE_URL
    SCAFACOS_URL
    MDI_URL
    N2P2_URL
)

###############################################################################
# generate proxy cmake file to trick CMake to download from local HTTP server
echo "# auto-generated proxy preset file" > "$HTTP_CACHE_DIR/proxy.cmake"

for t in ${TARBALLS[@]}
do
    FILENAME_VAR="${t/_URL/_FILENAME}"
    if [ -z "${!FILENAME_VAR}" ]
    then
        filename="$(basename ${!t})"
    else
        filename="${!FILENAME_VAR}"
    fi

    if [ ! -f "$HTTP_CACHE_DIR/thirdparty/$filename" ]
    then
        wget -O "$HTTP_CACHE_DIR/thirdparty/$filename" "${!t}"
    fi

    echo "set(${t} \"\${LAMMPS_DOWNLOADS_URL}/thirdparty/$filename\" CACHE STRING \"\" FORCE)" >> "$HTTP_CACHE_DIR/proxy.cmake"
done
