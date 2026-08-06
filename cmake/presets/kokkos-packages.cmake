# preset that turns on all packages that contain styles with KOKKOS support.
# It is intended to be combined with one of the kokkos-*.cmake presets that
# enable the KOKKOS package itself and select the KOKKOS backend(s), e.g.:
#
#   cmake -C ../cmake/presets/gcc.cmake -C ../cmake/presets/kokkos-openmp.cmake \
#         -C ../cmake/presets/kokkos-packages.cmake ../cmake
#
# The ML-IAP package also has KOKKOS support but requires the KOKKOS package
# compiled for double precision (-D KOKKOS_PREC=double), so it is not enabled here; add
# -D PKG_ML-IAP=on explicitly for double precision builds.
#
# The package list below can be regenerated with the following shell command:
#   cd src/KOKKOS ; for f in *_kokkos.h; do b=${f/_kokkos/}; \
#     [ -f ../$b ] || ls ../*/$b 2> /dev/null; done | cut -d/ -f2 | sort -u

set(ALL_PACKAGES
  ASPHERE
  CG-SPICA
  CLASS2
  COLLOID
  COLVARS
  DIPOLE
  DPD-BASIC
  DPD-REACT
  EXTRA-COMPUTE
  EXTRA-FIX
  EXTRA-MOLECULE
  EXTRA-PAIR
  GRANULAR
  INTERLAYER
  KSPACE
  MANYBODY
  MEAM
  ML-PACE
  ML-POD
  ML-SNAP
  ML-UF3
  MOFFF
  MOLECULE
  PHONON
  REAXFF
  RIGID
  SPIN
  YAFF)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
