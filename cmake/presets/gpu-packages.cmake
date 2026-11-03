# preset that turns on all packages that contain styles with GPU support.
# enable the GPU package itself and select the GPU backend(s), e.g.:
#
#   cmake -C ../cmake/presets/gcc.cmake  -D PKG_GPU=ON -D GPU_API=cuda -D GPU_PREC=mixed \
#         -C ../cmake/presets/gpu-packages.cmake ../cmake
#
# The package list below can be regenerated with the following shell command:
#   cd src/GPU ; for f in *_gpu.h; do b=${f/_gpu/}; \
#     [ -f ../$b ] || ls ../*/$b 2> /dev/null; done | cut -d/ -f2 | sort -u

set(ALL_PACKAGES
  AMOEBA
  ASPHERE
  CG-SPICA
  CLASS2
  COLLOID
  CORESHELL
  DIPOLE
  DPD-BASIC
  DPD-MESO
  EXTRA-PAIR
  FEP
  KSPACE
  MANYBODY
  MOLECULE
  SPH)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
