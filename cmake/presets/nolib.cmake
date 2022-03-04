# preset that turns off all packages that require some form of external
# library or special compiler (fortran or cuda) or equivalent.

set(PACKAGES_WITH_LIB
  ADIOS
  ATC
  AWPMD
  COMPRESS
  GPU
  H5MD
  KIM
  KOKKOS
  LATBOLTZ
  LATTE
  MACHDYN
  MDI
  MESONT
  MESSAGE
  ML-HDNNP
  ML-PACE
  ML-QUIP
  MOLFILE
  MPIIO
  MSCG
  NETCDF
  PLUMED
  PYTHON
  QMMM
  SCAFACOS
  VORONOI
  VTK)

foreach(PKG ${PACKAGES_WITH_LIB})
  set(PKG_${PKG} OFF CACHE BOOL "" FORCE)
endforeach()
