# preset that turns off all packages that require some form of external
# library or special compiler (fortran or cuda) or equivalent.

set(PACKAGES_WITH_LIB COMPRESS GPU KIM KOKKOS LATTE MESSAGE MPIIO MSCG
  PYTHON VORONOI
  ADIOS ATC AWPMD H5MD ML-HDNNP LATBOLTZ MOLFILE
  MESONT MDI NETCDF ML-PACE PLUMED QMMM ML-QUIP
  SCAFACOS MACHDYN VTK)

foreach(PKG ${PACKAGES_WITH_LIB})
  set(PKG_${PKG} OFF CACHE BOOL "" FORCE)
endforeach()
