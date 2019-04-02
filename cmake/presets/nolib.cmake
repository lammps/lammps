# preset that turns off all packages that require some form of external
# library or special compiler (fortran or cuda) or equivalent.

set(PACKAGES_WITH_LIB COMPRESS GPU KIM KOKKOS LATTE MPIIO MSCG PYTHON
        VORONOI USER-ADIOS USER-ATC USER-AWPMD USER-H5MD USER-LB
        USER-MOLFILE USER-NETCDF USER-PLUMED USER-QMMM USER-QUIP
        USER-SMD USER-VTK)

foreach(PKG ${PACKAGES_WITH_LIB})
  set(PKG_${PKG} OFF CACHE BOOL "" FORCE)
endforeach()
