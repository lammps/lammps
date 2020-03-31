if(PKG_USER-NETCDF)
  # USER-NETCDF can use NetCDF, Parallel NetCDF (PNetCDF), or both. At least one necessary.
  # NetCDF library enables dump style "netcdf", while PNetCDF enables dump style "netcdf/mpiio"
  find_package(NetCDF)
  if(NETCDF_FOUND)
    find_package(PNetCDF)
  else(NETCDF_FOUND)
    find_package(PNetCDF REQUIRED)
  endif(NETCDF_FOUND)

  if(NETCDF_FOUND)
    include_directories(${NETCDF_INCLUDE_DIRS})
    list(APPEND LAMMPS_LINK_LIBS ${NETCDF_LIBRARIES})
    add_definitions(-DLMP_HAS_NETCDF)
  endif(NETCDF_FOUND)

  if(PNETCDF_FOUND)
    include_directories(${PNETCDF_INCLUDES})
    list(APPEND LAMMPS_LINK_LIBS ${PNETCDF_LIBRARIES})
    add_definitions(-DLMP_HAS_PNETCDF)
  endif(PNETCDF_FOUND)

  add_definitions(-DNC_64BIT_DATA=0x0020)
endif()
