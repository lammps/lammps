# NETCDF can use NetCDF, Parallel NetCDF (PNetCDF), or both. At least one necessary.
# NetCDF library enables dump style "netcdf", while PNetCDF enables dump style "netcdf/mpiio"

# may use NetCDF or PNetCDF with MPI, but must have NetCDF without
if(NOT BUILD_MPI)
  find_package(NetCDF REQUIRED)
else()
  find_package(NetCDF)
  if(NETCDF_FOUND)
    find_package(PNetCDF)
  else()
    find_package(PNetCDF REQUIRED)
  endif()
endif()

if(NETCDF_FOUND)
  target_link_libraries(lammps PRIVATE NetCDF::NetCDF)
  target_compile_definitions(lammps PRIVATE -DLMP_HAS_NETCDF)
endif(NETCDF_FOUND)

if(PNETCDF_FOUND)
  target_link_libraries(lammps PRIVATE PNetCDF::PNetCDF)
  target_compile_definitions(lammps PRIVATE -DLMP_HAS_PNETCDF)
endif(PNETCDF_FOUND)

target_compile_definitions(lammps PRIVATE -DNC_64BIT_DATA=0x0020)
