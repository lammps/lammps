enable_language(C)

# we don't use the parallel i/o interface.
set(HDF5_PREFER_PARALLEL FALSE)

find_package(HDF5 REQUIRED)

# parallel HDF5 will import incompatible MPI headers with a serial build
if((NOT BUILD_MPI) AND HDF5_IS_PARALLEL)
  message(FATAL_ERROR "Serial LAMMPS build and parallel HDF5 library are not compatible")
endif()

target_link_libraries(h5md PRIVATE ${HDF5_LIBRARIES})
target_include_directories(h5md PUBLIC ${HDF5_INCLUDE_DIRS})
