# set-up MS-MPI library for Windows with MinGW compatibility
message(STATUS "Downloading and configuring MS-MPI 10.1 for Windows cross-compilation")
set(MPICH2_WIN64_DEVEL_URL "${LAMMPS_THIRDPARTY_URL}/msmpi-win64-devel.tar.gz" CACHE STRING "URL for MS-MPI (win64) tarball")
set(MPICH2_WIN64_DEVEL_MD5 "86314daf1bffb809f1fcbefb8a547490" CACHE STRING "MD5 checksum of MS-MPI (win64) tarball")
mark_as_advanced(MPICH2_WIN64_DEVEL_URL)
mark_as_advanced(MPICH2_WIN64_DEVEL_MD5)

include(ExternalProject)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  ExternalProject_Add(mpi4win_build
    URL     ${MPICH2_WIN64_DEVEL_URL}
    URL_MD5 ${MPICH2_WIN64_DEVEL_MD5}
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND ""
    BUILD_BYPRODUCTS <SOURCE_DIR>/lib/libmsmpi.a)
else()
  message(FATAL_ERROR "Only x86 64-bit builds are supported with MS-MPI")
endif()

ExternalProject_get_property(mpi4win_build SOURCE_DIR)
file(MAKE_DIRECTORY "${SOURCE_DIR}/include")
add_library(MPI::MPI_CXX UNKNOWN IMPORTED)
set_target_properties(MPI::MPI_CXX PROPERTIES
  IMPORTED_LOCATION "${SOURCE_DIR}/lib/libmsmpi.a"
  INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include"
  INTERFACE_COMPILE_DEFINITIONS "MPICH_SKIP_MPICXX=1")
add_dependencies(MPI::MPI_CXX mpi4win_build)

# set variables for status reporting at the end of CMake run
set(MPI_CXX_INCLUDE_PATH "${SOURCE_DIR}/include")
set(MPI_CXX_COMPILE_DEFINITIONS "MPICH_SKIP_MPICXX=1")
set(MPI_CXX_LIBRARIES "${SOURCE_DIR}/lib/libmsmpi.a")
