# Download and configure custom MPICH files for Windows
message(STATUS "Downloading and configuring MPICH-1.4.1 for Windows")
include(ExternalProject)
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  ExternalProject_Add(mpi4win_build
    URL https://download.lammps.org/thirdparty/mpich2-win64-devel.tar.gz
    URL_MD5 4939fdb59d13182fd5dd65211e469f14
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND ""
    BUILD_BYPRODUCTS <SOURCE_DIR>/lib/libmpi.a)
else()
  ExternalProject_Add(mpi4win_build
    URL https://download.lammps.org/thirdparty/mpich2-win32-devel.tar.gz
    URL_MD5 a61d153500dce44e21b755ee7257e031
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND ""
    BUILD_BYPRODUCTS <SOURCE_DIR>/lib/libmpi.a)
endif()

ExternalProject_get_property(mpi4win_build SOURCE_DIR)
file(MAKE_DIRECTORY "${SOURCE_DIR}/include")
add_library(MPI::MPI_CXX UNKNOWN IMPORTED)
set_target_properties(MPI::MPI_CXX PROPERTIES
  IMPORTED_LOCATION "${SOURCE_DIR}/lib/libmpi.a"
  INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include"
  INTERFACE_COMPILE_DEFINITIONS "MPICH_SKIP_MPICXX")
add_dependencies(MPI::MPI_CXX mpi4win_build)
