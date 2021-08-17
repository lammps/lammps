# Download and configure custom MPICH files for Windows
message(STATUS "Downloading and configuring MPICH-1.4.1 for Windows")
set(MPICH2_WIN64_DEVEL_URL "${LAMMPS_THIRDPARTY_URL}/mpich2-win64-devel.tar.gz" CACHE STRING "URL for MPICH2 (win64) tarball")
set(MPICH2_WIN32_DEVEL_URL "${LAMMPS_THIRDPARTY_URL}/mpich2-win32-devel.tar.gz" CACHE STRING "URL for MPICH2 (win32) tarball")
set(MPICH2_WIN64_DEVEL_MD5 "4939fdb59d13182fd5dd65211e469f14" CACHE STRING "MD5 checksum of MPICH2 (win64) tarball")
set(MPICH2_WIN32_DEVEL_MD5 "a61d153500dce44e21b755ee7257e031" CACHE STRING "MD5 checksum of MPICH2 (win32) tarball")
mark_as_advanced(MPICH2_WIN64_DEVEL_URL)
mark_as_advanced(MPICH2_WIN32_DEVEL_URL)
mark_as_advanced(MPICH2_WIN64_DEVEL_MD5)
mark_as_advanced(MPICH2_WIN32_DEVEL_MD5)

include(ExternalProject)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  ExternalProject_Add(mpi4win_build
    URL     ${MPICH2_WIN64_DEVEL_URL}
    URL_MD5 ${MPICH2_WIN64_DEVEL_MD5}
    CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND ""
    BUILD_BYPRODUCTS <SOURCE_DIR>/lib/libmpi.a)
else()
  ExternalProject_Add(mpi4win_build
    URL     ${MPICH2_WIN32_DEVEL_URL}
    URL_MD5 ${MPICH2_WIN32_DEVEL_MD5}
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

# set variables for status reporting at the end of CMake run
set(MPI_CXX_INCLUDE_PATH "${SOURCE_DIR}/include")
set(MPI_CXX_COMPILE_DEFINITIONS "MPICH_SKIP_MPICXX")
set(MPI_CXX_LIBRARIES "${SOURCE_DIR}/lib/libmpi.a")
