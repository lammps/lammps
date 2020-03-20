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
add_definitions(-DMPICH_SKIP_MPICXX)
include_directories("${SOURCE_DIR}/include")
set(MPI4WIN_LIBRARIES "${SOURCE_DIR}/lib/libmpi.a")
list(APPEND LAMMPS_DEPS mpi4win_build)
set(LAMMPS_USE_MPI4WIN ON)
