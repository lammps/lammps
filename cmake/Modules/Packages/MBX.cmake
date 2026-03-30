# MBX v1.3.0+ support for MBX package
# Based off of PLUMED.cmake

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

if((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND (NOT CMAKE_CROSSCOMPILING))
  message(FATAL_ERROR "Compiling the MBX package natively for Windows is currently not supported")
endif()

# for supporting multiple concurrent mbx installations for debugging and testing
set(MBX_SUFFIX "" CACHE STRING "Suffix for MBX library")
mark_as_advanced(MBX_SUFFIX)

set(MBX_CONFIG_CC  ${CMAKE_C_COMPILER})
set(MBX_CONFIG_CXX  ${CMAKE_CXX_COMPILER})
if(BUILD_MPI)
  set(MBX_CONFIG_FLAGS "--enable-mpi")
  set(MBX_CONFIG_CC  ${MPI_C_COMPILER})
  set(MBX_CONFIG_CXX  ${MPI_CXX_COMPILER})
else()
  set(MBX_CONFIG_FLAGS "--disable-mpi")
endif()

if(BUILD_OMP)
  list(APPEND MBX_CONFIG_FLAGS --enable-openmp)
else()
  list(APPEND MBX_CONFIG_FLAGS --disable-openmp)
endif()

if(CONFIGURE_REQUEST_PIC)
  list(APPEND MBX_CONFIG_FLAGS ${CONFIGURE_REQUEST_PIC})
endif()

set(MBXLIB_URL "https://github.com/paesanilab/MBX/releases/download/v1.3.12/mbx-1.3.12.tar.gz" CACHE STRING "URL for MBX tarball")
set(MBXLIB_SHA256 "0f3600d0841e1736abec3fe6a38c8ff673523fe4494169db2839ddaa9bb618bf" CACHE STRING "SHA256 checksum of MBX tarball")

mark_as_advanced(MBXLIB_URL)
mark_as_advanced(MBXLIB_SHA256)

set(MBX_LINK_LIBS)
find_package(FFTW3 REQUIRED)
if(FFTW3_FOUND)
  list(APPEND MBX_LINK_LIBS FFTW3::FFTW3)
endif()

find_package(PkgConfig QUIET)
set(DOWNLOAD_MBX_DEFAULT ON)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(MBX QUIET mbx)
  if(MBX_FOUND)
    set(DOWNLOAD_MBX_DEFAULT OFF)
  endif()
endif()

option(DOWNLOAD_MBX "Download MBX package instead of using an already installed one" ${DOWNLOAD_MBX_DEFAULT})
if(DOWNLOAD_MBX)
  message(STATUS "MBX download requested - we will build our own")
  set(MBX_BUILD_BYPRODUCTS "<INSTALL_DIR>/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mbx${CMAKE_STATIC_LIBRARY_SUFFIX}")

  message(STATUS "MBX_CONFIG_CC: ${MBX_CONFIG_CC}")
  message(STATUS "MBX_CONFIG_CXX: ${MBX_CONFIG_CXX}")
  message(STATUS "MBX_CONFIG_FLAGS: ${MBX_CONFIG_FLAGS}")

  include(ExternalProject)

  # hacks to make Linux-to-Windows cross-compilation work
  if((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND (CMAKE_CROSSCOMPILING))
    ExternalProject_Add(mbx_build
      URL     ${MBXLIB_URL}
      URL_HASH SHA256=${MBXLIB_SHA256}
      BUILD_IN_SOURCE TRUE
      CONFIGURE_COMMAND mingw64-configure
                        --prefix=<INSTALL_DIR>
                        --disable-i-pi-plugin
                        ${MBX_CONFIG_FLAGS}
      INSTALL_COMMAND make install prefix=-build includedir=/include bindir=/bin
                                 datadir=/share exec_prefix=/ infodir=/share/info
                                 libdir=/lib libexecdir=/libexec localstatedir=/var
                                 mandir=/share/man sbindir=/sbin sharedstatedir=/com
                                 sysconfdir=/etc DESTDIR=<INSTALL_DIR>
      BUILD_BYPRODUCTS ${MBX_BUILD_BYPRODUCTS}
    )
  else()
    ExternalProject_Add(mbx_build
      URL     ${MBXLIB_URL}
      URL_HASH SHA256=${MBXLIB_SHA256}
      CONFIGURE_COMMAND <SOURCE_DIR>/configure
                        --prefix=<INSTALL_DIR>
                        --disable-i-pi-plugin
                        ${MBX_CONFIG_FLAGS}
                        CXX=${MBX_CONFIG_CXX}
                        CC=${MBX_CONFIG_CC}
                        CPPFLAGS=-I${FFTW3_INCLUDE_DIRS}
      BUILD_BYPRODUCTS ${MBX_BUILD_BYPRODUCTS}
    )
  endif()
  ExternalProject_get_property(mbx_build INSTALL_DIR)
  add_library(LAMMPS::MBX UNKNOWN IMPORTED)
  add_dependencies(LAMMPS::MBX mbx_build)
  set_target_properties(LAMMPS::MBX PROPERTIES
    IMPORTED_LOCATION ${INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mbx${CMAKE_STATIC_LIBRARY_SUFFIX}
    INTERFACE_LINK_LIBRARIES "${MBX_LINK_LIBS};${CMAKE_DL_LIBS}"
  )

  set_target_properties(LAMMPS::MBX PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)
  file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
  if(CMAKE_PROJECT_NAME STREQUAL "lammps")
    target_link_libraries(lammps PRIVATE LAMMPS::MBX)
  endif()
else()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(MBX REQUIRED mbx${MBX_SUFFIX})
  add_library(LAMMPS::MBX INTERFACE IMPORTED)

  file(MAKE_DIRECTORY ${MBX_INCLUDE_DIRS})

  target_include_directories(LAMMPS::MBX INTERFACE ${MBX_INCLUDE_DIRS})
  target_link_directories(LAMMPS::MBX INTERFACE ${MBX_LIBRARY_DIRS})
  target_link_libraries(LAMMPS::MBX INTERFACE "${MBX_LINK_LIBS};${MBX_LIBRARIES};${CMAKE_DL_LIBS}")

  if(CMAKE_PROJECT_NAME STREQUAL "lammps")
    target_link_libraries(lammps PUBLIC LAMMPS::MBX)
  endif()
endif()

