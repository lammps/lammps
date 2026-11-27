# MBX v1.3.0+ support for MBX package
# Based off of PLUMED.cmake

# set policy to use the time of extraction as timestamps of files unpacked from downloaded
# archives, so that updating an archive version triggers rebuilding all dependent objects
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
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
  if(MPI_C_COMPILER)
    set(MPI_CONFIG_CC  ${MPI_C_COMPILER})
  endif()
  if(MPI_CXX_COMPILER)
    set(MBX_CONFIG_CXX  ${MPI_CXX_COMPILER})
  endif()
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
  SetDownloadSettings(MBXLIB "MBX"
    "https://github.com/paesanilab/MBX/releases/download/v1.4.0/mbx-1.4.0.tar.gz"
    "219de5af7cd81bdf8c3394c3aeee923cf4de7dd322769441a35c7cf4950bf912")

  # apply patches to the MBX sources of known versions (e.g. so that they compile without OpenMP)
  get_filename_component(_mbx_archive "${MBXLIB_URL}" NAME)
  string(REGEX REPLACE "\\.(tar\\.gz|tgz|tar\\.xz|txz|zip)$" "" _mbx_version "${_mbx_archive}")
  set(MBX_PATCH_FILE ${LAMMPS_DIR}/cmake/patches/${_mbx_version}-openmp.patch)
  set(MBX_PATCH_COMMAND)
  if(EXISTS ${MBX_PATCH_FILE})
    find_package(Patch REQUIRED)
    message(STATUS "Patching MBX sources with ${MBX_PATCH_FILE}")
    set(MBX_PATCH_COMMAND COMMAND ${Patch_EXECUTABLE} -p1 -i ${MBX_PATCH_FILE})
  else()
    message(STATUS "No patch file ${MBX_PATCH_FILE} found for MBX archive ${_mbx_archive}")
  endif()
  set(MBX_BUILD_BYPRODUCTS "<INSTALL_DIR>/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mbx${CMAKE_STATIC_LIBRARY_SUFFIX}")

  message(STATUS "MBX_CONFIG_CC: ${MBX_CONFIG_CC}")
  message(STATUS "MBX_CONFIG_CXX: ${MBX_CONFIG_CXX}")
  message(STATUS "MBX_CONFIG_FLAGS: ${MBX_CONFIG_FLAGS}")

  include(ExternalProject)

  # hacks to make testing Linux-to-Windows cross-compilation work.  Compiling some MBX sources
  # with optimization takes a long time, so without MPI a precompiled MBX is used where one
  # exists for this MBX version and the major version of the cross-compiler: an archive made with
  # cmake/packaging/build-mbx-mingw.sh that is only unpacked.  MBX_MINGW_SHA256_<name> below lists
  # those on download.lammps.org.  MBX_MINGW_FROM_SOURCE=on compiles MBX anyway (to make a new one),
  # and MBX_MINGW_PRECOMPILED tells cmake/packaging/cmake-win-on-linux.py which of the two happened.
  if((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND (CMAKE_CROSSCOMPILING))
    # <MBX version>-mingw64-gcc<major version>, e.g. set(MBX_MINGW_SHA256_mbx-1.4.0-mingw64-gcc16 "<SHA256>")
    option(MBX_MINGW_FROM_SOURCE "Compile MBX even where a precompiled one exists" OFF)
    mark_as_advanced(MBX_MINGW_FROM_SOURCE)
    string(REGEX MATCH "^[0-9]+" _mbx_gcc "${CMAKE_CXX_COMPILER_VERSION}")
    set(MBX_MINGW_ARCHIVE "${_mbx_version}-mingw64-gcc${_mbx_gcc}" CACHE INTERNAL "Precompiled MBX for this MBX version and compiler")
    set(MBX_MINGW_PRECOMPILED OFF CACHE INTERNAL "MBX unpacked from a precompiled archive")
    if((NOT BUILD_MPI) AND (NOT MBX_MINGW_FROM_SOURCE) AND (DEFINED MBX_MINGW_SHA256_${MBX_MINGW_ARCHIVE}))
      SetDownloadSettings(MBXMINGW "precompiled MBX"
        "https://download.lammps.org/thirdparty/${MBX_MINGW_ARCHIVE}.tar.gz"
        "${MBX_MINGW_SHA256_${MBX_MINGW_ARCHIVE}}")
      message(STATUS "Using the precompiled MBX ${MBX_MINGW_ARCHIVE}")
      ExternalProject_Add(mbx_build
        URL     ${MBXMINGW_URL}
        URL_HASH SHA256=${MBXMINGW_SHA256}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR> <INSTALL_DIR>
        BUILD_BYPRODUCTS ${MBX_BUILD_BYPRODUCTS}
      )
      set(MBX_MINGW_PRECOMPILED ON CACHE INTERNAL "MBX unpacked from a precompiled archive")
    else()
      ExternalProject_Add(mbx_build
        URL     ${MBXLIB_URL}
        URL_HASH SHA256=${MBXLIB_SHA256}
        BUILD_IN_SOURCE TRUE
        # restore the timestamp order that autotools generated files require (see AutotoolsTouch.cmake)
        PATCH_COMMAND ${CMAKE_COMMAND} -D SOURCE_DIR=<SOURCE_DIR> -P ${LAMMPS_DIR}/cmake/Modules/AutotoolsTouch.cmake
                      ${MBX_PATCH_COMMAND}
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
    endif()
  else()
    ExternalProject_Add(mbx_build
      URL     ${MBXLIB_URL}
      URL_HASH SHA256=${MBXLIB_SHA256}
      # restore the timestamp order that autotools generated files require (see AutotoolsTouch.cmake)
      PATCH_COMMAND ${CMAKE_COMMAND} -D SOURCE_DIR=<SOURCE_DIR> -P ${LAMMPS_DIR}/cmake/Modules/AutotoolsTouch.cmake
                    ${MBX_PATCH_COMMAND}
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

