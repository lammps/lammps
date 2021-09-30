########################################################################
# As of version 3.3.0 Kokkos requires C++14
if(CMAKE_CXX_STANDARD LESS 14)
  set(CMAKE_CXX_STANDARD 14)
endif()
########################################################################
# consistency checks and Kokkos options/settings required by LAMMPS
if(Kokkos_ENABLE_CUDA)
  message(STATUS "KOKKOS: Enabling CUDA LAMBDA function support")
  set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
endif()
# Adding OpenMP compiler flags without the checks done for
# BUILD_OMP can result in compile failures. Enforce consistency.
if(Kokkos_ENABLE_OPENMP AND NOT BUILD_OMP)
  message(FATAL_ERROR "Must enable BUILD_OMP with Kokkos_ENABLE_OPENMP")
endif()
########################################################################

option(EXTERNAL_KOKKOS "Build against external kokkos library" OFF)
option(DOWNLOAD_KOKKOS "Download the KOKKOS library instead of using the bundled one" OFF)
if(DOWNLOAD_KOKKOS)
  # extract Kokkos-related variables and values so we can forward them to the Kokkos library build
  get_cmake_property(_VARS VARIABLES)
  list(FILTER _VARS INCLUDE REGEX ^Kokkos_)
  foreach(_VAR IN LISTS _VARS)
    list(APPEND KOKKOS_LIB_BUILD_ARGS "-D${_VAR}=${${_VAR}}")
  endforeach()
  message(STATUS "KOKKOS download requested - we will build our own")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
  if(CMAKE_REQUEST_PIC)
    list(APPEND KOKKOS_LIB_BUILD_ARGS ${CMAKE_REQUEST_PIC})
  endif()
  # append other CMake variables that need to be forwarded to CMAKE_ARGS
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_INSTALL_LIBDIR=lib")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}")
  list(APPEND KOKKOS_LIB_BUILD_ARGS "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
  include(ExternalProject)
  set(KOKKOS_URL "https://github.com/kokkos/kokkos/archive/3.4.01.tar.gz" CACHE STRING "URL for KOKKOS tarball")
  set(KOKKOS_MD5 "4c84698917c93a18985b311bb6caf84f" CACHE STRING "MD5 checksum of KOKKOS tarball")
  mark_as_advanced(KOKKOS_URL)
  mark_as_advanced(KOKKOS_MD5)
  ExternalProject_Add(kokkos_build
    URL     ${KOKKOS_URL}
    URL_MD5 ${KOKKOS_MD5}
    CMAKE_ARGS ${KOKKOS_LIB_BUILD_ARGS}
    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/libkokkoscore.a
  )
  ExternalProject_get_property(kokkos_build INSTALL_DIR)
  file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
  add_library(LAMMPS::KOKKOS UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::KOKKOS PROPERTIES
    IMPORTED_LOCATION "${INSTALL_DIR}/lib/libkokkoscore.a"
    INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/include"
    INTERFACE_LINK_LIBRARIES ${CMAKE_DL_LIBS})
  target_link_libraries(lammps PRIVATE LAMMPS::KOKKOS)
  target_link_libraries(lmp PRIVATE LAMMPS::KOKKOS)
  add_dependencies(LAMMPS::KOKKOS kokkos_build)
elseif(EXTERNAL_KOKKOS)
  find_package(Kokkos 3.4.01 REQUIRED CONFIG)
  target_link_libraries(lammps PRIVATE Kokkos::kokkos)
  target_link_libraries(lmp PRIVATE Kokkos::kokkos)
else()
  set(LAMMPS_LIB_KOKKOS_SRC_DIR ${LAMMPS_LIB_SOURCE_DIR}/kokkos)
  set(LAMMPS_LIB_KOKKOS_BIN_DIR ${LAMMPS_LIB_BINARY_DIR}/kokkos)
  add_subdirectory(${LAMMPS_LIB_KOKKOS_SRC_DIR} ${LAMMPS_LIB_KOKKOS_BIN_DIR})

  set(Kokkos_INCLUDE_DIRS ${LAMMPS_LIB_KOKKOS_SRC_DIR}/core/src
                          ${LAMMPS_LIB_KOKKOS_SRC_DIR}/containers/src
                          ${LAMMPS_LIB_KOKKOS_SRC_DIR}/algorithms/src
                          ${LAMMPS_LIB_KOKKOS_BIN_DIR})
  target_include_directories(lammps PRIVATE ${Kokkos_INCLUDE_DIRS})
  target_link_libraries(lammps PRIVATE kokkos)
  target_link_libraries(lmp PRIVATE kokkos)
endif()
target_compile_definitions(lammps PUBLIC $<BUILD_INTERFACE:LMP_KOKKOS>)

set(KOKKOS_PKG_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/KOKKOS)
set(KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/atom_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/atom_map_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/atom_vec_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/comm_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/comm_tiled_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/min_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/min_linesearch_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/neighbor_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/neigh_list_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/neigh_bond_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/fix_nh_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/nbin_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/npair_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/npair_halffull_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/domain_kokkos.cpp
                       ${KOKKOS_PKG_SOURCES_DIR}/modify_kokkos.cpp)

if(PKG_KSPACE)
  list(APPEND KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/fft3d_kokkos.cpp
                                 ${KOKKOS_PKG_SOURCES_DIR}/gridcomm_kokkos.cpp
                                 ${KOKKOS_PKG_SOURCES_DIR}/remap_kokkos.cpp)
  if(Kokkos_ENABLE_CUDA)
    if(NOT (FFT STREQUAL "KISS"))
      target_compile_definitions(lammps PRIVATE -DFFT_CUFFT)
      target_link_libraries(lammps PRIVATE cufft)
    endif()
  endif()
endif()

set_property(GLOBAL PROPERTY "KOKKOS_PKG_SOURCES" "${KOKKOS_PKG_SOURCES}")

# detects styles which have KOKKOS version
RegisterStylesExt(${KOKKOS_PKG_SOURCES_DIR} kokkos KOKKOS_PKG_SOURCES)

# register kokkos-only styles
RegisterNBinStyle(${KOKKOS_PKG_SOURCES_DIR}/nbin_kokkos.h)
RegisterNPairStyle(${KOKKOS_PKG_SOURCES_DIR}/npair_kokkos.h)
RegisterNPairStyle(${KOKKOS_PKG_SOURCES_DIR}/npair_halffull_kokkos.h)

if(PKG_DPD-REACT)
  get_property(KOKKOS_PKG_SOURCES GLOBAL PROPERTY KOKKOS_PKG_SOURCES)
  list(APPEND KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/npair_ssa_kokkos.cpp)
  RegisterNPairStyle(${KOKKOS_PKG_SOURCES_DIR}/npair_ssa_kokkos.h)
  set_property(GLOBAL PROPERTY "KOKKOS_PKG_SOURCES" "${KOKKOS_PKG_SOURCES}")
endif()

get_property(KOKKOS_PKG_SOURCES GLOBAL PROPERTY KOKKOS_PKG_SOURCES)

target_sources(lammps PRIVATE ${KOKKOS_PKG_SOURCES})
target_include_directories(lammps PUBLIC $<BUILD_INTERFACE:${KOKKOS_PKG_SOURCES_DIR}>)
