if(PKG_KOKKOS)
  find_package(Kokkos 3 QUIET)
  if(Kokkos_FOUND)
    set(DOWNLOAD_KOKKOS_DEFAULT OFF)
  else()
    set(DOWNLOAD_KOKKOS_DEFAULT ON)
  endif()
  option(DOWNLOAD_KOKKOS "Download the KOKKOS library instead of using an already installed one" ${DOWNLOAD_KOKKOS_DEFAULT})
  if(DOWNLOAD_KOKKOS)
    message(STATUS "KOKKOS download requested - we will build our own")
    # Workaround for cross compilation with MinGW where ${CMAKE_INSTALL_LIBDIR}
    # is a full path, so we need to remove the prefix
    string(REPLACE ${CMAKE_INSTALL_PREFIX} "" _KOKKOS_LIBDIR ${CMAKE_INSTALL_LIBDIR})
    include(ExternalProject)
    ExternalProject_Add(kokkos_build
      URL https://github.com/kokkos/kokkos/archive/3.0.00.tar.gz
      URL_MD5 281c7093aa3a603276e93abdf4be23b9
      PATCH_COMMAND patch -p3 < ${CMAKE_CURRENT_SOURCE_DIR}/2826.patch
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ${CMAKE_REQUEST_PIC}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM} -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
      BUILD_BYPRODUCTS <INSTALL_DIR>/${_KOKKOS_LIBDIR}/libkokkoscore.a
    )
    list(APPEND LAMMPS_DEPS kokkos_build)
    ExternalProject_get_property(kokkos_build INSTALL_DIR)
    set(KOKKOS_LIBRARIES ${INSTALL_DIR}/${_KOKKOS_LIBDIR}/libkokkoscore.a) 
    set(KOKKOS_INCLUDE_DIRS ${INSTALL_DIR}/include)
    include_directories(${KOKKOS_INCLUDE_DIRS})
    list(APPEND LAMMPS_LINK_LIBS ${KOKKOS_LIBRARIES})
  else()
    find_package(Kokkos 3)
    if(NOT Kokkos_FOUND)
      message(FATAL_ERROR "KOKKOS library not found, help CMake to find it by setting KOKKOS_LIBRARY, or set DOWNLOAD_KOKKOS=ON to download it")
    endif()
    list(APPEND LAMMPS_LINK_LIBS Kokkos::kokkos)
  endif()
  add_definitions(-DLMP_KOKKOS)

  set(KOKKOS_PKG_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/KOKKOS)
  set(KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/atom_kokkos.cpp
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
    if(KOKKOS_ENABLE_CUDA)
      if(NOT ${FFT} STREQUAL "KISS")
        add_definitions(-DFFT_CUFFT)
        list(APPEND LAMMPS_LINK_LIBS cufft)
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

  if(PKG_USER-DPD)
    get_property(KOKKOS_PKG_SOURCES GLOBAL PROPERTY KOKKOS_PKG_SOURCES)
    list(APPEND KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/npair_ssa_kokkos.cpp)
    RegisterNPairStyle(${KOKKOS_PKG_SOURCES_DIR}/npair_ssa_kokkos.h)
    set_property(GLOBAL PROPERTY "KOKKOS_PKG_SOURCES" "${KOKKOS_PKG_SOURCES}")
  endif()

  get_property(KOKKOS_PKG_SOURCES GLOBAL PROPERTY KOKKOS_PKG_SOURCES)

  list(APPEND LIB_SOURCES ${KOKKOS_PKG_SOURCES})
  include_directories(${KOKKOS_PKG_SOURCES_DIR})
endif()
