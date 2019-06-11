if(PKG_KOKKOS)
  set(LAMMPS_LIB_KOKKOS_SRC_DIR ${LAMMPS_LIB_SOURCE_DIR}/kokkos)
  set(LAMMPS_LIB_KOKKOS_BIN_DIR ${LAMMPS_LIB_BINARY_DIR}/kokkos)
  add_definitions(-DLMP_KOKKOS)
  add_subdirectory(${LAMMPS_LIB_KOKKOS_SRC_DIR} ${LAMMPS_LIB_KOKKOS_BIN_DIR})

  set(Kokkos_INCLUDE_DIRS ${LAMMPS_LIB_KOKKOS_SRC_DIR}/core/src
                          ${LAMMPS_LIB_KOKKOS_SRC_DIR}/containers/src
                          ${LAMMPS_LIB_KOKKOS_SRC_DIR}/algorithms/src
                          ${LAMMPS_LIB_KOKKOS_BIN_DIR})
  include_directories(${Kokkos_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS kokkos)

  set(KOKKOS_PKG_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/KOKKOS)
  set(KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/atom_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/atom_vec_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/comm_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/comm_tiled_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/neighbor_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/neigh_list_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/neigh_bond_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/fix_nh_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/nbin_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/npair_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/domain_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/modify_kokkos.cpp)

  if(PKG_KSPACE)
    list(APPEND KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/gridcomm_kokkos.cpp)
  endif()

  set_property(GLOBAL PROPERTY "KOKKOS_PKG_SOURCES" "${KOKKOS_PKG_SOURCES}")

  # detects styles which have KOKKOS version
  RegisterStylesExt(${KOKKOS_PKG_SOURCES_DIR} kokkos KOKKOS_PKG_SOURCES)

  # register kokkos-only styles
  RegisterNBinStyle(${KOKKOS_PKG_SOURCES_DIR}/nbin_kokkos.h)
  RegisterNPairStyle(${KOKKOS_PKG_SOURCES_DIR}/npair_kokkos.h)

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
