  set(OPENMP_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/OPENMP)
  set(OPENMP_SOURCES ${OPENMP_SOURCES_DIR}/thr_data.cpp
                       ${OPENMP_SOURCES_DIR}/thr_omp.cpp
                       ${OPENMP_SOURCES_DIR}/fix_omp.cpp
                       ${OPENMP_SOURCES_DIR}/fix_nh_omp.cpp
                       ${OPENMP_SOURCES_DIR}/fix_nh_sphere_omp.cpp
                       ${OPENMP_SOURCES_DIR}/domain_omp.cpp)
  target_compile_definitions(lammps PRIVATE -DLMP_USER_OMP)
  set_property(GLOBAL PROPERTY "OMP_SOURCES" "${OPENMP_SOURCES}")

  # detects styles which have OPENMP version
  RegisterStylesExt(${OPENMP_SOURCES_DIR} omp OMP_SOURCES)
  RegisterFixStyle(${OPENMP_SOURCES_DIR}/fix_omp.h)

  get_property(OPENMP_SOURCES GLOBAL PROPERTY OMP_SOURCES)

  # manually add package dependent source files from OPENMP that do not provide styles

  if(PKG_ASPHERE)
    list(APPEND OPENMP_SOURCES ${OPENMP_SOURCES_DIR}/fix_nh_asphere_omp.cpp)
  endif()

  if(PKG_RIGID)
    list(APPEND OPENMP_SOURCES ${OPENMP_SOURCES_DIR}/fix_rigid_nh_omp.cpp)
  endif()

  if(PKG_USER-REAXC)
    list(APPEND OPENMP_SOURCES ${OPENMP_SOURCES_DIR}/reaxc_bond_orders_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxc_hydrogen_bonds_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxc_nonbonded_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxc_bonds_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxc_init_md_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxc_torsion_angles_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxc_forces_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxc_multi_body_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxc_valence_angles_omp.cpp)
  endif()

  target_sources(lammps PRIVATE ${OPENMP_SOURCES})
  target_include_directories(lammps PRIVATE ${OPENMP_SOURCES_DIR})
