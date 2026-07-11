  set(OPENMP_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/OPENMP)
  set(OPENMP_SOURCES ${OPENMP_SOURCES_DIR}/thr_data.cpp
                       ${OPENMP_SOURCES_DIR}/thr_omp.cpp
                       ${OPENMP_SOURCES_DIR}/fix_omp.cpp
                       ${OPENMP_SOURCES_DIR}/fix_nh_omp.cpp
                       ${OPENMP_SOURCES_DIR}/fix_nh_sphere_omp.cpp
                       ${OPENMP_SOURCES_DIR}/domain_omp.cpp)
  target_compile_definitions(lammps PRIVATE -DLMP_OPENMP)
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

  if(PKG_REAXFF)
    list(APPEND OPENMP_SOURCES ${OPENMP_SOURCES_DIR}/reaxff_bond_orders_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxff_hydrogen_bonds_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxff_nonbonded_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxff_bonds_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxff_init_md_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxff_torsion_angles_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxff_forces_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxff_multi_body_omp.cpp
                                 ${OPENMP_SOURCES_DIR}/reaxff_valence_angles_omp.cpp)
  endif()

  target_sources(lammps PRIVATE ${OPENMP_SOURCES})
  target_include_directories(lammps PRIVATE ${OPENMP_SOURCES_DIR})
