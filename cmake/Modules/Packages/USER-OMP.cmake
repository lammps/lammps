if(PKG_USER-OMP)
    set(USER-OMP_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/USER-OMP)
    set(USER-OMP_SOURCES ${USER-OMP_SOURCES_DIR}/thr_data.cpp
                         ${USER-OMP_SOURCES_DIR}/thr_omp.cpp
                         ${USER-OMP_SOURCES_DIR}/fix_omp.cpp
                         ${USER-OMP_SOURCES_DIR}/fix_nh_omp.cpp
                         ${USER-OMP_SOURCES_DIR}/fix_nh_sphere_omp.cpp
                         ${USER-OMP_SOURCES_DIR}/domain_omp.cpp)
    add_definitions(-DLMP_USER_OMP)
    set_property(GLOBAL PROPERTY "OMP_SOURCES" "${USER-OMP_SOURCES}")

    # detects styles which have USER-OMP version
    RegisterStylesExt(${USER-OMP_SOURCES_DIR} omp OMP_SOURCES)
    RegisterFixStyle(${USER-OMP_SOURCES_DIR}/fix_omp.h)

    get_property(USER-OMP_SOURCES GLOBAL PROPERTY OMP_SOURCES)

    # manually add package dependent source files from USER-OMP that do not provide styles

    if(PKG_ASPHERE)
      list(APPEND USER-OMP_SOURCES ${USER-OMP_SOURCES_DIR}/fix_nh_asphere_omp.cpp)
    endif()

    if(PKG_RIGID)
      list(APPEND USER-OMP_SOURCES ${USER-OMP_SOURCES_DIR}/fix_rigid_nh_omp.cpp)
    endif()

    if(PKG_USER-REAXC)
      list(APPEND USER-OMP_SOURCES ${USER-OMP_SOURCES_DIR}/reaxc_bond_orders_omp.cpp
                                   ${USER-OMP_SOURCES_DIR}/reaxc_hydrogen_bonds_omp.cpp
                                   ${USER-OMP_SOURCES_DIR}/reaxc_nonbonded_omp.cpp
                                   ${USER-OMP_SOURCES_DIR}/reaxc_bonds_omp.cpp
                                   ${USER-OMP_SOURCES_DIR}/reaxc_init_md_omp.cpp
                                   ${USER-OMP_SOURCES_DIR}/reaxc_torsion_angles_omp.cpp
                                   ${USER-OMP_SOURCES_DIR}/reaxc_forces_omp.cpp
                                   ${USER-OMP_SOURCES_DIR}/reaxc_multi_body_omp.cpp
                                   ${USER-OMP_SOURCES_DIR}/reaxc_valence_angles_omp.cpp)
    endif()

    list(APPEND LIB_SOURCES ${USER-OMP_SOURCES})
    include_directories(${USER-OMP_SOURCES_DIR})
endif()
