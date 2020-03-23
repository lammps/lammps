if(PKG_USER-COLVARS)

  set(COLVARS_SOURCE_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars)

  file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/[^.]*.cpp)

  # Build Lepton by default
  option(COLVARS_LEPTON "Build and link the Lepton library" ON)

  if(COLVARS_LEPTON)
    set(LEPTON_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars/lepton)
    file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/[^.]*.cpp)
    add_library(lepton STATIC ${LEPTON_SOURCES})
    target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)
  endif()

  add_library(colvars STATIC ${COLVARS_SOURCES})
  target_include_directories(colvars PUBLIC ${LAMMPS_LIB_SOURCE_DIR}/colvars)
  target_link_libraries(lammps PRIVATE colvars)

  if(COLVARS_LEPTON)
    target_link_libraries(lammps PRIVATE lepton)
    target_compile_options(colvars PRIVATE -DLEPTON)
    target_include_directories(colvars PUBLIC ${LEPTON_DIR}/include)
  endif()

endif()
