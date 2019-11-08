if(PKG_USER-COLVARS)

  set(COLVARS_SOURCE_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars)

  file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/[^.]*.cpp)

  # Build Lepton by default
  set(COLVARS_LEPTON_DEFAULT ON)
  # but not if C++11 is disabled per user request
  if(DEFINED DISABLE_CXX11_REQUIREMENT)
    if(DISABLE_CXX11_REQUIREMENT)
      set(COLVARS_LEPTON_DEFAULT OFF)
    endif()
  endif()

  option(COLVARS_LEPTON "Build and link the Lepton library" ${COLVARS_LEPTON_DEFAULT})

  # Verify that the user's choice is consistent
  if(DEFINED DISABLE_CXX11_REQUIREMENT)
    if((DISABLE_CXX11_REQUIREMENT) AND (COLVARS_LEPTON))
      message(FATAL_ERROR "Building the Lepton library requires C++11 or later.")
    endif()
  endif()

  if(COLVARS_LEPTON)
    set(LEPTON_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars/lepton)
    file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/[^.]*.cpp)
    add_library(lepton STATIC ${LEPTON_SOURCES})
    target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)
  endif()

  add_library(colvars STATIC ${COLVARS_SOURCES})
  target_include_directories(colvars PUBLIC ${LAMMPS_LIB_SOURCE_DIR}/colvars)
  list(APPEND LAMMPS_LINK_LIBS colvars)

  if(COLVARS_LEPTON)
    list(APPEND LAMMPS_LINK_LIBS lepton)
    target_compile_options(colvars PRIVATE -DLEPTON)
    target_include_directories(colvars PUBLIC ${LEPTON_DIR}/include)
  endif()

endif()
