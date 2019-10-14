if(PKG_USER-COLVARS)

  set(COLVARS_SOURCE_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars)

  file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/[^.]*.cpp)

  if(DEFINED CMAKE_CXX_STANDARD)
    # The condition below should match most modern C++ standards
    if((${CMAKE_CXX_STANDARD} GREATER 0) AND (${CMAKE_CXX_STANDARD} LESS 70))
      set(ENABLE_LEPTON ON)
      set(LEPTON_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars/lepton)
      file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/[^.]*.cpp)
      add_library(lepton STATIC ${LEPTON_SOURCES})
      target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)
    endif()
  endif()

  add_library(colvars STATIC ${COLVARS_SOURCES})
  target_include_directories(colvars PUBLIC ${LAMMPS_LIB_SOURCE_DIR}/colvars)
  list(APPEND LAMMPS_LINK_LIBS colvars)

  if(ENABLE_LEPTON)
    list(APPEND LAMMPS_LINK_LIBS lepton)
    target_compile_options(colvars PRIVATE -DLEPTON)
    target_include_directories(colvars PUBLIC ${LEPTON_DIR}/include)
  endif()

endif()
