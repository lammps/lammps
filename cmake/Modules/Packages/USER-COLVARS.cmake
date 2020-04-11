set(COLVARS_SOURCE_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars)

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/[^.]*.cpp)

# Build Lepton by default
option(COLVARS_LEPTON "Build and link the Lepton library" ON)

if(COLVARS_LEPTON)
  set(LEPTON_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars/lepton)
  file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/[^.]*.cpp)
  add_library(lepton STATIC ${LEPTON_SOURCES})
  if(BUILD_LIB AND NOT BUILD_SHARED_LIBS)
    install(TARGETS lepton EXPORT LAMMPS_Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()
  set_target_properties(lepton PROPERTIES OUTPUT_NAME lammps_lepton${LAMMPS_LIB_SUFFIX})
  target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)
endif()

add_library(colvars STATIC ${COLVARS_SOURCES})
if(BUILD_LIB AND NOT BUILD_SHARED_LIBS)
  install(TARGETS colvars EXPORT LAMMPS_Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_compile_definitions(colvars PRIVATE -DLAMMPS_${LAMMPS_SIZES})
set_target_properties(colvars PROPERTIES OUTPUT_NAME lammps_colvars${LAMMPS_LIB_SUFFIX})
target_include_directories(colvars PUBLIC ${LAMMPS_LIB_SOURCE_DIR}/colvars)
target_link_libraries(lammps PRIVATE colvars)

if(COLVARS_LEPTON)
  target_link_libraries(lammps PRIVATE lepton)
  target_compile_options(colvars PRIVATE -DLEPTON)
  target_include_directories(colvars PUBLIC ${LEPTON_DIR}/include)
endif()
