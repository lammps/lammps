set(COLVARS_SOURCE_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars)

file(GLOB COLVARS_SOURCES CONFIGURE_DEPENDS ${COLVARS_SOURCE_DIR}/[^.]*.cpp)

option(COLVARS_DEBUG "Enable debugging messages for Colvars (quite verbose)" OFF)

option(COLVARS_LEPTON "Use the Lepton library for custom expressions" ON)

if(COLVARS_LEPTON)
  if(NOT LEPTON_SOURCE_DIR)
    include(Packages/LEPTON)
  endif()
endif()

add_library(colvars STATIC ${COLVARS_SOURCES})
target_compile_definitions(colvars PRIVATE -DCOLVARS_LAMMPS)
separate_arguments(CMAKE_TUNE_FLAGS)
foreach(_FLAG ${CMAKE_TUNE_FLAGS})
  target_compile_options(colvars PRIVATE ${_FLAG})
endforeach()
set_target_properties(colvars PROPERTIES OUTPUT_NAME lammps_colvars${LAMMPS_MACHINE})
target_include_directories(colvars PUBLIC ${LAMMPS_LIB_SOURCE_DIR}/colvars)
# The line below is needed to locate math_eigen_impl.h
target_include_directories(colvars PRIVATE ${LAMMPS_SOURCE_DIR})
target_link_libraries(lammps PRIVATE colvars)

if(COLVARS_DEBUG)
  # Need to export the define publicly to be valid in interface code
  target_compile_definitions(colvars PUBLIC -DCOLVARS_DEBUG)
endif()

if(COLVARS_LEPTON)
  target_compile_definitions(colvars PRIVATE -DLEPTON)
  target_link_libraries(colvars PUBLIC lepton)
endif()
