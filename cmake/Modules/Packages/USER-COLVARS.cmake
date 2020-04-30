set(COLVARS_SOURCE_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars)

file(GLOB COLVARS_SOURCES ${COLVARS_SOURCE_DIR}/[^.]*.cpp)

# Build Lepton by default
option(COLVARS_LEPTON "Build and link the Lepton library" ON)

if(COLVARS_LEPTON)
  set(LEPTON_DIR ${LAMMPS_LIB_SOURCE_DIR}/colvars/lepton)
  file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/[^.]*.cpp)
  add_library(lepton STATIC ${LEPTON_SOURCES})
  # Change the define below to LEPTON_BUILDING_SHARED_LIBRARY when linking Lepton as a DLL with MSVC
  target_compile_definitions(lepton PRIVATE -DLEPTON_BUILDING_STATIC_LIBRARY)
  set_target_properties(lepton PROPERTIES OUTPUT_NAME lammps_lepton${LAMMPS_MACHINE})
  target_include_directories(lepton PRIVATE ${LEPTON_DIR}/include)
endif()

add_library(colvars STATIC ${COLVARS_SOURCES})
target_compile_definitions(colvars PRIVATE -DLAMMPS_${LAMMPS_SIZES})
set_target_properties(colvars PROPERTIES OUTPUT_NAME lammps_colvars${LAMMPS_MACHINE})
target_include_directories(colvars PUBLIC ${LAMMPS_LIB_SOURCE_DIR}/colvars)
target_link_libraries(lammps PRIVATE colvars)

if(COLVARS_LEPTON)
  target_link_libraries(lammps PRIVATE lepton)
  target_compile_definitions(colvars PRIVATE -DLEPTON)
  # Disable the line below when linking Lepton as a DLL with MSVC
  target_compile_definitions(colvars PRIVATE -DLEPTON_USE_STATIC_LIBRARIES)
  target_include_directories(colvars PUBLIC ${LEPTON_DIR}/include)
endif()
