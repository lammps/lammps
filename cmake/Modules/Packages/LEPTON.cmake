set(LEPTON_SOURCE_DIR ${LAMMPS_LIB_SOURCE_DIR}/lepton)

file(GLOB LEPTON_SOURCES ${LEPTON_SOURCE_DIR}/src/[^.]*.cpp)
add_library(lmplepton STATIC ${LEPTON_SOURCES})
target_compile_definitions(lmplepton PRIVATE -DLEPTON_BUILDING_STATIC_LIBRARY)
set_target_properties(lmplepton PROPERTIES OUTPUT_NAME lammps_lmplepton${LAMMPS_MACHINE})
target_include_directories(lmplepton PUBLIC ${LEPTON_SOURCE_DIR}/include)
target_link_libraries(lammps PRIVATE lmplepton)
