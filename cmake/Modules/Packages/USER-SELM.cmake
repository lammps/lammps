find_package(FFTW3 REQUIRED)

set(SELM_SOURCE_DIR ${LAMMPS_LIB_SOURCE_DIR}/selm)
file(GLOB SELM_SOURCES ${SELM_SOURCE_DIR}/[^.]*.cpp)

get_target_property(LAMMPS_GLOBAL_DEFS lammps COMPILE_DEFINITIONS)
message(STATUS "LAMMPS defs: ${LAMMPS_GLOBAL_DEFS}")

add_library(selm STATIC ${SELM_SOURCES})
set_target_properties(selm PROPERTIES OUTPUT_NAME lammps_selm${LAMMPS_MACHINE})
target_compile_definitions(selm PRIVATE -DFFT_FFTW)
target_include_directories(selm PUBLIC ${LAMMPS_LIB_SOURCE_DIR}/selm)
target_include_directories(selm PRIVATE ${LAMMPS_SOURCE_DIR} ${LAMMPS_SOURCE_DIR}/USER-SELM)
target_link_libraries(selm PRIVATE MPI::MPI_CXX FFTW3::FFTW3)
target_link_libraries(lammps PRIVATE selm)
