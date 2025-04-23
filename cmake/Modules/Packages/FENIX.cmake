find_package(fenix REQUIRED)

#set(FENIX_PKG_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/FENIX)
#set(FENIX_PKG_SOURCES ${FENIX_PKG_SOURCES_DIR}/universe_fenix.cpp
#                      ${FENIX_PKG_SOURCES_DIR}/fail_command.cpp
#                      ${FENIX_PKG_SOURCES_DIR}/fix_fail.cpp)

target_link_libraries(lammps PRIVATE fenix)
target_compile_definitions(lammps PRIVATE -DLAMMPS_FENIX)
target_link_libraries(lmp PRIVATE fenix)
target_compile_definitions(lmp PRIVATE -DLAMMPS_FENIX)
