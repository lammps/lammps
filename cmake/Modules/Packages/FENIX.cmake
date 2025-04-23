find_package(fenix REQUIRED)

target_link_libraries(lammps PRIVATE fenix)
target_compile_definitions(lammps PRIVATE -DLAMMPS_FENIX)
target_link_libraries(lmp PRIVATE fenix)
target_compile_definitions(lmp PRIVATE -DLAMMPS_FENIX)
