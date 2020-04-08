enable_language(Fortran)
find_package(QUIP REQUIRED)
target_link_libraries(lammps PRIVATE QUIP::QUIP ${LAPACK_LIBRARIES})
