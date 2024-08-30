find_package(GSL 2.6 REQUIRED)
target_link_libraries(lammps PRIVATE GSL::gsl)
