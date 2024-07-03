find_package(GSL 2.7 REQUIRED)
target_link_libraries(lammps PRIVATE GSL::gsl)
