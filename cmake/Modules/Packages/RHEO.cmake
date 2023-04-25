find_package(GSL REQUIRED)
target_link_libraries(lammps PRIVATE GSL::gsl)
