# unit tests for the LAMMPS fortran module require a real MPI library
if(NOT BUILD_MPI)
  message(STATUS "Skipping Tests for the LAMMPS Fortran Module: BUILD_MPI=on is required")
  return()
endif()

include(CheckGeneratorSupport)
if(NOT CMAKE_GENERATOR_SUPPORT_FORTRAN)
  message(STATUS "Skipping Tests for the LAMMPS Fortran Module: no Fortran support in build tool")
  return()
endif()

include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
  find_package(MPI REQUIRED)
  get_filename_component(LAMMPS_FORTRAN_MODULE ${LAMMPS_SOURCE_DIR}/../examples/COUPLE/fortran/lammps.f90 ABSOLUTE)

  add_library(flammps STATIC ${LAMMPS_FORTRAN_MODULE})

  add_executable(fortran-create wrap-create.cpp fortran-create.f90)
  target_link_libraries(fortran-create PRIVATE flammps lammps "MPI::MPI_Fortran;GTest::GTest")
  add_test(FortranOpen fortran-create)

  add_executable(fortran-commands wrap-commands.cpp fortran-commands.f90)
  target_link_libraries(fortran-commands PRIVATE flammps lammps "MPI::MPI_Fortran;GTest::GTest")
  add_test(FortranCommands fortran-commands)
else()
  message(STATUS "Skipping Tests for the LAMMPS Fortran Module: no Fortran compiler")
endif()
