include(CheckGeneratorSupport)
if(NOT CMAKE_GENERATOR_SUPPORT_FORTRAN)
  message(STATUS "Skipping Tests for the LAMMPS Fortran Module: no Fortran support in build tool")
  return()
endif()

include(CheckLanguage)
check_language(Fortran)

if(CMAKE_Fortran_COMPILER)
  enable_language(C)
  enable_language(Fortran)
  if(NOT CMAKE_Fortran_COMPILER_ID)
    message(STATUS "Skipping Tests for the LAMMPS Fortran Module: cannot identify Fortran compiler")
    return()
  endif()
  # GNU Fortran 4.8.x on RHEL/CentOS 7.x is not sufficient to compile the Fortran module.
  # Also GNU Fortran 7.x on Ubuntu18.04LTS is not sufficient. Work around flang being detected as GNU
  get_filename_component(_tmp_fc ${CMAKE_Fortran_COMPILER} NAME)
  if((CMAKE_Fortran_COMPILER_ID STREQUAL "GNU") AND (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9.0) AND NOT (_tmp_fc STREQUAL "flang"))
    message(STATUS "Need GNU Fortran compiler version 9.x or later for Fortran unit testing. Skipping...")
    return()
  endif()

  get_filename_component(LAMMPS_FORTRAN_MODULE ${LAMMPS_SOURCE_DIR}/../fortran/lammps.f90 ABSOLUTE)

  add_library(flammps STATIC ${LAMMPS_FORTRAN_MODULE} keepstuff.f90)
  get_filename_component(_tmp_fc ${CMAKE_Fortran_COMPILER} NAME)
  if (_tmp_fc STREQUAL "flang")
    target_link_libraries(flammps PUBLIC gfortran)
  endif()

  add_executable(test_fortran_create wrap_create.cpp test_fortran_create.f90)
  target_link_libraries(test_fortran_create PRIVATE flammps lammps GTest::GTestMain)
  target_include_directories(test_fortran_create PRIVATE "${LAMMPS_SOURCE_DIR}/../fortran")
  add_test(NAME FortranOpen COMMAND test_fortran_create)

  add_executable(test_fortran_commands wrap_commands.cpp test_fortran_commands.f90)
  target_link_libraries(test_fortran_commands PRIVATE flammps lammps GTest::GTestMain)
  add_test(NAME FortranCommands COMMAND test_fortran_commands)

  add_executable(test_fortran_get_thermo wrap_get_thermo.cpp test_fortran_get_thermo.f90)
  target_link_libraries(test_fortran_get_thermo PRIVATE flammps lammps GTest::GTestMain)
  add_test(NAME FortranGetThermo COMMAND test_fortran_get_thermo)

  add_executable(test_fortran_box wrap_box.cpp test_fortran_box.f90)
  target_link_libraries(test_fortran_box PRIVATE flammps lammps GTest::GTestMain)
  add_test(NAME FortranBox COMMAND test_fortran_box)

  add_executable(test_fortran_properties wrap_properties.cpp test_fortran_properties.f90 test_fortran_commands.f90)
  target_link_libraries(test_fortran_properties PRIVATE flammps lammps GTest::GMockMain)
  add_test(NAME FortranProperties COMMAND test_fortran_properties)

  add_executable(test_fortran_extract_global wrap_extract_global.cpp test_fortran_extract_global.f90)
  target_link_libraries(test_fortran_extract_global PRIVATE flammps lammps GTest::GTestMain)
  add_test(NAME FortranExtractGlobal COMMAND test_fortran_extract_global)

  add_executable(test_fortran_extract_atom wrap_extract_atom.cpp test_fortran_extract_atom.f90)
  target_link_libraries(test_fortran_extract_atom PRIVATE flammps lammps GTest::GTestMain)
  add_test(NAME FortranExtractAtom COMMAND test_fortran_extract_atom)

  add_executable(test_fortran_extract_compute wrap_extract_compute.cpp test_fortran_extract_compute.f90)
  target_link_libraries(test_fortran_extract_compute PRIVATE flammps lammps GTest::GTestMain)
  add_test(NAME FortranExtractCompute COMMAND test_fortran_extract_compute)

  add_executable(test_fortran_extract_fix wrap_extract_fix.cpp test_fortran_extract_fix.f90)
  target_link_libraries(test_fortran_extract_fix PRIVATE flammps lammps GTest::GTestMain)
  add_test(NAME FortranExtractFix COMMAND test_fortran_extract_fix)

  add_executable(test_fortran_extract_variable wrap_extract_variable.cpp test_fortran_extract_variable.f90)
  target_compile_definitions(test_fortran_extract_variable PRIVATE -DTEST_INPUT_FOLDER=${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(test_fortran_extract_variable PRIVATE flammps lammps GTest::GTestMain)
  add_test(NAME FortranExtractVariable COMMAND test_fortran_extract_variable)

  add_executable(test_fortran_gather_scatter wrap_gather_scatter.cpp test_fortran_gather_scatter.f90)
  target_link_libraries(test_fortran_gather_scatter PRIVATE flammps lammps GTest::GTestMain)
  add_test(NAME FortranGatherScatter COMMAND test_fortran_gather_scatter)

  add_executable(test_fortran_create_atoms wrap_create_atoms.cpp test_fortran_create_atoms.f90)
  target_link_libraries(test_fortran_create_atoms PRIVATE flammps lammps GTest::GTestMain)
  add_test(NAME FortranCreateAtoms COMMAND test_fortran_create_atoms)

  add_executable(test_fortran_configuration wrap_configuration.cpp test_fortran_configuration.f90 test_fortran_commands.f90)
  target_link_libraries(test_fortran_configuration PRIVATE flammps lammps GTest::GMockMain)
  add_test(NAME FortranConfiguration COMMAND test_fortran_configuration)

  add_executable(test_fortran_neighlist wrap_neighlist.cpp test_fortran_neighlist.f90)
  target_link_libraries(test_fortran_neighlist PRIVATE flammps lammps GTest::GMockMain)
  add_test(NAME FortranNeighlist COMMAND test_fortran_neighlist)

  # gfortran from GCC 9 on ubuntu has an ICE compiling this test.
  if ((CMAKE_Fortran_COMPILER_ID STREQUAL "GNU") AND (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 9.9))
    add_executable(test_fortran_fixexternal wrap_fixexternal.cpp test_fortran_fixexternal.f90)
    target_link_libraries(test_fortran_fixexternal PRIVATE flammps lammps GTest::GMockMain)
    add_test(NAME FortranFixExternal COMMAND test_fortran_fixexternal)
  endif()

else()
  message(STATUS "Skipping Tests for the LAMMPS Fortran Module: no Fortran compiler")
endif()
