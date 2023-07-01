# -*- CMake -*- file for tests of utily functions and classes in LAMMPS

# we use python 3's subprocess module to run the tools and check the output
find_package(Python 3.6 COMPONENTS Interpreter)

get_property(BUILD_IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILD_IS_MULTI_CONFIG)
  set(LAMMPS_SHELL_EXE_DIR ${CMAKE_BINARY_DIR}/$<CONFIG>)
else()
  set(LAMMPS_SHELL_EXE_DIR ${CMAKE_BINARY_DIR})
endif()

if(Python_EXECUTABLE)
  if(BUILD_LAMMPS_SHELL)
    add_test(NAME LammpsShell
      COMMAND ${Python_EXECUTABLE} -u ${CMAKE_CURRENT_SOURCE_DIR}/test_lammps_shell.py -v
      WORKING_DIRECTORY ${LAMMPS_SHELL_EXE_DIR})
  endif()
else()
  message(STATUS "Skipping Tests for LAMMPS tools: no suitable Python interpreter")
endif()
