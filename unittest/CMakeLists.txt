########################################
# CMake build for automated testing
# This file is part of LAMMPS
# Created by Axel Kohlmeyer and Richard Berger
########################################
# download and build googletest framework

# cannot compile googletest anymore with the default GCC on RHEL 7.x
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0))
  message(FATAL_ERROR "Need GNU C++ compiler version 6.x or later for unit testing")
endif()

message(STATUS "Downloading and building googletest framework")
set(GTEST_URL "https://github.com/google/googletest/archive/release-1.12.1.tar.gz" CACHE STRING "URL of googletest source")
set(GTEST_MD5 "e82199374acdfda3f425331028eb4e2a" CACHE STRING "MD5 sum for googletest source")
mark_as_advanced(GTEST_URL)
mark_as_advanced(GTEST_MD5)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

include(ExternalCMakeProject)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  set(cxx_strict_flags "-Wno-error -Wno-unused-command-line-argument")
endif()
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(INSTALL_GMOCK OFF CACHE BOOL "" FORCE)
ExternalCMakeProject(googletest ${GTEST_URL} ${GTEST_MD5} googletest . "")
add_library(GTest::GTest ALIAS gtest)
add_library(GTest::GMock ALIAS gmock)
add_library(GTest::GTestMain ALIAS gtest_main)
add_library(GTest::GMockMain ALIAS gmock_main)

# import
foreach(_FLAG ${CMAKE_TUNE_FLAGS})
  add_compile_options(${_FLAG})
endforeach()

# must repeat handling coverage for older CMake
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND ENABLE_COVERAGE)
  if(CMAKE_VERSION VERSION_LESS 3.13)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_${CMAKE_BUILD_TYPE}_FLAGS} --coverage")
  endif()
endif()

########################################
# General tests using the LAMMPS executable itself
########################################

# check if we can run the compiled executable and whether it prints
# the LAMMPS version header in the output for an empty input
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/in.empty "")
add_test(NAME RunLammps
         COMMAND $<TARGET_FILE:lmp> -log none -echo none -in in.empty)
set_tests_properties(RunLammps PROPERTIES
        ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1;HWLOC_HIDE_ERRORS=2"
        PASS_REGULAR_EXPRESSION "LAMMPS \\([0-9]+ [A-Za-z]+ 2[0-9][0-9][0-9]( - Update [0-9]+)?\\)")

# check if the compiled executable will print the help message
add_test(NAME HelpMessage
         COMMAND $<TARGET_FILE:lmp> -h)
set_tests_properties(HelpMessage PROPERTIES
        ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1;HWLOC_HIDE_ERRORS=2"
         PASS_REGULAR_EXPRESSION ".*Large-scale Atomic/Molecular Massively Parallel Simulator -.*Usage example:.*")

# check if the compiled executable will error out on an invalid command line flag
add_test(NAME InvalidFlag
         COMMAND $<TARGET_FILE:lmp> -xxx)
set_tests_properties(InvalidFlag PROPERTIES
        ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1;HWLOC_HIDE_ERRORS=2"
         PASS_REGULAR_EXPRESSION "ERROR: Invalid command-line argument.*")

# convenience function for adding tests requiring to be run in parallel with MPI
if(BUILD_MPI)
  function(add_mpi_test)
    set(MPI_TEST_NUM_PROCS 1)
    set(MPI_TEST_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    cmake_parse_arguments(MPI_TEST "" "NAME;NUM_PROCS;WORKING_DIRECTORY" "COMMAND" ${ARGN})
    # Do not add test when oversubscribing
    if(MPI_TEST_NUMPROCS GREATER MPIEXEC_MAX_NUMPROCS)
      return()
    endif()
    list(GET MPI_TEST_COMMAND 0 EXECUTABLE)
    list(REMOVE_AT MPI_TEST_COMMAND 0)
    set(ARGS ${MPI_TEST_COMMAND})
    add_test(NAME ${MPI_TEST_NAME}
             WORKING_DIRECTORY ${MPI_TEST_WORKING_DIRECTORY}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPI_TEST_NUM_PROCS}
                     ${MPIEXEC_PREFLAGS} ${EXECUTABLE} ${MPIEXEC_POSTFLAGS} ${ARGS})
  endfunction()
else()
  function(add_mpi_test)
    cmake_parse_arguments(MPI_TEST "" "NAME;NUM_PROCS;WORKING_DIRECTORY" "COMMAND" ${ARGN})
    message(STATUS "Skipping test ${NAME} on non-MPI compilation")
  endfunction()
endif()

# incorporate categories of specific tests from subdirectories

add_subdirectory(utils)
add_subdirectory(formats)
add_subdirectory(commands)
add_subdirectory(c-library)
add_subdirectory(cplusplus)
add_subdirectory(fortran)
add_subdirectory(python)
add_subdirectory(tools)
add_subdirectory(force-styles)

# clang-format support for test sources

find_package(ClangFormat 8.0)

if(ClangFormat_FOUND)
  add_custom_target(format-tests
    COMMAND ${ClangFormat_EXECUTABLE} --verbose -i -style=file */*.cpp */*.h
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
