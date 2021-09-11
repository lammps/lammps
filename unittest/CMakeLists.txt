########################################
# CMake build for automated testing
# This file is part of LAMMPS
# Created by Axel Kohlmeyer and Richard Berger
########################################
# download and build googletest framework
message(STATUS "Downloading and building googletest framework")
set(GTEST_URL "https://github.com/google/googletest/archive/release-1.11.0.tar.gz" CACHE STRING "URL of googletest source")
set(GTEST_MD5 "e8a8df240b6938bb6384155d4c37d937" CACHE STRING "MD5 sum for googletest source")
mark_as_advanced(GTEST_URL)
mark_as_advanced(GTEST_MD5)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# workaround for older CMake versions (tested on Ubuntu 18.04 with CMake 3.10)
if(CMAKE_VERSION VERSION_LESS 3.14)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_deps)
  file(DOWNLOAD ${GTEST_URL} ${CMAKE_BINARY_DIR}/_deps/googletest.tar.gz EXPECTED_HASH MD5=${GTEST_MD5})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_BINARY_DIR}/_deps/googletest.tar.gz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/_deps)
  file(GLOB GTEST_SOURCE "${CMAKE_BINARY_DIR}/_deps/googletest-*")
  # sanity check. do not allow to have multiple downloaded and extracted versions of the source
  list(LENGTH GTEST_SOURCE _num)
  if(_num GREATER 1)
    message(FATAL_ERROR "Inconsistent googletest sources. Please delete ${CMAKE_BINARY_DIR}/_deps and re-run cmake")
  endif()
  add_subdirectory(${GTEST_SOURCE} ${CMAKE_BINARY_DIR}/googletest)
else()
  include(FetchContent)
  FetchContent_Declare(googletest URL ${GTEST_URL} URL_HASH MD5=${GTEST_MD5})
  FetchContent_MakeAvailable(googletest)
endif()

########################################
# General tests using the LAMMPS executable itself
########################################

# check if we can run the compiled executable and whether it prints
# the LAMMPS version header in the output for an empty input
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/in.empty "")
add_test(NAME RunLammps
         COMMAND $<TARGET_FILE:lmp> -log none -echo none -in in.empty
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(RunLammps PROPERTIES
        ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1"
        PASS_REGULAR_EXPRESSION "^LAMMPS \\([0-9]+ [A-Za-z]+ 2[0-9][0-9][0-9]\\)")

# check if the compiled executable will print the help message
add_test(NAME HelpMessage
         COMMAND $<TARGET_FILE:lmp> -h
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(HelpMessage PROPERTIES
         ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1"
         PASS_REGULAR_EXPRESSION ".*Large-scale Atomic/Molecular Massively Parallel Simulator -.*Usage example:.*")

# check if the compiled executable will error out on an invalid command line flag
add_test(NAME InvalidFlag
         COMMAND $<TARGET_FILE:lmp> -xxx
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(InvalidFlag PROPERTIES
         ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1"
         PASS_REGULAR_EXPRESSION "ERROR: Invalid command-line argument.*")

# convenience function for adding tests requiring to be run in parallel with MPI
if(BUILD_MPI)
  function(add_mpi_test)
    set(MPI_TEST_NUM_PROCS 1)
    set(MPI_TEST_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    cmake_parse_arguments(MPI_TEST "" "NAME;NUM_PROCS;WORKING_DIRECTORY" "COMMAND" ${ARGN})
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
