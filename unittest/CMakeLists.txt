include(GTest)

# check if we can run the compiled executable and whether it prints
# the LAMMPS version header in the output for an empty input
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/in.empty "")
add_test(NAME RunLammps
         COMMAND $<TARGET_FILE:lmp> -log none -echo none -in in.empty
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(RunLammps PROPERTIES
	ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1;HWLOC_HIDE_ERRORS=1"
        PASS_REGULAR_EXPRESSION "^LAMMPS \\([0-9]+ [A-Za-z]+ 2[0-9][0-9][0-9]\\)")

# check if the compiled executable will print the help message
add_test(NAME HelpMessage
         COMMAND $<TARGET_FILE:lmp> -h
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(HelpMessage PROPERTIES
	ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1;HWLOC_HIDE_ERRORS=1"
         PASS_REGULAR_EXPRESSION ".*Large-scale Atomic/Molecular Massively Parallel Simulator -.*Usage example:.*")

# check if the compiled executable will error out on an invalid command line flag
add_test(NAME InvalidFlag
         COMMAND $<TARGET_FILE:lmp> -xxx
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(InvalidFlag PROPERTIES
	ENVIRONMENT "TSAN_OPTIONS=ignore_noninstrumented_modules=1;HWLOC_HIDE_ERRORS=1"
         PASS_REGULAR_EXPRESSION "ERROR: Invalid command-line argument.*")

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
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPI_TEST_NUM_PROCS} ${MPIEXEC_PREFLAGS}
                     ${EXECUTABLE} ${MPIEXEC_POSTFLAGS} ${ARGS}
    )
  endfunction()
endif()

add_subdirectory(utils)
add_subdirectory(formats)
add_subdirectory(commands)
add_subdirectory(c-library)
add_subdirectory(cplusplus)
add_subdirectory(fortran)
add_subdirectory(python)
add_subdirectory(tools)
add_subdirectory(force-styles)

find_package(ClangFormat 8.0)

if(ClangFormat_FOUND)
  add_custom_target(format-tests
    COMMAND ${ClangFormat_EXECUTABLE} --verbose -i -style=file */*.cpp */*.h
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
