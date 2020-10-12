include(GTest)

if(BUILD_MPI)
  function(add_mpi_test name binary nproc)
    add_test(NAME ${name} COMMAND ${MPIEXEC_EXECUTABLE} -np ${nproc} ${binary})
  endfunction()
endif()

add_subdirectory(utils)
add_subdirectory(formats)
add_subdirectory(commands)
add_subdirectory(c-library)
add_subdirectory(cplusplus)
add_subdirectory(fortran)
add_subdirectory(python)
add_subdirectory(force-styles)

find_package(ClangFormat 8.0)

if(ClangFormat_FOUND)
  set(UNITTEST_SOURCES)
  file(GLOB_RECURSE UNITTEST_SOURCES *.cpp *.h)
  add_custom_target(format-tests
    COMMAND ${ClangFormat_EXECUTABLE} --verbose -i -style=file ${UNITTEST_SOURCES}
    DEPENDS ${UNITTEST_SOURCES})
endif()
