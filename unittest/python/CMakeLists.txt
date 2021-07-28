# Test calling Python from LAMMPS (and importing the LAMMPS module
# inside those functions). This can do an "anonymous" import of symbols
# from the executable, so the shared library is not needed. The
# availability of the PYTHON package is tested for inside the tester.

set(TEST_INPUT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(test_python_package test_python_package.cpp)
target_link_libraries(test_python_package PRIVATE lammps GTest::GMock GTest::GTest)
target_compile_definitions(test_python_package PRIVATE -DTEST_INPUT_FOLDER=${TEST_INPUT_FOLDER})
add_test(NAME PythonPackage COMMAND test_python_package WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(PythonPackage PROPERTIES ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR};PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${LAMMPS_PYTHON_DIR}:$ENV{PYTHONPATH};PYTHONUNBUFFERED=1")

# we must have shared libraries enabled for testing the python module
if(NOT BUILD_SHARED_LIBS)
  message(STATUS "Skipping Tests for the LAMMPS Python Module: must enable BUILD_SHARED_LIBS")
  return()
endif()

if(CMAKE_VERSION VERSION_LESS 3.12)
  find_package(PythonInterp 3.5) # Deprecated since version 3.12
  if(PYTHONINTERP_FOUND)
    set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
  endif()
else()
  find_package(Python3 COMPONENTS Interpreter)
endif()

if(Python_EXECUTABLE)
  # prepare to augment the environment so that the LAMMPS python module and the shared library is found.
  set(PYTHON_TEST_ENVIRONMENT PYTHONPATH=${LAMMPS_PYTHON_DIR}:$ENV{PYTHONPATH})
  list(APPEND PYTHON_TEST_ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
  list(APPEND PYTHON_TEST_ENVIRONMENT "TEST_INPUT_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
  if(APPLE)
    list(APPEND PYTHON_TEST_ENVIRONMENT "DYLD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$ENV{DYLD_LIBRARY_PATH};LAMMPS_CMAKE_CACHE=${CMAKE_BINARY_DIR}/CMakeCache.txt")
  else()
    list(APPEND PYTHON_TEST_ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$ENV{LD_LIBRARY_PATH};LAMMPS_CMAKE_CACHE=${CMAKE_BINARY_DIR}/CMakeCache.txt")
  endif()
  if(LAMMPS_MACHINE)
    # convert from '_machine' to 'machine'
    string(SUBSTRING ${LAMMPS_MACHINE} 1 -1 LAMMPS_MACHINE_NAME)
    list(APPEND PYTHON_TEST_ENVIRONMENT LAMMPS_MACHINE_NAME=${LAMMPS_MACHINE_NAME})
  endif()

  if(ENABLE_COVERAGE)
    find_program(COVERAGE_BINARY coverage)
    find_package_handle_standard_args(COVERAGE DEFAULT_MSG COVERAGE_BINARY)

    if(COVERAGE_FOUND)
        set(PYTHON_TEST_RUNNER ${Python_EXECUTABLE} -u ${COVERAGE_BINARY} run --parallel-mode --include=${LAMMPS_PYTHON_DIR}/lammps/*.py --omit=${LAMMPS_PYTHON_DIR}/install.py)
    else()
        set(PYTHON_TEST_RUNNER ${Python_EXECUTABLE} -u)
    endif()
  else()
    set(PYTHON_TEST_RUNNER ${Python_EXECUTABLE} -u)
  endif()

  add_test(NAME PythonOpen
           COMMAND ${PYTHON_TEST_RUNNER} ${CMAKE_CURRENT_SOURCE_DIR}/python-open.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonOpen PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")

  # some of the tests in this file will fail without exceptions enabled
  if(LAMMPS_EXCEPTIONS)
    add_test(NAME PythonCommands
             COMMAND ${PYTHON_TEST_RUNNER} ${CMAKE_CURRENT_SOURCE_DIR}/python-commands.py -v
             WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
    set_tests_properties(PythonCommands PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")
  endif()

  add_test(NAME PythonNumpy
           COMMAND ${PYTHON_TEST_RUNNER} ${CMAKE_CURRENT_SOURCE_DIR}/python-numpy.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonNumpy PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")

  add_test(NAME PythonCapabilities
           COMMAND ${PYTHON_TEST_RUNNER} ${CMAKE_CURRENT_SOURCE_DIR}/python-capabilities.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonCapabilities PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")

  add_test(NAME PythonPyLammps
           COMMAND ${PYTHON_TEST_RUNNER} ${CMAKE_CURRENT_SOURCE_DIR}/python-pylammps.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonPyLammps PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")

  add_test(NAME PythonFormats
           COMMAND ${PYTHON_TEST_RUNNER} ${CMAKE_CURRENT_SOURCE_DIR}/python-formats.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonFormats PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")

  add_test(NAME PythonFixExternal
           COMMAND ${PYTHON_TEST_RUNNER} ${CMAKE_CURRENT_SOURCE_DIR}/python-fix-external.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonFixExternal PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")

  add_test(NAME PythonScatterGather
           COMMAND ${PYTHON_TEST_RUNNER} ${CMAKE_CURRENT_SOURCE_DIR}/python-scatter-gather.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonScatterGather PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")

else()
  message(STATUS "Skipping Tests for the LAMMPS Python Module: no suitable Python interpreter")
endif()
