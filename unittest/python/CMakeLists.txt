# Test calling Python from LAMMPS (and importing the LAMMPS module
# inside those functions). This can do an "anonymous" import of symbols
# from the executable, so the shared library is not needed. The
# availability of the PYTHON package is tested for inside the tester.

set(TEST_INPUT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})

# we must have shared libraries enabled for testing the python module
if(NOT BUILD_SHARED_LIBS)
  message(STATUS "Skipping Tests for the LAMMPS Python Module: must enable BUILD_SHARED_LIBS")
  return()
endif()

if(CMAKE_VERSION VERSION_LESS 3.12)
  find_package(PythonInterp 3.6) # Deprecated since version 3.12
  if(PYTHONINTERP_FOUND)
    set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
  endif()
else()
  find_package(Python COMPONENTS Interpreter Development)
endif()

add_executable(test_python_package test_python_package.cpp)
target_link_libraries(test_python_package PRIVATE lammps GTest::GMock)
target_compile_definitions(test_python_package PRIVATE -DTEST_INPUT_FOLDER=${TEST_INPUT_FOLDER})
# this requires CMake 3.12. don't care to add backward compatibility for this.
if(Python_Development_FOUND)
  target_compile_definitions(test_python_package PRIVATE -DTEST_HAVE_PYTHON_DEVELOPMENT=1)
  target_link_libraries(test_python_package PRIVATE Python::Python)
endif()
add_test(NAME PythonPackage COMMAND test_python_package)

# build list of environment variables for testing python functionality
if(WIN32)
  set(PYTHON_TEST_ENVIRONMENT PYTHONPATH=${LAMMPS_PYTHON_DIR})
else()
  set(PYTHON_TEST_ENVIRONMENT PYTHONPATH=${LAMMPS_PYTHON_DIR}:$ENV{PYTHONPATH})
endif()
get_property(BUILD_IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILD_IS_MULTI_CONFIG)
  set(LAMMPS_LIB_PATH ${CMAKE_BINARY_DIR}/$<CONFIG>)
else()
  set(LAMMPS_LIB_PATH ${CMAKE_BINARY_DIR})
endif()
list(APPEND PYTHON_TEST_ENVIRONMENT "LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS_DIR}")
list(APPEND PYTHON_TEST_ENVIRONMENT "PYTHONUNBUFFERED=1")
if(APPLE)
  list(APPEND PYTHON_TEST_ENVIRONMENT "DYLD_LIBRARY_PATH=${LAMMPS_LIB_PATH}:$ENV{DYLD_LIBRARY_PATH}")
elseif(WIN32)
  list(APPEND PYTHON_TEST_ENVIRONMENT "LAMMPSDLLPATH=${LAMMPS_LIB_PATH}")
else()
  list(APPEND PYTHON_TEST_ENVIRONMENT "LD_LIBRARY_PATH=${LAMMPS_LIB_PATH}:$ENV{LD_LIBRARY_PATH}")
endif()
set_tests_properties(PythonPackage PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")

if(Python_EXECUTABLE)
  list(APPEND PYTHON_TEST_ENVIRONMENT "TEST_INPUT_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
  list(APPEND PYTHON_TEST_ENVIRONMENT "LAMMPS_CMAKE_CACHE=${CMAKE_BINARY_DIR}/CMakeCache.txt")
  if(LAMMPS_MACHINE)
    # convert from '_machine' to 'machine'
    string(SUBSTRING ${LAMMPS_MACHINE} 1 -1 LAMMPS_MACHINE_NAME)
    list(APPEND PYTHON_TEST_ENVIRONMENT LAMMPS_MACHINE_NAME=${LAMMPS_MACHINE_NAME})
  endif()

  if(ENABLE_COVERAGE)
    find_program(COVERAGE_BINARY coverage)
    find_package_handle_standard_args(COVERAGE DEFAULT_MSG COVERAGE_BINARY)

    if(COVERAGE_FOUND)
        set(PYTHON_TEST_RUNNER ${Python_EXECUTABLE} -u ${COVERAGE_BINARY} run --rcfile=${CMAKE_BINARY_DIR}/.coveragerc)
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

  add_test(NAME PythonPizza
           COMMAND ${PYTHON_TEST_RUNNER} ${CMAKE_CURRENT_SOURCE_DIR}/python-pizza.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonPizza PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")

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
