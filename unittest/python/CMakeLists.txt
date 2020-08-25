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

if (Python_EXECUTABLE)
  # prepare to augment the environment so that the LAMMPS python module and the shared library is found.
  set(PYTHON_TEST_ENVIRONMENT PYTHONPATH=${LAMMPS_PYTHON_DIR}:$ENV{PYTHONPATH})
  if(APPLE)
    list(APPEND PYTHON_TEST_ENVIRONMENT DYLD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$ENV{DYLD_LIBRARY_PATH})
  else()
    list(APPEND PYTHON_TEST_ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$ENV{LD_LIBRARY_PATH})
  endif()
  if(LAMMPS_MACHINE)
    # convert from '_machine' to 'machine'
    string(SUBSTRING ${LAMMPS_MACHINE} 1 -1 LAMMPS_MACHINE_NAME)
    list(APPEND PYTHON_TEST_ENVIRONMENT LAMMPS_MACHINE_NAME=${LAMMPS_MACHINE_NAME})
  endif()

  add_test(NAME PythonOpen
           COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python-open.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonOpen PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")
  add_test(NAME PythonCommands
           COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python-commands.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonCommands PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")
else()
  message(STATUS "Skipping Tests for the LAMMPS Python Module: no suitable Python interpreter")
endif()
