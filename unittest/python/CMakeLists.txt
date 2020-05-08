# we must have shared libraries enabled for testing the python module
if(NOT BUILD_SHARED_LIBS)
  message(STATUS "Skipping Tests for the LAMMPS Python Module: must enable BUILD_SHARED_LIBS")
  return()
endif()

if(CMAKE_VERSION VERSION_LESS 3.12)
  find_package(PythonInterp) # Deprecated since version 3.12
  if(PYTHONINTERP_FOUND)
    set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
  endif()
else()
  find_package(Python COMPONENTS Interpreter)
endif()

if (Python_EXECUTABLE)
  add_test(NAME PythonOpen
           COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python-open.py -v
           WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
  set_tests_properties(PythonOpen PROPERTIES ENVIRONMENT
           PYTHONPATH=${LAMMPS_PYTHON_DIR}:$ENV{PYTHONPATH})
else()
  message(STATUS "Skipping Tests for the LAMMPS Python Module: no suitable Python interpreter")
endif()
