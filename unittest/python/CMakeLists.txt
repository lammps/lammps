
if(CMAKE_VERSION VERSION_LESS 3.12)
  find_package(PythonInterp) # Deprecated since version 3.12
  if(PYTHONINTERP_FOUND)
    set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
  endif()
else()
  find_package(Python COMPONENTS Interpreter)
endif()

if (Python_EXECUTABLE AND BUILD_SHARED_LIBS)
   add_test(NAME PythonOpen
            COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python-open.py ${LAMMPS_PYTHON_DIR}
            WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
else()
  message(STATUS "Skipping Tests for the LAMMPS Python Module")
endif()
