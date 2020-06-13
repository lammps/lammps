if(CMAKE_VERSION VERSION_LESS 3.12)
    find_package(PythonInterp) # Deprecated since version 3.12
    if(PYTHONINTERP_FOUND)
        set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
    endif()
else()
    find_package(Python COMPONENTS Interpreter)
endif()

if (Python_EXECUTABLE)
    add_custom_target(
      check-whitespace
      ${Python_EXECUTABLE} ${LAMMPS_TOOLS_DIR}/coding_standard/whitespace.py .
      WORKING_DIRECTORY  ${LAMMPS_DIR}
      COMMENT "Check for whitespace errors")
    add_custom_target(
      fix-whitespace
      ${Python_EXECUTABLE} ${LAMMPS_TOOLS_DIR}/coding_standard/whitespace.py -f .
      WORKING_DIRECTORY  ${LAMMPS_DIR}
      COMMENT "Fix whitespace errors")
endif()
