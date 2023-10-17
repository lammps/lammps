
if(NOT Python_INTERPRETER)
  # backward compatibility
  if(PYTHON_EXECUTABLE)
    set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
  endif()
  find_package(Python COMPONENTS Interpreter)
endif()
find_package(Python REQUIRED COMPONENTS Interpreter Development)
target_link_libraries(lammps PRIVATE Python::Python)
target_compile_definitions(lammps PRIVATE -DLMP_PYTHON)
