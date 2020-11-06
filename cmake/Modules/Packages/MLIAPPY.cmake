if(CMAKE_VERSION VERSION_LESS 3.12)
  #This block was not tested, mimmicks PYTHON.cmake.
  find_package(PythonLibs REQUIRED) # Deprecated since version 3.12
  target_include_directories(lammps PRIVATE ${PYTHON_INCLUDE_DIR})
  target_link_libraries(lammps PRIVATE ${PYTHON_LIBRARY})
  target_include_directories(lammps PRIVATE ${Python_NumPy_INCLUDE_DIRS})
else()
  find_package(Python REQUIRED COMPONENTS NumPy)
  target_include_directories(lammps PRIVATE ${Python_NumPy_INCLUDE_DIRS})
endif()
target_compile_definitions(lammps PRIVATE -DLMP_MLIAPPY)

