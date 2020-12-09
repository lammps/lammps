# if PYTHON package is included we may also include Python support in MLIAP
set(MLIAP_ENABLE_PYTHON_DEFAULT OFF)
if(PKG_PYTHON)
  find_package(Cythonize)
  if(Cythonize_FOUND)
    set(MLIAP_ENABLE_PYTHON_DEFAULT ON)
  endif()
endif()

option(MLIAP_ENABLE_PYTHON "Build MLIAP package with Python support" ${MLIAP_ENABLE_PYTHON_DEFAULT})

if(MLIAP_ENABLE_PYTHON)
  find_package(Cythonize REQUIRED)
  if(NOT PKG_PYTHON)
    message(FATAL_ERROR "Must enable PYTHON package for including Python support in MLIAP")
  endif()

  set(MLIAP_CYTHON_DIR ${CMAKE_BINARY_DIR}/cython)
  file(MAKE_DIRECTORY ${MLIAP_CYTHON_DIR})
  add_custom_command(OUTPUT  ${MLIAP_CYTHON_DIR}/mliap_model_python_couple.cpp ${MLIAP_CYTHON_DIR}/mliap_model_python_couple.h
            COMMAND            ${CMAKE_COMMAND} -E copy ${LAMMPS_SOURCE_DIR}/MLIAP/mliap_model_python_couple.pyx ${MLIAP_CYTHON_DIR}/mliap_model_python_couple.pyx
            COMMAND            ${Cythonize_EXECUTABLE} -3 ${MLIAP_CYTHON_DIR}/mliap_model_python_couple.pyx
            WORKING_DIRECTORY  ${MLIAP_CYTHON_DIR}
            MAIN_DEPENDENCY    ${LAMMPS_SOURCE_DIR}/MLIAP/mliap_model_python_couple.pyx
            COMMENT "Generating C++ sources with cythonize...")
  target_compile_definitions(lammps PRIVATE -DMLIAP_PYTHON)
  target_sources(lammps PRIVATE ${MLIAP_CYTHON_DIR}/mliap_model_python_couple.cpp)
  target_include_directories(lammps PRIVATE ${MLIAP_CYTHON_DIR})
endif()
