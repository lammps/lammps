# if PYTHON package is included we may also include Python support in ML-IAP
set(MLIAP_ENABLE_PYTHON_DEFAULT OFF)
if(PKG_PYTHON)
  find_package(Cythonize QUIET)
  if(Cythonize_FOUND)
    set(MLIAP_ENABLE_PYTHON_DEFAULT ON)
  endif()
endif()

option(MLIAP_ENABLE_PYTHON "Build ML-IAP package with Python support" ${MLIAP_ENABLE_PYTHON_DEFAULT})

if(MLIAP_ENABLE_PYTHON)
  find_package(Cythonize REQUIRED)
  if(NOT PKG_PYTHON)
    message(FATAL_ERROR "Must enable PYTHON package for including Python support in ML-IAP")
  endif()
  if(CMAKE_VERSION VERSION_LESS 3.12)
    if(PYTHONLIBS_VERSION_STRING VERSION_LESS 3.6)
      message(FATAL_ERROR "Python support in ML-IAP requires Python 3.6 or later")
    endif()
  else()
    if(Python_VERSION VERSION_LESS 3.6)
      message(FATAL_ERROR "Python support in ML-IAP requires Python 3.6 or later")
    endif()
  endif()

  set(MLIAP_BINARY_DIR ${CMAKE_BINARY_DIR}/cython)
  set(MLIAP_CYTHON_SRC ${LAMMPS_SOURCE_DIR}/ML-IAP/mliap_model_python_couple.pyx)
  get_filename_component(MLIAP_CYTHON_BASE ${MLIAP_CYTHON_SRC} NAME_WE)
  file(MAKE_DIRECTORY ${MLIAP_BINARY_DIR})
  add_custom_command(OUTPUT  ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.cpp ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.h
          COMMAND            ${CMAKE_COMMAND} -E copy_if_different ${MLIAP_CYTHON_SRC} ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.pyx
          COMMAND            ${Cythonize_EXECUTABLE} -3 ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.pyx
          WORKING_DIRECTORY  ${MLIAP_BINARY_DIR}
          MAIN_DEPENDENCY    ${MLIAP_CYTHON_SRC}
          COMMENT "Generating C++ sources with cythonize...")
  target_compile_definitions(lammps PRIVATE -DMLIAP_PYTHON)
  target_sources(lammps PRIVATE ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.cpp)
  target_include_directories(lammps PRIVATE ${MLIAP_BINARY_DIR})
endif()
